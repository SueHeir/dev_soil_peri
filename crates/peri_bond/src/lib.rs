//! # peri_bond — Bond-based peridynamics
//!
//! Bond-based peridynamics (PD) on the SOIL substrate. Each material point
//! bonds, once and for all in the **reference configuration**, to every other
//! point inside the horizon δ. Those bonds are stored in SOIL's
//! [`BondStore`](soil_core::BondStore) (partner addressed by global tag, so the
//! family survives migration and ghost exchange). Every step, each point sums a
//! pairwise force over its surviving bonds; a bond whose stretch exceeds the
//! critical stretch breaks irreversibly, and the fraction of broken bonds is the
//! point's **damage**.
//!
//! ## Force model
//!
//! For a bond with reference length `|ξ| = r₀` and current length `|y|`, the
//! stretch is `s = (|y| − r₀) / r₀`. The pairwise force *density* is
//! `f = c · s · ŷ`, with constant micromodulus `c = 18K/(πδ⁴)`. Discretised
//! over point volumes, the force on point `i` is
//!
//! ```text
//! F_i = V_i · Σ_j  c · s · (y_j − y_i)/|y_j − y_i| · V_j
//! ```
//!
//! applied to point `i` only — every point independently accumulates the force
//! *on itself* from each of its bonds. This avoids any Newton's-third-law /
//! reverse-communication bookkeeping and is embarrassingly parallel: a bond
//! `i–j` is evaluated twice (once by `i`, once by `j`), each endpoint using the
//! same ghost-synchronised geometry, so breakage stays consistent without
//! communication.
//!
//! ## Breakage state
//!
//! A broken bond is flagged in place by setting its [`BondEntry::bond_type`] to
//! [`BOND_BROKEN`]; the initial bond count (the damage denominator) is stored as
//! the [`PeriAtom::n0`](peri_atom::PeriAtom) column. Both live in registered
//! `AtomData`, so they are permuted and migrated *with* the atoms by the
//! substrate — no side resource to fall out of sync when SOIL spatially sorts
//! or migrates points.
//!
//! ## Status
//!
//! Family construction ([`build_families`]) walks SOIL's binned neighbour/halo
//! data in the reference configuration. Partners are recorded by global tag, so
//! bonds can cross MPI rank boundaries while the [`BondStore`] still migrates
//! and restarts with the owning atom.

use std::collections::{HashMap, HashSet};

use grass_app::prelude::*;
use grass_scheduler::prelude::*;

use peri_atom::{PeriAtom, PeriMaterialTable};
use soil_core::{
    Atom, AtomDataRegistry, BondEntry, BondStore, CommResource, Domain, Neighbor,
    ParticleSimScheduleSet,
};

/// `bond_type` value marking an intact bond.
pub const BOND_INTACT: u32 = 0;
/// `bond_type` value marking a broken bond (set in place when the bond fails).
pub const BOND_BROKEN: u32 = 1;

/// Bond-based peridynamics: builds the reference-config family and integrates
/// the constant-micromodulus bond force with critical-stretch breakage.
pub struct BondBasedPlugin;

impl Plugin for BondBasedPlugin {
    fn requires(&self) -> Vec<&str> {
        vec!["pd_points"]
    }

    fn provides(&self) -> Vec<&str> {
        vec!["pd_bond_forces"]
    }

    fn build(&self, app: &mut App) {
        // Registers BondStore into the AtomDataRegistry.
        app.add_plugins(soil_core::BondPlugin);

        app.add_update_system(
            build_families.label("peri_build_families"),
            ParticleSimScheduleSet::Neighbor,
        );
        app.add_update_system(bond_based_force, ParticleSimScheduleSet::Force);
    }
}

/// Setup system: forms the peridynamic family in the reference configuration.
///
/// For every local point, scans SOIL's binned neighbour/halo data for points
/// within the horizon, pushes a [`BondEntry`] addressed by global partner tag,
/// and records the point's initial bond count in
/// [`PeriAtom::n0`](peri_atom::PeriAtom). Cross-rank partners are stored only in
/// the local endpoint's list; the partner rank builds the reciprocal endpoint
/// from its own halo.
pub fn build_families(
    atoms: Res<Atom>,
    registry: Res<AtomDataRegistry>,
    table: Res<PeriMaterialTable>,
    domain: Res<Domain>,
    neighbor: Res<Neighbor>,
    comm: Res<CommResource>,
) {
    let mut bond_store = registry.expect_mut::<BondStore>("build_families");
    let mut pd = registry.expect_mut::<PeriAtom>("build_families");
    let nlocal = atoms.nlocal as usize;

    // Idempotent: the reference family is formed once, after SOIL has built the
    // initial halo and bin data.
    if pd.n0.iter().take(nlocal).any(|&n| n > 0.0) {
        return;
    }
    if neighbor.bin_start.is_empty() || neighbor.bin_sorted_atoms.is_empty() {
        return;
    }

    let nbonds = construct_families_from_neighbor_bins(
        &atoms,
        &mut bond_store,
        &mut pd,
        &neighbor,
        &domain,
        table.horizon,
    );

    let total = comm.all_reduce_sum_f64(nbonds as f64) as u64;
    if comm.rank() == 0 {
        println!("BondBased: formed {total} peridynamic family entries (reference configuration)");
    }
}

fn construct_families_from_neighbor_bins(
    atoms: &Atom,
    bond_store: &mut BondStore,
    pd: &mut PeriAtom,
    neighbor: &Neighbor,
    domain: &Domain,
    horizon: f64,
) -> u64 {
    let nlocal = atoms.nlocal as usize;
    let total = atoms.len();
    bond_store.bonds.clear();
    bond_store.bonds.resize_with(nlocal, Vec::new);

    let horizon_sq = horizon * horizon;
    let pflags = domain.periodic_flags();
    let box_size = domain.size;
    let mut family_entries = 0u64;

    for i in 0..nlocal {
        let mut seen: HashSet<u32> = HashSet::new();
        let my_cell = match neighbor.bin_atom_cell.get(i) {
            Some(&cell) => cell as usize,
            None => continue,
        };

        scan_family_cell(
            i,
            my_cell,
            atoms,
            bond_store,
            &mut seen,
            neighbor,
            horizon_sq,
            pflags,
            box_size,
            &mut family_entries,
        );

        for &offset in &neighbor.bin_stencil {
            if offset == 0 {
                continue;
            }
            let Some(cell) = (my_cell as i32).checked_add(offset) else {
                continue;
            };
            if cell < 0 {
                continue;
            }
            let cell = cell as usize;
            if cell + 1 >= neighbor.bin_start.len() {
                continue;
            }
            scan_family_cell(
                i,
                cell,
                atoms,
                bond_store,
                &mut seen,
                neighbor,
                horizon_sq,
                pflags,
                box_size,
                &mut family_entries,
            );
        }

        pd.n0[i] = bond_store.bonds[i].len() as f64;
    }

    debug_assert!(bond_store.bonds.len() >= nlocal);
    debug_assert!(total >= nlocal);
    family_entries
}

#[allow(clippy::too_many_arguments)]
fn scan_family_cell(
    i: usize,
    cell: usize,
    atoms: &Atom,
    bond_store: &mut BondStore,
    seen: &mut HashSet<u32>,
    neighbor: &Neighbor,
    horizon_sq: f64,
    pflags: [bool; 3],
    box_size: [f64; 3],
    family_entries: &mut u64,
) {
    let Some(&start) = neighbor.bin_start.get(cell) else {
        return;
    };
    let Some(&end) = neighbor.bin_start.get(cell + 1) else {
        return;
    };

    let xi = atoms.pos[i];
    for m in start as usize..end as usize {
        let Some(&j_u32) = neighbor.bin_sorted_atoms.get(m) else {
            continue;
        };
        let j = j_u32 as usize;
        if j == i || j >= atoms.len() {
            continue;
        }
        let tag = atoms.tag[j];
        if !seen.insert(tag) {
            continue;
        }

        let mut d = [
            atoms.pos[j][0] as f64 - xi[0] as f64,
            atoms.pos[j][1] as f64 - xi[1] as f64,
            atoms.pos[j][2] as f64 - xi[2] as f64,
        ];
        for k in 0..3 {
            if pflags[k] && box_size[k] > 0.0 {
                d[k] -= box_size[k] * (d[k] / box_size[k]).round();
            }
        }
        let r2 = d[0].mul_add(d[0], d[1].mul_add(d[1], d[2] * d[2]));
        if r2 > 0.0 && r2 <= horizon_sq {
            bond_store.bonds[i].push(BondEntry {
                partner_tag: tag,
                bond_type: BOND_INTACT,
                r0: r2.sqrt(),
            });
            *family_entries += 1;
        }
    }
}

/// Force system: accumulates the bond-based PD force on every local point and
/// updates its damage, breaking bonds that exceed the critical stretch.
pub fn bond_based_force(
    mut atoms: ResMut<Atom>,
    registry: Res<AtomDataRegistry>,
    table: Res<PeriMaterialTable>,
) {
    let mut bond_store = registry.expect_mut::<BondStore>("bond_based_force");
    let mut pd = registry.expect_mut::<PeriAtom>("bond_based_force");

    let nlocal = atoms.nlocal as usize;
    let ntotal = atoms.len();
    if bond_store.bonds.len() < nlocal {
        return;
    }

    // tag → current index (locals + ghosts), so a bond partner that is currently
    // a ghost still resolves to its position and volume.
    let mut tag_to_index: HashMap<u32, usize> = HashMap::with_capacity(ntotal);
    for idx in 0..ntotal {
        tag_to_index.insert(atoms.tag[idx], idx);
    }

    for i in 0..nlocal {
        let n_bonds = bond_store.bonds[i].len();
        if n_bonds == 0 {
            continue;
        }
        let vi = pd.volume[i];
        let ti = atoms.atom_type[i] as usize;
        let xi = atoms.pos[i];

        let mut f = [0.0f64; 3];
        let mut nbroken = 0.0f64;

        for k in 0..n_bonds {
            // Read the bond fields by value (BondEntry is not Copy); the
            // mutable borrow for breakage happens later via direct indexing.
            let bond = &bond_store.bonds[i][k];
            if bond.bond_type == BOND_BROKEN {
                nbroken += 1.0;
                continue;
            }
            let partner_tag = bond.partner_tag;
            let r0 = bond.r0;

            let j = match tag_to_index.get(&partner_tag) {
                Some(&j) => j,
                None => continue, // partner not currently present
            };

            let dx = atoms.pos[j][0] - xi[0];
            let dy = atoms.pos[j][1] - xi[1];
            let dz = atoms.pos[j][2] - xi[2];
            let ylen = (dx * dx + dy * dy + dz * dz).sqrt();
            if ylen == 0.0 {
                continue;
            }
            let s = (ylen - r0) / r0;

            let tj = atoms.atom_type[j] as usize;
            if s > table.critical_stretch_ij(ti, tj) {
                bond_store.bonds[i][k].bond_type = BOND_BROKEN;
                nbroken += 1.0;
                continue;
            }

            // F_i contribution = c · s · ŷ · V_i · V_j.
            let c = table.micromodulus_ij(ti, tj);
            let vj = pd.volume[j];
            let coeff = c * s * vi * vj / ylen;
            f[0] += coeff * dx;
            f[1] += coeff * dy;
            f[2] += coeff * dz;
        }

        atoms.force[i][0] += f[0];
        atoms.force[i][1] += f[1];
        atoms.force[i][2] += f[2];

        let n0 = pd.n0[i];
        pd.damage[i] = if n0 > 0.0 { nbroken / n0 } else { 0.0 };
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use soil_core::BoundaryType;

    fn atom_with_points(points: &[(u32, [f64; 3], bool)]) -> Atom {
        let mut atoms = Atom::new();
        atoms.nlocal = points.iter().filter(|(_, _, ghost)| !*ghost).count() as u32;
        atoms.nghost = points.iter().filter(|(_, _, ghost)| *ghost).count() as u32;
        atoms.natoms = atoms.nlocal as u64;
        for &(tag, pos, ghost) in points {
            atoms.tag.push(tag);
            atoms.atom_type.push(0);
            atoms.origin_index.push(if ghost { 0 } else { -1 });
            atoms.is_ghost.push(ghost);
            atoms.pos.push(pos);
            atoms.vel.push([0.0; 3]);
            atoms.force.push([0.0; 3]);
            atoms.cutoff_radius.push(0.6);
            atoms.mass.push(1.0);
            atoms.inv_mass.push(1.0);
            atoms.image.push([0, 0, 0]);
        }
        atoms
    }

    fn single_cell_neighbor(n: usize) -> Neighbor {
        let mut neighbor = Neighbor::new();
        neighbor.bin_atom_cell = vec![0; n];
        neighbor.bin_start = vec![0, n as u32];
        neighbor.bin_sorted_atoms = (0..n as u32).collect();
        neighbor.bin_stencil.clear();
        neighbor
    }

    fn fixed_domain() -> Domain {
        let mut domain = Domain::new();
        domain.boundary_type = [BoundaryType::Fixed; 3];
        domain
    }

    #[test]
    fn neighbor_backed_family_includes_each_local_row_and_ghost_tags() {
        let atoms = atom_with_points(&[
            (10, [0.0, 0.0, 0.0], false),
            (20, [0.5, 0.0, 0.0], false),
            (30, [0.9, 0.0, 0.0], true),
        ]);
        let neighbor = single_cell_neighbor(atoms.len());
        let domain = fixed_domain();
        let mut store = BondStore::new();
        let mut pd = PeriAtom {
            ref_pos: atoms.pos.clone(),
            volume: vec![1.0; atoms.len()],
            n0: vec![0.0; atoms.len()],
            damage: vec![0.0; atoms.len()],
        };

        let entries = construct_families_from_neighbor_bins(
            &atoms, &mut store, &mut pd, &neighbor, &domain, 0.6,
        );

        assert_eq!(entries, 3);
        assert_eq!(pd.n0[..2], [1.0, 2.0]);
        assert!(store.has_bond(0, 20));
        assert!(store.has_bond(1, 10));
        assert!(store.has_bond(1, 30));
    }

    #[test]
    fn neighbor_backed_family_filters_outside_horizon() {
        let atoms = atom_with_points(&[
            (10, [0.0, 0.0, 0.0], false),
            (20, [0.5, 0.0, 0.0], false),
            (30, [0.7, 0.0, 0.0], false),
        ]);
        let neighbor = single_cell_neighbor(atoms.len());
        let domain = fixed_domain();
        let mut store = BondStore::new();
        let mut pd = PeriAtom {
            ref_pos: atoms.pos.clone(),
            volume: vec![1.0; atoms.len()],
            n0: vec![0.0; atoms.len()],
            damage: vec![0.0; atoms.len()],
        };

        let entries = construct_families_from_neighbor_bins(
            &atoms, &mut store, &mut pd, &neighbor, &domain, 0.6,
        );

        assert_eq!(entries, 4);
        assert_eq!(pd.n0[..3], [1.0, 2.0, 1.0]);
        assert!(!store.has_bond(0, 30));
    }
}
