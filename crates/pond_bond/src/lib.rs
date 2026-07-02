//! # pond_bond — Bond-based peridynamics
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
//! the [`PdAtom::n0`](pond_atom::PdAtom) column. Both live in registered
//! `AtomData`, so they are permuted and migrated *with* the atoms by the
//! substrate — no side resource to fall out of sync when SOIL spatially sorts
//! or migrates points.
//!
//! ## Status
//!
//! Family construction ([`build_families`]) is an `O(N²)` reference-config pass
//! over **local** points, so the current slice is single-process
//! (`--no-default-features`). Cross-rank bonds and a neighbour-list-accelerated
//! build are the next step.

use std::collections::HashMap;

use grass_app::prelude::*;
use grass_scheduler::prelude::*;

use pond_atom::{PdAtom, PdMaterialTable};
use soil_core::{
    Atom, AtomDataRegistry, BondEntry, BondStore, CommResource, Domain, ParticleSimScheduleSet,
    ScheduleSetupSet,
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

        app.add_setup_system(
            build_families.label("pond_build_families"),
            ScheduleSetupSet::PostSetup,
        );
        app.add_update_system(bond_based_force, ParticleSimScheduleSet::Force);
    }
}

/// Setup system: forms the peridynamic family in the reference configuration.
///
/// For every local point pair within the horizon, pushes a [`BondEntry`] onto
/// both endpoints (so each point's family is complete) recording the reference
/// bond length `r₀`, then records each point's initial bond count in
/// [`PdAtom::n0`](pond_atom::PdAtom).
pub fn build_families(
    atoms: Res<Atom>,
    registry: Res<AtomDataRegistry>,
    table: Res<PdMaterialTable>,
    domain: Res<Domain>,
    comm: Res<CommResource>,
) {
    let mut bond_store = registry.expect_mut::<BondStore>("build_families");
    let mut pd = registry.expect_mut::<PdAtom>("build_families");
    let nlocal = atoms.nlocal as usize;

    // Idempotent: only build once (setup may be revisited across stages).
    if bond_store.bonds.iter().take(nlocal).any(|b| !b.is_empty()) {
        return;
    }
    while bond_store.bonds.len() < nlocal {
        bond_store.bonds.push(Vec::new());
    }

    let horizon = table.horizon;
    let pflags = domain.periodic_flags();
    let box_size = domain.size;
    let mut nbonds = 0u64;

    for i in 0..nlocal {
        for j in (i + 1)..nlocal {
            let mut d = [
                atoms.pos[j][0] - atoms.pos[i][0],
                atoms.pos[j][1] - atoms.pos[i][1],
                atoms.pos[j][2] - atoms.pos[i][2],
            ];
            // Minimum image on periodic axes (reference configuration).
            for k in 0..3 {
                if pflags[k] {
                    d[k] -= box_size[k] * (d[k] / box_size[k]).round();
                }
            }
            let dist = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
            if dist > 0.0 && dist <= horizon {
                bond_store.bonds[i].push(BondEntry {
                    partner_tag: atoms.tag[j],
                    bond_type: BOND_INTACT,
                    r0: dist,
                });
                bond_store.bonds[j].push(BondEntry {
                    partner_tag: atoms.tag[i],
                    bond_type: BOND_INTACT,
                    r0: dist,
                });
                nbonds += 1;
            }
        }
    }

    // Record each point's initial bond count as the damage denominator.
    for i in 0..nlocal {
        pd.n0[i] = bond_store.bonds[i].len() as f64;
    }

    let total = comm.all_reduce_sum_f64(nbonds as f64) as u64;
    if comm.rank() == 0 {
        println!("BondBased: formed {total} peridynamic bonds (reference configuration)");
    }
}

/// Force system: accumulates the bond-based PD force on every local point and
/// updates its damage, breaking bonds that exceed the critical stretch.
pub fn bond_based_force(
    mut atoms: ResMut<Atom>,
    registry: Res<AtomDataRegistry>,
    table: Res<PdMaterialTable>,
) {
    let mut bond_store = registry.expect_mut::<BondStore>("bond_based_force");
    let mut pd = registry.expect_mut::<PdAtom>("bond_based_force");

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
