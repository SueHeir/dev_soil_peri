//! Cross-rank family construction consistency check.
//!
//! This small validation case builds one reference peridynamic family, carves a
//! deterministic crack-plane notch by breaking reference bonds, and reports
//! rank-reduced topology/damage/path metrics. The companion `sweep.py` runs the
//! same declarative lattice on one and two MPI ranks and compares the metrics.

use std::collections::HashMap;

use peri_core::prelude::*;

const Y_MID: f64 = 0.25;
const NOTCH_X: f64 = 0.30;
const DAMAGE_THRESHOLD: f64 = 0.0;

#[derive(Default)]
struct FamilyAuditState {
    done: bool,
}

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(BondBasedDefaultPlugins);
    app.add_resource(FamilyAuditState::default());
    app.add_update_system(
        audit_family
            .after("peri_build_families")
            .label("peri_family_consistency_audit"),
        ParticleSimScheduleSet::Neighbor,
    );
    app.start();
}

fn audit_family(
    atoms: Res<Atom>,
    registry: Res<AtomDataRegistry>,
    comm: Res<CommResource>,
    mut state: ResMut<FamilyAuditState>,
) {
    if state.done {
        return;
    }
    state.done = true;

    let mut bond_store = registry.expect_mut::<BondStore>("audit_family");
    let mut pd = registry.expect_mut::<PeriAtom>("audit_family");
    let nlocal = atoms.nlocal() as usize;
    let ntotal = atoms.len();

    let mut tag_to_index: HashMap<u32, usize> = HashMap::with_capacity(ntotal);
    for i in 0..ntotal {
        tag_to_index.insert(atoms.tag[i], i);
    }

    let mut local_entries = 0.0;
    let mut local_checksum = 0.0;
    let mut local_broken = 0.0;

    for i in 0..nlocal {
        let xi = pd.ref_pos[i];
        let n_bonds = bond_store.bonds[i].len();
        local_entries += n_bonds as f64;
        for k in 0..n_bonds {
            let partner_tag = bond_store.bonds[i][k].partner_tag;
            local_checksum += (atoms.tag[i] as f64) * 1_000_003.0 + (partner_tag as f64) * 97.0;

            let Some(&j) = tag_to_index.get(&partner_tag) else {
                continue;
            };
            let xj = pd.ref_pos[j];
            if (xi[1] - Y_MID) * (xj[1] - Y_MID) >= 0.0 {
                continue;
            }
            let t = (Y_MID - xi[1]) / (xj[1] - xi[1]);
            let x_cross = xi[0] + t * (xj[0] - xi[0]);
            if x_cross <= NOTCH_X && bond_store.bonds[i][k].bond_type != BOND_BROKEN {
                bond_store.bonds[i][k].bond_type = BOND_BROKEN;
                local_broken += 1.0;
            }
        }
    }

    let mut local_damage_sum = 0.0;
    let mut local_crack_count = 0.0;
    let mut local_crack_min = f64::INFINITY;
    let mut local_crack_max = f64::NEG_INFINITY;
    for i in 0..nlocal {
        let n0 = pd.n0[i];
        let nb = bond_store.bonds[i]
            .iter()
            .filter(|b| b.bond_type == BOND_BROKEN)
            .count() as f64;
        pd.damage[i] = if n0 > 0.0 { nb / n0 } else { 0.0 };
        local_damage_sum += pd.damage[i];
        if pd.damage[i] > DAMAGE_THRESHOLD {
            local_crack_count += 1.0;
            local_crack_min = local_crack_min.min(pd.ref_pos[i][0]);
            local_crack_max = local_crack_max.max(pd.ref_pos[i][0]);
        }
    }

    let entries = comm.all_reduce_sum_f64(local_entries);
    let checksum = comm.all_reduce_sum_f64(local_checksum);
    let broken = comm.all_reduce_sum_f64(local_broken);
    let damage_sum = comm.all_reduce_sum_f64(local_damage_sum);
    let crack_count = comm.all_reduce_sum_f64(local_crack_count);
    let crack_min = comm.all_reduce_min_f64(local_crack_min);
    let crack_max = -comm.all_reduce_min_f64(-local_crack_max);
    let crack_span = if crack_count > 0.0 {
        crack_max - crack_min
    } else {
        0.0
    };

    if comm.rank() == 0 {
        println!(
            "FAMILYCHECK ranks={} entries={:.0} checksum={:.6} broken={:.0} damage_sum={:.12} crack_count={:.0} crack_span={:.12}",
            comm.size(),
            entries,
            checksum,
            broken,
            damage_sum,
            crack_count,
            crack_span
        );
    }
}
