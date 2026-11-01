//! Bond-based peridynamic bar in uniaxial tension.
//!
//! A 3-D bar is gripped at both ends: the left grip is held fixed, the right
//! grip is pulled at a constant velocity. The nonlocal bond force carries the
//! resulting stress wave down the bar; bonds that exceed the critical stretch
//! break, and the damage field records where.
//!
//! ```bash
//! cargo run --release --example bb_bar_tension --no-default-features -- examples/bb_bar_tension/config.toml
//! ```

use peri_core::prelude::*;

/// Reference x-bounds of the bar [m] (must match `[[dev_soil_peri.lattice]]`).
const BAR_X_MIN: f64 = 0.0;
const BAR_X_MAX: f64 = 0.04;
/// Grip width at each end [m].
const GRIP: f64 = 0.004;
/// Constant pull velocity applied to the right grip [m/s].
const PULL_VEL: f64 = 2.0;

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(BondBasedDefaultPlugins);

    // Prescribe grip kinematics after the bond force is computed and before the
    // final velocity update, so the grips move at exactly the imposed velocity.
    app.add_update_system(apply_grips, ParticleSimScheduleSet::PostForce);

    // Periodic damage / energy diagnostic.
    app.add_update_system(report, ParticleSimScheduleSet::PostFinalIntegration);

    app.start();
}

/// Hold the left grip fixed and pull the right grip at constant velocity. Both
/// grips have their force zeroed so the imposed velocity is not perturbed by
/// the bond force during the final integration half-step.
fn apply_grips(mut atoms: ResMut<Atom>, registry: Res<AtomDataRegistry>) {
    let pd = registry.expect::<PeriAtom>("apply_grips");
    let nlocal = atoms.nlocal() as usize;
    for i in 0..nlocal {
        let x0 = pd.ref_pos[i][0];
        if x0 < BAR_X_MIN + GRIP {
            atoms.force[i] = [0.0; 3];
            atoms.vel[i] = [0.0; 3];
        } else if x0 > BAR_X_MAX - GRIP {
            atoms.force[i] = [0.0; 3];
            atoms.vel[i] = [PULL_VEL, 0.0, 0.0];
        }
    }
}

/// Print kinetic energy and peak damage every 1000 steps.
fn report(
    atoms: Res<Atom>,
    registry: Res<AtomDataRegistry>,
    run_state: Res<RunState>,
    comm: Res<CommResource>,
) {
    let step = run_state.total_cycle;
    if step % 1000 != 0 {
        return;
    }

    let pd = registry.expect::<PeriAtom>("report");
    let nlocal = atoms.nlocal() as usize;
    let mut ke = 0.0f64;
    let mut dmax = 0.0f64;
    let mut nbroken = 0.0f64;
    for i in 0..nlocal {
        let v = atoms.vel[i];
        ke += 0.5 * atoms.mass[i] * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        dmax = dmax.max(pd.damage[i]);
        nbroken += pd.damage[i];
    }
    let ke = comm.all_reduce_sum_f64(ke);
    let damage_sum = comm.all_reduce_sum_f64(nbroken);

    if comm.rank() == 0 {
        println!(
            "step {step:6}  KE = {ke:.3e} J   max_damage = {dmax:.3}   Σdamage = {damage_sum:.2}"
        );
    }
}
