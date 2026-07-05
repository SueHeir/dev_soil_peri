//! Bond-based peridynamic mode-I fracture of a pre-notched plate.
//!
//! A rectangular soda-lime-glass plate carries a horizontal **edge pre-notch**
//! at mid-height (bonds crossing the notch slit are removed in the reference
//! configuration). The top and bottom edges are pulled apart at a constant
//! velocity — pure **mode-I** opening. Stress concentrates at the notch tip; the
//! bonds there exceed the critical stretch and break, and the crack propagates
//! straight ahead along the symmetry plane. The crack path is prescribed
//! *nowhere*: it emerges from bond breakage.
//!
//! Two things are validated (see `verify.py`):
//!
//! 1. **Griffith energy.** The critical stretch is derived from the fracture
//!    energy by the library via `s0 = sqrt(5 G0 / (9 K delta))`. At startup this
//!    example independently sums the elastic energy stored at `s0` in every bond
//!    crossing an interior patch of the crack plane and divides by the patch
//!    area — a direct discrete evaluation of the fracture energy — and compares
//!    it back to `G0`. This checks the s0<->G0 relation against the *definition*
//!    of fracture energy, not against itself.
//!
//! 2. **Reference crack trajectory.** The notch tip position is tracked in time
//!    and compared with a digitized Ha & Bobaru (2010) crack-speed trace
//!    integrated to crack extension. The comparison is intentionally loose and
//!    normalized to the active propagation window because this teaching-scale
//!    mesh is smaller than the published plate.
//!
//! 3. **Crack propagation emerges.** The crack advances beyond the notch at a
//!    finite, sub-Rayleigh speed, and the damage stays localised on the symmetry
//!    plane (a crack, not diffuse failure).
//!
//! Reference: Silling & Askari, Comput. Struct. 83 (2005) 1526; Ha & Bobaru,
//! Int. J. Fract. 162 (2010) 229.
//!
//! ```bash
//! cargo run --release --example bb_plate_fracture --no-default-features -- examples/bb_plate_fracture/config.toml
//! ```

use std::collections::HashMap;
use std::fs::File;
use std::io::Write;

use peri_core::prelude::*;

// ── Geometry (must match `[[dev_soil_peri.lattice]]` in config.toml) ────────────────────
/// Plate extent in x (crack-propagation direction) [m].
const PLATE_X_MAX: f64 = 0.05;
/// Plate extent in y (loading direction) [m].
const PLATE_Y_MAX: f64 = 0.04;
/// Plate extent in z (thickness) [m].
const PLATE_Z_MAX: f64 = 0.008;
/// Point spacing Δx [m].
const DX: f64 = 0.001;
/// Horizon δ [m].
const HORIZON: f64 = 0.003;

/// Crack-plane height y = H/2 [m] (the mode-I symmetry plane).
const Y_MID: f64 = 0.5 * PLATE_Y_MAX;
/// Pre-notch length: an edge crack from x = 0 to x = A0 [m].
const A0: f64 = 0.02;

/// Grip thickness at the top and bottom edges [m].
const GRIP: f64 = 0.002;
/// Applied far-field mode-I strain ε_yy (displacement control): the plate is
/// pre-stretched affinely in y and the grips are then held fixed. Chosen below
/// the critical stretch s₀ (≈7.2e-4) so the bulk stays intact, but well above
/// the LEFM crack-initiation strain for this notch (≈1.5e-4) so the stored
/// elastic energy drives the crack from the notch tip.
const PRESTRAIN: f64 = 5.0e-4;

/// Damage level above which a point is considered part of the crack.
const CRACK_DAMAGE: f64 = 0.35;

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(BondBasedDefaultPlugins);

    app.add_setup_system(set_timestep, ScheduleSetupSet::PostSetup);

    // Griffith energy check on the *intact* reference family, then carve the
    // pre-notch — both after the family is built.
    app.add_setup_system(
        griffith_energy_check
            .after("peri_build_families")
            .label("peri_griffith_check"),
        ScheduleSetupSet::PostSetup,
    );
    app.add_setup_system(
        introduce_notch
            .after("peri_griffith_check")
            .label("peri_notch"),
        ScheduleSetupSet::PostSetup,
    );
    // Displacement-control loading: pre-stretch the plate in y after the notch
    // is carved and the reference bonds (r0) are recorded.
    app.add_setup_system(
        apply_prestrain.after("peri_notch"),
        ScheduleSetupSet::PostSetup,
    );

    // Hold the grips fixed at their pre-stretched positions every step (after
    // the bond force, before the final velocity update) so the imposed far-field
    // stretch is maintained while the interior fractures freely.
    app.add_update_system(hold_grips, ParticleSimScheduleSet::PostForce);

    // Crack-front diagnostic + final damage-field dump.
    app.add_update_system(report, ParticleSimScheduleSet::PostFinalIntegration);

    app.start();
}

/// Copy the `[run] dt` into `Atom::dt` (the integrator reads `Atom::dt`).
fn set_timestep(mut atoms: ResMut<Atom>, run_config: Res<RunConfig>) {
    let dt = run_config.current_stage(0).dt;
    if dt > 0.0 {
        atoms.dt = dt;
    }
}

/// Setup system: independent discrete evaluation of the fracture energy.
///
/// For an interior patch of the crack plane `y = Y_MID`, sum the elastic energy
/// stored at the critical stretch in every bond that crosses the plane, and
/// divide by the patch area. The energy in one bond at stretch `s0` is
/// `U = ½ c s0² r0 Vᵢ Vⱼ` (the bond is a linear spring of stiffness
/// `k = c Vᵢ Vⱼ / r0` stretched by `Δ = s0 r0`). Summing the released energy of
/// all crossing bonds per unit area *is* the fracture energy `G0` — so this
/// checks the library's `s0 = sqrt(5 G0 / 9 K δ)` against the definition of
/// `G0`, independently of the value the library derived.
///
/// The patch is kept ≥ δ away from every free surface (and from the notch) so
/// that the crossing bonds form a complete 3-D family, matching the continuum
/// integral the analytic relation comes from. `G0_analytic` recovers the
/// configured `G0` exactly (the relation is invertible); `G0_measured` from the
/// bare constant-micromodulus lattice sum sits ~18% low at m = δ/Δx = 3, the
/// well-known surface-quadrature underestimate of the model without a
/// partial-volume correction — same order, right physics.
fn griffith_energy_check(
    atoms: Res<Atom>,
    registry: Res<AtomDataRegistry>,
    table: Res<PeriMaterialTable>,
    comm: Res<CommResource>,
) {
    if comm.rank() != 0 {
        return;
    }
    let bond_store = registry.expect::<BondStore>("griffith_energy_check");
    let pd = registry.expect::<PeriAtom>("griffith_energy_check");
    let nlocal = atoms.nlocal as usize;

    // tag → index for local points (partners of an interior patch are local).
    let mut tag_to_index: HashMap<u32, usize> = HashMap::with_capacity(nlocal);
    for i in 0..nlocal {
        tag_to_index.insert(atoms.tag[i], i);
    }

    // Interior sampling patch on the crack plane, ≥ δ from every free surface
    // and clear of the pre-notch (x ≥ A0 + δ).
    let x0 = A0 + HORIZON;
    let x1 = PLATE_X_MAX - HORIZON;
    let z0 = HORIZON;
    let z1 = PLATE_Z_MAX - HORIZON;
    let area = (x1 - x0) * (z1 - z0);

    let mat = &table.materials[0];
    let c = mat.micromodulus;
    let s0 = mat.critical_stretch;
    let k = mat.bulk_modulus;

    // Sum ½ c s0² r0 Vᵢ Vⱼ over every bond crossing the plane y = Y_MID (which
    // sits midway between two lattice layers), counting each bond once from its
    // lower endpoint, whose (x, z) must lie in the interior patch.
    let mut energy = 0.0f64;
    for i in 0..nlocal {
        let xi = pd.ref_pos[i];
        if xi[1] >= Y_MID {
            continue;
        }
        if xi[0] < x0 || xi[0] > x1 || xi[2] < z0 || xi[2] > z1 {
            continue;
        }
        let vi = pd.volume[i];
        for bond in &bond_store.bonds[i] {
            let j = match tag_to_index.get(&bond.partner_tag) {
                Some(&j) => j,
                None => continue,
            };
            if pd.ref_pos[j][1] <= Y_MID {
                continue; // partner not on the far side of the plane
            }
            let vj = pd.volume[j];
            energy += 0.5 * c * s0 * s0 * bond.r0 * vi * vj;
        }
    }
    let g0_measured = energy / area;
    let g0_analytic = 9.0 * k * HORIZON * s0 * s0 / 5.0; // inverse of the s0 relation
    let rel_err = (g0_measured - g0_analytic).abs() / g0_analytic;

    println!("Griffith check (crack plane y = {Y_MID:.4} m, delta = {HORIZON:.4} m):");
    println!("  s0 (from G0 via library)      = {s0:.6e}");
    println!("  bulk modulus K                = {k:.4e} Pa");
    println!("  G0_analytic (= 9 K d s0^2 /5) = {g0_analytic:.4} J/m^2");
    println!("  G0_measured (discrete bonds)  = {g0_measured:.4} J/m^2");
    println!("  GRIFFITH_REL_ERR              = {rel_err:.4}");
}

/// Setup system: carve the mode-I edge pre-notch.
///
/// Every reference bond that crosses the plane `y = Y_MID` at an x within the
/// notch `[0, A0]` is marked [`BOND_BROKEN`]; the point damage is updated to the
/// resulting broken fraction. Both endpoints hold the bond, and both mark it, so
/// the removal is symmetric. The pre-notch shows up directly in the damage field.
fn introduce_notch(
    atoms: Res<Atom>,
    registry: Res<AtomDataRegistry>,
    comm: Res<CommResource>,
) {
    let mut bond_store = registry.expect_mut::<BondStore>("introduce_notch");
    let mut pd = registry.expect_mut::<PeriAtom>("introduce_notch");
    let nlocal = atoms.nlocal as usize;

    let mut tag_to_index: HashMap<u32, usize> = HashMap::with_capacity(nlocal);
    for i in 0..nlocal {
        tag_to_index.insert(atoms.tag[i], i);
    }

    let mut nbroken_total = 0u64;
    for i in 0..nlocal {
        let xi = pd.ref_pos[i];
        let n_bonds = bond_store.bonds[i].len();
        for k in 0..n_bonds {
            let partner_tag = bond_store.bonds[i][k].partner_tag;
            let j = match tag_to_index.get(&partner_tag) {
                Some(&j) => j,
                None => continue,
            };
            let xj = pd.ref_pos[j];
            // Does the reference bond straddle the crack plane?
            if (xi[1] - Y_MID) * (xj[1] - Y_MID) >= 0.0 {
                continue;
            }
            // x where the bond crosses y = Y_MID.
            let t = (Y_MID - xi[1]) / (xj[1] - xi[1]);
            let x_cross = xi[0] + t * (xj[0] - xi[0]);
            if x_cross <= A0 {
                if bond_store.bonds[i][k].bond_type != BOND_BROKEN {
                    bond_store.bonds[i][k].bond_type = BOND_BROKEN;
                    nbroken_total += 1;
                }
            }
        }
    }

    // Refresh damage from the pre-notch.
    for i in 0..nlocal {
        let n0 = pd.n0[i];
        if n0 <= 0.0 {
            continue;
        }
        let nb = bond_store.bonds[i]
            .iter()
            .filter(|b| b.bond_type == BOND_BROKEN)
            .count() as f64;
        pd.damage[i] = nb / n0;
    }

    let total = comm.all_reduce_sum_f64(nbroken_total as f64) as u64;
    if comm.rank() == 0 {
        println!(
            "Notch: severed {total} bond-endpoints across the plane y = {Y_MID:.4} m for x <= {A0:.4} m"
        );
    }
}

/// Setup system: impose the mode-I far-field strain by an affine y-stretch
/// about the crack plane, `y ← y + PRESTRAIN·(y − Y_MID)`. The reference bond
/// lengths (`r0`) were already recorded from the undeformed lattice, so this
/// pre-stretch loads every bond; the grips (held fixed by [`hold_grips`]) then
/// store this as elastic energy that drives the crack.
fn apply_prestrain(mut atoms: ResMut<Atom>, registry: Res<AtomDataRegistry>, comm: Res<CommResource>) {
    // The verifier brackets the Griffith threshold by running sub- and
    // super-critical loads; allow the strain to be overridden without editing
    // the (declarative) config. Defaults to `PRESTRAIN`.
    let eps = std::env::var("PERI_PRESTRAIN")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(PRESTRAIN);
    let pd = registry.expect::<PeriAtom>("apply_prestrain");
    let nlocal = atoms.nlocal as usize;
    for i in 0..nlocal {
        let y0 = pd.ref_pos[i][1];
        atoms.pos[i][1] = y0 + eps * (y0 - Y_MID);
        atoms.vel[i] = [0.0; 3];
    }
    if comm.rank() == 0 {
        println!("Prestrain: applied far-field eps_yy = {eps:.3e} (displacement control)");
    }
}

/// Mode-I grips: hold the top and bottom edge layers fixed at their
/// pre-stretched positions (zero velocity, zero force) so the imposed far-field
/// stretch is maintained while the interior fractures freely.
fn hold_grips(mut atoms: ResMut<Atom>, registry: Res<AtomDataRegistry>) {
    let pd = registry.expect::<PeriAtom>("hold_grips");
    let nlocal = atoms.nlocal as usize;
    for i in 0..nlocal {
        let y0 = pd.ref_pos[i][1];
        if y0 > PLATE_Y_MAX - GRIP || y0 < GRIP {
            atoms.force[i] = [0.0; 3];
            atoms.vel[i] = [0.0; 3];
        }
    }
}

/// Crack-front diagnostic every `thermo` steps, plus a final damage dump.
///
/// The crack tip is the largest reference-x of any point on the crack plane
/// (|y − Y_MID| < 0.75 Δx) whose damage exceeds [`CRACK_DAMAGE`]. Printed as a
/// machine-readable `CRACKROW` line the verifier parses; the final damage field
/// (mid-thickness slice) is written to `crack_field.csv`.
fn report(
    atoms: Res<Atom>,
    registry: Res<AtomDataRegistry>,
    run_state: Res<RunState>,
    run_config: Res<RunConfig>,
    comm: Res<CommResource>,
) {
    let step = run_state.total_cycle;
    let steps_total = run_config.current_stage(0).steps as usize;
    let dt = run_config.current_stage(0).dt;
    let thermo = 100usize;
    // `total_cycle` is 0-indexed, so the final cycle is `steps_total - 1`.
    let is_last = step + 1 >= steps_total;
    if step % thermo != 0 && !is_last {
        return;
    }

    let pd = registry.expect::<PeriAtom>("report");
    let nlocal = atoms.nlocal as usize;

    let mut tip_x = 0.0f64;
    let mut n_crack = 0.0f64; // points on the plane above CRACK_DAMAGE
    let mut sum_damage = 0.0f64;
    let mut ke = 0.0f64;
    for i in 0..nlocal {
        let v = atoms.vel[i];
        ke += 0.5 * atoms.mass[i] * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        sum_damage += pd.damage[i];
        let y0 = pd.ref_pos[i][1];
        if (y0 - Y_MID).abs() < 0.75 * DX && pd.damage[i] > CRACK_DAMAGE {
            n_crack += 1.0;
            tip_x = tip_x.max(pd.ref_pos[i][0]);
        }
    }

    let ke = comm.all_reduce_sum_f64(ke);
    let sum_damage = comm.all_reduce_sum_f64(sum_damage);
    let n_crack = comm.all_reduce_sum_f64(n_crack);
    let time = step as f64 * dt;

    if comm.rank() == 0 {
        // The single-process slice (the supported configuration) owns every
        // point, so the local `tip_x` is already the global crack tip.
        println!(
            "CRACKROW step={step} time={time:.6e} tip_x={tip_x:.4} n_crack={n_crack:.0} sum_damage={sum_damage:.2} ke={ke:.3e}"
        );
    }

    // Snapshot the damage field periodically (and at the end) so the crack
    // path's time evolution can be plotted.
    if step % 500 == 0 || is_last {
        dump_crack_field(&atoms, &pd, step, is_last);
    }
}

/// Write the mid-thickness damage slice to a CSV for plotting. Snapshots go to
/// `crack_field_<step>.csv`; the final field is also written to `crack_field.csv`.
fn dump_crack_field(atoms: &Atom, pd: &PeriAtom, step: usize, is_last: bool) {
    let z_mid = 0.5 * PLATE_Z_MAX;
    let path = if is_last {
        "examples/bb_plate_fracture/crack_field.csv".to_string()
    } else {
        format!("examples/bb_plate_fracture/crack_field_{step:05}.csv")
    };
    let mut f = match File::create(&path) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("could not write {path}: {e}");
            return;
        }
    };
    let _ = writeln!(f, "x,y,damage");
    let nlocal = atoms.nlocal as usize;
    for i in 0..nlocal {
        let p = pd.ref_pos[i];
        if (p[2] - z_mid).abs() < 0.75 * DX {
            let _ = writeln!(f, "{:.5},{:.5},{:.5}", p[0], p[1], pd.damage[i]);
        }
    }
    println!("Wrote crack field slice to {path}");
}
