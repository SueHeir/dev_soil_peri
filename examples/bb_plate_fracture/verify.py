#!/usr/bin/env python3
"""Validate the bond-based mode-I pre-notched plate fracture example.

Builds and runs the dev_soil_peri `bb_plate_fracture` example (single-process,
`--no-default-features`) and checks that:

  1. GRIFFITH RELATION (s0 <-> G0).  The library derives the critical stretch
     from the fracture energy by  s0 = sqrt(5 G0 / (9 K delta)).  We check that
     this inverts exactly (the analytic G0 recovers the configured value) and,
     as corroboration, that a direct discrete sum of the elastic energy released
     by every bond crossing an interior patch of the crack plane at s0 gives the
     fracture energy to the coarse-mesh quadrature accuracy of the bare
     constant-micromodulus model (~18% low at m = delta/dx = 3).

  2. DIGITIZED REFERENCE CRACK CURVE.  Under a super-critical mode-I load the
     crack-tip trajectory is compared to a digitized Ha & Bobaru (2010) crack
     propagation trace.  Their published Figure 11(d) reports crack speed versus
     time for the soda-lime-glass plate; integrating the digitized speed trace
     gives the reference crack-extension-versus-time curve used here.  Because
     this example uses a smaller teaching mesh, the check compares normalized
     extension over the active propagation window, plus the dimensional mean
     speed, within loose qualitative tolerances.

  3. CRACK PROPAGATION EMERGES.  The same super-critical run must advance
     monotonically from the notch across the ligament at a finite, sub-Rayleigh
     speed.  The path is prescribed nowhere -- it comes out of the
     critical-stretch bond breakage.

  4. CORRECT (MODE-I) CRACK PATH.  In the final damage field the crack is a
     straight line on the symmetry plane y = H/2, spanning most of the width.

  5. EMERGENT GRIFFITH THRESHOLD.  A sub-critical applied strain (stored energy
     release rate below G0) leaves the notch dormant -- no propagation -- while
     the super-critical load drives it across.  The crack respects an energy
     threshold; it does not run at any load.

Reference: Silling & Askari, Comput. Struct. 83 (2005) 1526 (3-D bond-based PD,
critical stretch, the s0<->G0 relation); Ha & Bobaru, Int. J. Fract. 162 (2010)
229 (mode-I notch -> straight crack).

Run from the repo root:  $BENCH_PYTHON examples/bb_plate_fracture/verify.py
Exit code 0 = ALL CHECKS PASSED.
"""
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
EXAMPLE = "bb_plate_fracture"
CONFIG = os.path.join("examples", EXAMPLE, "config.toml")

# Material / geometry constants (must match config.toml + main.rs).
G0 = 135.0             # J/m^2  (configured fracture energy)
Y_MID = 0.02           # crack-plane height [m]
A0 = 0.02              # pre-notch length [m]
PLATE_X_MAX = 0.05     # plate width [m]
C_RAYLEIGH = 3160.0    # Rayleigh speed for soda-lime glass ~ 0.92 c_s [m/s]
SUBCRIT_STRAIN = 1.5e-4  # applied far-field strain below the Griffith threshold

# Ha & Bobaru (2010), Fig. 11(d), orange curve: "Crack speed by
# peridynamics (horizon = 0.5 mm)" for soda-lime glass, digitized manually from
# the published plot.  Units are microseconds and m/s.  The paper publishes a
# speed trace rather than a direct a(t); we integrate it below to obtain the
# crack-extension-vs-time reference curve.  This is intentionally static
# reference evidence, not fitted to this example's output.
REFERENCE_SPEED_DIGITIZED = [
    (0.0, 0.0),
    (2.0, 0.0),
    (4.0, 0.0),
    (6.0, 150.0),
    (8.0, 900.0),
    (10.0, 1350.0),
    (12.0, 1375.0),
    (14.0, 1250.0),
    (16.0, 1000.0),
    (18.0, 1180.0),
    (20.0, 1530.0),
    (22.0, 1300.0),
    (24.0, 1680.0),
    (26.0, 1600.0),
    (28.0, 1500.0),
    (30.0, 1380.0),
    (32.0, 1280.0),
    (34.0, 1300.0),
    (36.0, 1330.0),
    (38.0, 1380.0),
    (40.0, 1450.0),
]


def run(cmd, env=None, **kw):
    e = dict(os.environ)
    if env:
        e.update(env)
    return subprocess.run(cmd, cwd=REPO, text=True, capture_output=True, env=e, **kw)


def build():
    print("building bb_plate_fracture (release, --no-default-features)...")
    r = run(["cargo", "build", "--release", "--no-default-features",
             "--example", EXAMPLE])
    if r.returncode != 0:
        print(r.stdout[-2000:]); print(r.stderr[-2000:])
        sys.exit("BUILD FAILED")


def run_case(env=None, timeout=900):
    r = run(["cargo", "run", "--release", "--no-default-features",
             "--example", EXAMPLE, "--", CONFIG], env=env, timeout=timeout)
    if r.returncode != 0:
        print(r.stdout[-3000:]); print(r.stderr[-2000:])
        sys.exit("RUN FAILED")
    return r.stdout


def parse_griffith(out):
    # Match the value right before "J/m^2" (the line also contains a '=' inside
    # the "(= 9 K d s0^2 /5)" annotation, so anchor on the unit, not the '=').
    ana = float(re.search(r"G0_analytic[^\n]*?([\d.eE+-]+)\s*J/m\^2", out).group(1))
    mea = float(re.search(r"G0_measured[^\n]*?([\d.eE+-]+)\s*J/m\^2", out).group(1))
    return ana, mea


def parse_crackrows(out):
    rows = []
    for m in re.finditer(
        r"CRACKROW step=(\d+) time=([\d.eE+-]+) tip_x=([\d.eE+-]+) "
        r"n_crack=([\d.]+) sum_damage=([\d.]+) ke=([\d.eE+-]+)", out):
        rows.append(dict(step=int(m.group(1)), t=float(m.group(2)),
                         tip=float(m.group(3)), ncrack=float(m.group(4)),
                         sumd=float(m.group(5)), ke=float(m.group(6))))
    return rows


def digitized_reference_crack_curve():
    """Return (time_s, extension_m) from the integrated digitized reference."""
    out = []
    extension = 0.0
    last_t_us, last_v = REFERENCE_SPEED_DIGITIZED[0]
    out.append((last_t_us * 1e-6, extension))
    for t_us, v in REFERENCE_SPEED_DIGITIZED[1:]:
        dt = (t_us - last_t_us) * 1e-6
        extension += 0.5 * (last_v + v) * dt
        out.append((t_us * 1e-6, extension))
        last_t_us, last_v = t_us, v
    return out


def interp_curve(points, x):
    if x <= points[0][0]:
        return points[0][1]
    if x >= points[-1][0]:
        return points[-1][1]
    for (x0, y0), (x1, y1) in zip(points, points[1:]):
        if x0 <= x <= x1:
            a = (x - x0) / (x1 - x0)
            return y0 + a * (y1 - y0)
    return points[-1][1]


def normalized_active_curve(points, min_extension):
    active = [(t, x) for t, x in points if x >= min_extension]
    if len(active) < 2:
        return []
    t0, x0 = active[0]
    t1, x1 = active[-1]
    return [((t - t0) / (t1 - t0), (x - x0) / (x1 - x0)) for t, x in active]


def load_field():
    xs, ys, ds = [], [], []
    with open(os.path.join(REPO, "examples", EXAMPLE, "crack_field.csv")) as f:
        next(f)
        for line in f:
            x, y, d = line.split(",")
            xs.append(float(x)); ys.append(float(y)); ds.append(float(d))
    return xs, ys, ds


def maybe_plot(rows, field):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return None
    xs, ys, ds = field
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4))
    sc = a1.scatter(xs, ys, c=ds, s=6, cmap="inferno", vmin=0, vmax=1)
    a1.axhline(Y_MID, color="c", lw=0.6, ls="--")
    a1.set_title("final damage (mid-thickness slice)")
    a1.set_xlabel("x [m]"); a1.set_ylabel("y [m]"); a1.set_aspect("equal")
    fig.colorbar(sc, ax=a1, label="damage")
    a2.plot([r["t"] * 1e6 for r in rows], [r["tip"] * 1e3 for r in rows], "-o", ms=3)
    a2.axhline(A0 * 1e3, color="k", ls=":", label="notch tip")
    a2.set_title("crack tip vs time"); a2.set_xlabel("t [us]")
    a2.set_ylabel("tip x [mm]"); a2.legend()
    fig.tight_layout()
    out = os.path.join(REPO, "examples", EXAMPLE, "crack_path.png")
    fig.savefig(out, dpi=110)
    return out


def main():
    build()

    print("running super-critical case (default eps = 5e-4)...")
    out = run_case()
    ana, mea = parse_griffith(out)
    rows = parse_crackrows(out)
    if not rows:
        sys.exit("no CRACKROW output parsed")
    xs, ys, ds = load_field()   # snapshot the super-critical field before rerun

    print("running sub-critical case (eps = %.1e, below Griffith threshold)..."
          % SUBCRIT_STRAIN)
    sub = run_case(env={"PERI_PRESTRAIN": repr(SUBCRIT_STRAIN)})
    sub_rows = parse_crackrows(sub)
    sub_tip = max(r["tip"] for r in sub_rows) if sub_rows else 0.0

    checks = []

    # --- Check 1: Griffith relation ---
    round_trip = abs(ana - G0) / G0 < 1e-3
    corrob = 0.70 < mea / G0 < 1.30
    checks.append(("Griffith s0=sqrt(5G0/9Kd): G0_analytic=%.3f==%.1f exactly; "
                   "discrete G0_measured=%.1f (%.0f%% of G0, coarse-mesh m=3)"
                   % (ana, G0, mea, 100 * mea / G0), round_trip and corrob))

    # --- Check 2: digitized reference crack-extension curve ---
    ref_curve = digitized_reference_crack_curve()
    ref_norm = normalized_active_curve(ref_curve, min_extension=1.0e-3)
    sim_curve = [
        (r["t"], max(0.0, r["tip"] - A0))
        for r in rows
        if r["tip"] < PLATE_X_MAX - 1.0e-3
    ]
    sim_norm = normalized_active_curve(sim_curve, min_extension=1.0e-3)
    ref_avg_speed = None
    track_rmse = None
    track_max = None
    track_ok = False
    speed_ratio = None
    if len(ref_norm) >= 2 and len(sim_norm) >= 2:
        # Compare the shape of a(t) over each active propagation window.  The
        # reference plate is twice as long in x, so absolute crack length is not
        # the right quantity for this teaching-scale mesh; the dimensional speed
        # check below keeps the comparison tied to physical time.
        err2 = []
        abs_err = []
        for tau, a_sim in sim_norm:
            if 0.0 <= tau <= 1.0:
                a_ref = interp_curve(ref_norm, tau)
                e = a_sim - a_ref
                err2.append(e * e)
                abs_err.append(abs(e))
        if err2:
            track_rmse = (sum(err2) / len(err2)) ** 0.5
            track_max = max(abs_err)
        ref_t0, ref_x0 = ref_curve[0]
        for t, x in ref_curve:
            if x >= 1.0e-3:
                ref_t0, ref_x0 = t, x
                break
        ref_t1, ref_x1 = ref_curve[-1]
        ref_avg_speed = (ref_x1 - ref_x0) / (ref_t1 - ref_t0)

    # --- Check 3: crack propagation emerges (monotone, finite sub-Rayleigh) ---
    tips = [r["tip"] for r in rows]
    monotone = all(b >= a - 1e-9 for a, b in zip(tips, tips[1:]))
    advanced = max(tips) - A0
    prop = [r for r in rows if A0 + 1e-4 < r["tip"] < PLATE_X_MAX - 1e-3]
    speed = None
    if len(prop) >= 2:
        speed = (prop[-1]["tip"] - prop[0]["tip"]) / (prop[-1]["t"] - prop[0]["t"])
    ok = monotone and advanced > 0.015 and speed is not None and 0 < speed < C_RAYLEIGH
    checks.append(("crack propagates: tip %.3f->%.3f m, monotone=%s, "
                   "speed=%.0f m/s (0<v<c_R=%.0f)"
                   % (tips[0], max(tips), monotone, speed or -1, C_RAYLEIGH), ok))

    if speed is not None and ref_avg_speed is not None and track_rmse is not None:
        speed_ratio = speed / ref_avg_speed
        # Qualitative fracture-validation tolerances: <=25% RMS normalized
        # trajectory error, <=40% worst-point deviation, and mean speed within
        # +/-35% of the digitized Ha-Bobaru trace.
        track_ok = track_rmse <= 0.25 and track_max <= 0.40 and 0.65 <= speed_ratio <= 1.35
    checks.append(("digitized Ha-Bobaru Fig. 11(d) crack curve: normalized "
                   "RMSE=%.3f (tol 0.25), max=%.3f (tol 0.40), mean speed "
                   "%.0f vs ref %.0f m/s (ratio %.2f, tol 0.65-1.35)"
                   % (track_rmse or -1, track_max or -1, speed or -1,
                      ref_avg_speed or -1, speed_ratio or -1), track_ok))

    # --- Check 4: correct mode-I path (straight crack on the symmetry plane) ---
    hi = [(x, y) for x, y, d in zip(xs, ys, ds) if d > 0.35]
    on = [(x, y) for x, y in hi if abs(y - Y_MID) < 1.1e-3]
    frac_on = len(on) / max(1, len(hi))
    span = (max(x for x, _ in on) - min(x for x, _ in on)) if on else 0.0
    ok = frac_on > 0.55 and span > 0.6 * PLATE_X_MAX
    checks.append(("mode-I path: %.0f%% of damage on plane y=%.3f, spanning "
                   "%.3f m of width" % (100 * frac_on, Y_MID, span), ok))

    # --- Check 5: emergent Griffith threshold ---
    #  sub-critical load must NOT advance the crack (tip stays at the notch),
    #  while the super-critical load did.
    sub_static = sub_tip <= A0 + 2e-3
    super_ran = max(tips) > A0 + 0.015
    checks.append(("Griffith threshold: sub-critical eps=%.1e leaves tip at "
                   "%.3f m (dormant); super-critical eps=5e-4 runs to %.3f m"
                   % (SUBCRIT_STRAIN, sub_tip, max(tips)),
                   sub_static and super_ran))

    png = maybe_plot(rows, (xs, ys, ds))
    if png:
        print("wrote", os.path.relpath(png, REPO))

    print("\n--- checks ---")
    npass = 0
    for desc, ok in checks:
        print(("  PASS  " if ok else "  FAIL  ") + desc)
        npass += ok
    print("\n%d/%d checks passed" % (npass, len(checks)))
    if npass == len(checks):
        print("ALL CHECKS PASSED")
        sys.exit(0)
    print("CHECKS FAILED")
    sys.exit(1)


if __name__ == "__main__":
    main()
