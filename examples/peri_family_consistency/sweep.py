#!/usr/bin/env python3
"""Validate one-rank vs two-rank peridynamic family consistency."""

import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
EXAMPLE = "peri_family_consistency"
PLOTS = os.path.join(HERE, "plots")

TOL = {
    "entries": 0.0,
    "checksum": 1.0e-6,
    "broken": 0.0,
    "damage_sum": 1.0e-10,
    "crack_count": 0.0,
    "crack_span": 1.0e-12,
}


def run(cmd, **kw):
    return subprocess.run(cmd, cwd=REPO, text=True, capture_output=True, **kw)


def parse(out):
    m = re.search(
        r"FAMILYCHECK ranks=(\d+) entries=([\d.]+) checksum=([\d.eE+-]+) "
        r"broken=([\d.]+) damage_sum=([\d.eE+-]+) crack_count=([\d.]+) "
        r"crack_span=([\d.eE+-]+)",
        out,
    )
    if not m:
        print(out[-3000:])
        sys.exit("missing FAMILYCHECK line")
    return {
        "ranks": int(m.group(1)),
        "entries": float(m.group(2)),
        "checksum": float(m.group(3)),
        "broken": float(m.group(4)),
        "damage_sum": float(m.group(5)),
        "crack_count": float(m.group(6)),
        "crack_span": float(m.group(7)),
    }


def run_case(ranks, config):
    cmd = [
        "mpirun",
        "-np",
        str(ranks),
        "target/release/examples/" + EXAMPLE,
        os.path.join("examples", EXAMPLE, config),
    ]
    r = run(cmd, timeout=180)
    if r.returncode != 0:
        print(r.stdout[-3000:])
        print(r.stderr[-3000:])
        sys.exit(f"{ranks}-rank run failed")
    return parse(r.stdout + r.stderr)


def maybe_plot(one, two, diffs, passed):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except Exception:
        return None

    os.makedirs(PLOTS, exist_ok=True)
    keys = ["entries", "broken", "damage_sum", "crack_count", "crack_span"]
    x = np.arange(len(keys))
    width = 0.36
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))

    axes[0].bar(x - width / 2, [one[k] for k in keys], width, label="1 rank")
    axes[0].bar(x + width / 2, [two[k] for k in keys], width, label="2 ranks")
    axes[0].set_xticks(x, keys, rotation=25, ha="right")
    axes[0].set_title("rank-reduced family and crack metrics")
    axes[0].legend(fontsize=8)

    diff_keys = list(TOL)
    axes[1].bar(diff_keys, [diffs[k] for k in diff_keys], color="#4c78a8", label="abs diff")
    axes[1].plot(diff_keys, [TOL[k] for k in diff_keys], "o", color="#d62728", label="tolerance")
    axes[1].set_yscale("symlog", linthresh=1.0e-12)
    axes[1].set_title("1-rank vs 2-rank differences")
    axes[1].tick_params(axis="x", rotation=25)
    axes[1].legend(fontsize=8)
    axes[1].text(
        0.02,
        0.95,
        "PASS" if passed else "FAIL",
        transform=axes[1].transAxes,
        va="top",
        fontsize=11,
        color="#2ca02c" if passed else "#d62728",
    )

    fig.tight_layout()
    out = os.path.join(PLOTS, "family_consistency.png")
    fig.savefig(out, dpi=160)
    return out


def main():
    print("building peri_family_consistency (release, MPI enabled)...")
    r = run(["cargo", "build", "--release", "--example", EXAMPLE])
    if r.returncode != 0:
        print(r.stdout[-3000:])
        print(r.stderr[-3000:])
        sys.exit("build failed")

    one = run_case(1, "config_1rank.toml")
    two = run_case(2, "config_2rank.toml")
    diffs = {k: abs(one[k] - two[k]) for k in TOL}
    passed = all(diffs[k] <= TOL[k] for k in TOL)
    plot = maybe_plot(one, two, diffs, passed)
    if plot:
        print("wrote", os.path.relpath(plot, REPO))

    print("\n--- family consistency ---")
    for k in TOL:
        status = "PASS" if diffs[k] <= TOL[k] else "FAIL"
        print(f"  {status} {k}: one={one[k]:.12g} two={two[k]:.12g} diff={diffs[k]:.3g} tol={TOL[k]:.3g}")
    if not passed:
        sys.exit("FAMILY CONSISTENCY FAILED")
    print("ALL FAMILY CONSISTENCY CHECKS PASSED")


if __name__ == "__main__":
    main()
