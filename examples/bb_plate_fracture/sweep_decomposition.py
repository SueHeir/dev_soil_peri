#!/usr/bin/env python3
"""Fresh external-reference gate plus a rank-decomposition evidence plot."""
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
REPORT = HERE / "decomposition_report.json"
PLOT = HERE / "plots" / "decomposition_invariance.svg"


def run(*command):
    result = subprocess.run(command, cwd=REPO, text=True)
    if result.returncode:
        raise SystemExit(result.returncode)


def main():
    # Delete stale outputs first: the SVG is evidence only when both gates ran now.
    REPORT.unlink(missing_ok=True)
    PLOT.unlink(missing_ok=True)
    run(sys.executable, str(HERE / "verify.py"))
    run(sys.executable, str(HERE / "verify_decomposition.py"),
        "--reference-checked", "--report", str(REPORT))
    evidence = json.loads(REPORT.read_text(encoding="utf-8"))
    damage = evidence["max_abs_difference"]["sum_damage"]
    energy = evidence["max_abs_difference"]["ke"]
    damage_ratio = damage["value"] / damage["bound"]
    energy_ratio = energy["value"] / energy["bound"]
    # A bar past the vertical line is a failed check; exact fields are separately stated.
    PLOT.parent.mkdir(exist_ok=True)
    PLOT.write_text(f'''<svg xmlns="http://www.w3.org/2000/svg" width="760" height="350" viewBox="0 0 760 350">
<rect width="100%" height="100%" fill="white"/><style>text{{font-family:sans-serif;fill:#222}}</style>
<text x="35" y="38" font-size="22" font-weight="bold">dev_soil_peri fracture: 1 vs 2 MPI ranks</text>
<text x="35" y="65" font-size="15">Fresh Ha--Bobaru/Griffith external-reference gate passed before this comparison.</text>
<text x="35" y="93" font-size="15">Every one of {len(evidence['rows'])} full-precision rows: time, crack tip, and cracked-point count exact.</text>
<text x="35" y="145" font-size="16">max |Δ sum damage| / stated f64 reduction bound</text><rect x="35" y="155" width="600" height="28" fill="#eee"/><rect x="35" y="155" width="{min(damage_ratio,1.2)*600:.2f}" height="28" fill="#2878b5"/><line x1="635" y1="145" x2="635" y2="193" stroke="#111" stroke-width="3"/>
<text x="35" y="230" font-size="16">max |Δ kinetic energy| / stated f64 reduction bound</text><rect x="35" y="240" width="600" height="28" fill="#eee"/><rect x="35" y="240" width="{min(energy_ratio,1.2)*600:.2f}" height="28" fill="#d95f02"/><line x1="635" y1="230" x2="635" y2="278" stroke="#111" stroke-width="3"/>
<text x="35" y="320" font-size="14">Black line = pass limit (ratio 1); values: {damage['value']:.3e}/{damage['bound']:.3e}, {energy['value']:.3e}/{energy['bound']:.3e}.</text></svg>\n''', encoding="utf-8")
    print(f"PASS wrote {PLOT.relative_to(REPO)}")


if __name__ == "__main__":
    main()
