#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

: "${BENCH_PYTHON:=python3}"

echo "=== dev_soil_peri validation set ==="
"$BENCH_PYTHON" examples/bb_plate_fracture/verify.py
echo "ALL VALIDATIONS PASSED"
