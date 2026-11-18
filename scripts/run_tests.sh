#!/usr/bin/env bash
# Run the QC self-test suite. Run from the project root.
set -euo pipefail
cd "$(dirname "$0")/.."

PY=.venv/bin/python
[ -x "$PY" ] || PY=python3

echo "=== paired-end QC self-tests (logic + HTTP, 60 checks) ==="
"$PY" -I tests/test_qc.py
