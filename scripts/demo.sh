#!/usr/bin/env bash
# End-to-end demo against a running QC review service (default
# http://localhost:8000). Generates paired FASTQ + candidate JSON, then shows:
#   1. a clean calibrated run
#   2. cross-reference / orientation / span / undetermined anomalies
#   3. whole-request failures (corrupted pair file, bad candidate,
#      insufficient calibration pairs)
set -euo pipefail

BASE_URL="${BASE_URL:-http://localhost:8000}"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

need() { command -v "$1" >/dev/null || { echo "missing tool: $1" >&2; exit 1; }; }
need curl
need python3

echo "== health =="
curl -fsS "$BASE_URL/health"
echo

python3 - "$WORK" <<'PY'
import json
import os
import sys

work = sys.argv[1]
READ_LEN = 100

def fastq(path, side, ids):
    with open(path, "wb") as fh:
        for rid in ids:
            seq = (b"A" if side == 1 else b"C") * READ_LEN
            fh.write(b"@%s/%d desc\n" % (rid.encode(), side) + seq
                     + b"\n+\n" + b"~" * READ_LEN + b"\n")

def c(rid, side, ref, start, strand="+", score=50):
    return {"read_id": rid, "side": side, "reference": ref,
            "start": start, "strand": strand, "score": score}

# 24 clean inward calibration pairs, span 300 (R1 start 1000,
# R2 reverse 5' end at 1200 -> rightmost base 1299).
ids = [f"r{i:04d}" for i in range(24)]
fastq(os.path.join(work, "R1.fastq"), 1, ids)
fastq(os.path.join(work, "R2.fastq"), 2, ids)
cands = []
for rid in ids:
    cands.append(c(rid, 1, "chr1", 1000))
    cands.append(c(rid, 2, "chr1", 1200, "-"))

# Anomalous pairs appended after the calibration set.
extra = ["xref", "orient", "spanout", "tied"]
fastq(os.path.join(work, "R1.fastq"), 1, ids + extra)
fastq(os.path.join(work, "R2.fastq"), 2, ids + extra)
cands += [
    c("xref", 1, "chr1", 1000), c("xref", 2, "chr2", 1200, "-"),
    c("orient", 1, "chr1", 2000, "-"), c("orient", 2, "chr1", 1000, "+"),
    # spanout carries a low-scoring second R1 candidate so it cannot enter
    # the single-candidate calibration set; the best combination (60+50) is
    # inward with span 900 vs calibrated median 300 -> span anomaly.
    c("spanout", 1, "chr1", 1000, "+", 60),
    c("spanout", 1, "chr1", 1100, "+", 5),
    c("spanout", 2, "chr1", 1700, "-", 50),
    # Tied R1 choices: inward normal on chr1 vs cross-reference to chr2,
    # equal combined score -> undetermined.
    c("tied", 1, "chr1", 1000, "+", 25),
    c("tied", 1, "chr2", 1000, "+", 25),
    c("tied", 2, "chr1", 1200, "-", 25),
]
with open(os.path.join(work, "candidates.json"), "w") as fh:
    json.dump({"candidates": cands}, fh)

# Corrupted pair file (truncated trailing record).
with open(os.path.join(work, "R1_bad.fastq"), "wb") as fh:
    fh.write(open(os.path.join(work, "R1.fastq"), "rb").read())
    fh.write(b"@tail/1\nACGT\n")

# Candidate document with a non-positive start.
bad_cands = json.loads(json.dumps(cands))
bad_cands[0]["start"] = 0
with open(os.path.join(work, "candidates_bad.json"), "w") as fh:
    json.dump({"candidates": bad_cands}, fh)

# Only 19 pairs -> calibration failure.
few_ids = ids[:19]
fastq(os.path.join(work, "few_R1.fastq"), 1, few_ids)
fastq(os.path.join(work, "few_R2.fastq"), 2, few_ids)
few_cands = [x for x in cands if x["read_id"] in set(few_ids)]
with open(os.path.join(work, "few_candidates.json"), "w") as fh:
    json.dump({"candidates": few_cands}, fh)
PY

post() { # label r1 r2 cands expected_status
  local label="$1" r1="$2" r2="$3" cands="$4" expected="$5"
  echo "== $label (expect HTTP $expected) =="
  curl -s -o "$WORK/resp.json" -w 'HTTP %{http_code}\n' \
    -X POST "$BASE_URL/qc/analyze" \
    -F "task_id=demo-paired-task-001" \
    -F "r1=@$r1;filename=R1.fastq" \
    -F "r2=@$r2;filename=R2.fastq" \
    -F "candidates=@$cands;filename=candidates.json;type=application/json"
  python3 -m json.tool "$WORK/resp.json"
}

post "full QC run with anomalies" \
     "$WORK/R1.fastq" "$WORK/R2.fastq" "$WORK/candidates.json" 200

post "corrupted pair file -> whole request fails" \
     "$WORK/R1_bad.fastq" "$WORK/R2.fastq" "$WORK/candidates.json" 400

post "non-positive candidate start -> whole request fails" \
     "$WORK/R1.fastq" "$WORK/R2.fastq" "$WORK/candidates_bad.json" 400

post "fewer than 20 calibration pairs -> calibration failure" \
     "$WORK/few_R1.fastq" "$WORK/few_R2.fastq" "$WORK/few_candidates.json" 400

echo "ALL DEMO STEPS DONE"
