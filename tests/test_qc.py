"""Self-tests for the paired-end fragment QC review service.

Covers both the pure QC logic (app.qc / app.fastq) and the full HTTP surface
with FastAPI's in-process test client. Run:

    python -I tests/test_qc.py
"""

from __future__ import annotations

import io
import json
import sys
from pathlib import Path

# Make the project root importable when executed as a script.
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fastapi.testclient import TestClient  # noqa: E402

from app import qc  # noqa: E402
from app.fastq import parse_fastq, validate_paired_files  # noqa: E402
from app.main import app  # noqa: E402

RESULTS: list[tuple[str, bool, str]] = []
client = TestClient(app)


def record(name: str, ok: bool, detail: str = "") -> None:
    RESULTS.append((name, ok, detail))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f" - {detail}" if detail else ""))


def check(name: str, condition: bool, detail: str = "") -> None:
    record(name, bool(condition), detail)


# --------------------------------------------------------------------------- #
# Fixture builders
# --------------------------------------------------------------------------- #

def make_fastq(records: list[tuple[str, str]], side: int) -> bytes:
    out = io.BytesIO()
    for read_id, seq in records:
        out.write(
            b"@" + read_id.encode() + b"/" + str(side).encode() + b" desc\n"
            + seq.encode() + b"\n+\n" + b"~" * len(seq) + b"\n"
        )
    return out.getvalue()


def cand(read_id: str, side: int, reference: str, start: int,
         strand: str = "+", score: int = 10) -> dict:
    return {
        "read_id": read_id,
        "side": side,
        "reference": reference,
        "start": start,
        "strand": strand,
        "score": score,
    }


def pairs_from(r1: bytes, r2: bytes) -> list[qc.PairInfo]:
    tuples = validate_paired_files(parse_fastq(r1, 1), parse_fastq(r2, 2))
    return [qc.PairInfo(rid, l1, l2) for rid, l1, l2 in tuples]


def make_clean_dataset(
    n: int = 25, read_len: int = 100, gap: int = 100, ref: str = "chr1"
) -> tuple[list[qc.PairInfo], list[dict]]:
    """n inward pairs with constant span = 2*read_len + gap, all normal."""
    ids = [f"r{i:04d}" for i in range(n)]
    r1 = make_fastq([(rid, "A" * read_len) for rid in ids], 1)
    r2 = make_fastq([(rid, "C" * read_len) for rid in ids], 2)
    pairs = pairs_from(r1, r2)
    candidates: list[dict] = []
    for rid in ids:
        r1_start = 1000
        r2_start = r1_start + read_len + gap  # reverse mate 5' end
        candidates.append(cand(rid, 1, ref, r1_start, "+", 50))
        candidates.append(cand(rid, 2, ref, r2_start, "-", 50))
    return pairs, candidates


# --------------------------------------------------------------------------- #
# 1. Pure-logic tests
# --------------------------------------------------------------------------- #

def test_calibration_median_mad_tolerance() -> None:
    pairs, candidates = make_clean_dataset(25, read_len=100, gap=100)
    report = qc.analyze(pairs, qc.parse_candidates(json.dumps(
        {"candidates": candidates}).encode()))
    cal = report["calibration"]
    # span = 100 + 100 + 100 = 300 for every pair
    check("calibration: median span 300", cal["median_span"] == 300.0,
          str(cal))
    check("calibration: MAD 0 on constant spans", cal["mad_span"] == 0.0)
    check("calibration: tolerance max(50, 0) = 50",
          cal["span_tolerance"] == 50.0)
    check("calibration: 25 calibration pairs", cal["pair_count"] == 25)
    check("calibration: all 25 normal", report["counts"][qc.NORMAL] == 25)
    check("calibration: no anomalies", report["anomalies"] == [])


def test_span_anomaly_inside_and_outside_tolerance() -> None:
    pairs, candidates = make_clean_dataset(24, read_len=50, gap=200)
    # 24 clean pairs, span = 300; append one span-anomaly pair.
    rid = "span-out"
    pairs.append(qc.PairInfo(rid, 50, 50))
    candidates.append(cand(rid, 1, "chr1", 1000, "+", 60))
    candidates.append(cand(rid, 2, "chr1", 1000 + 50 + 400, "-", 60))  # span 500
    report = qc.analyze(pairs, qc.parse_candidates(json.dumps(
        {"candidates": candidates}).encode()))
    check("span: deviation 200 > tolerance 50 -> span anomaly",
          report["counts"][qc.SPAN] == 1 and len(report["anomalies"]) == 1
          and report["anomalies"][0]["category"] == qc.SPAN)

    # A deviation exactly at the tolerance boundary is normal (<=).
    rid2 = "span-edge"
    pairs2 = pairs[:24] + [qc.PairInfo(rid2, 50, 50)]
    cands2 = candidates[:48] + [
        cand(rid2, 1, "chr1", 1000, "+", 60),
        cand(rid2, 2, "chr1", 1000 + 50 + 250, "-", 60),  # span 350, |delta|=50
    ]
    report2 = qc.analyze(pairs2, qc.parse_candidates(json.dumps(
        {"candidates": cands2}).encode()))
    check("span: |delta| == tolerance (50) counts as normal",
          report2["counts"][qc.NORMAL] == 25 and report2["anomalies"] == [])


def test_mad_based_tolerance() -> None:
    # 21 pairs: 20 with span 300, one with span 600 -> median 300, MAD 0,
    # tolerance stays at the absolute floor of 50.
    pairs, candidates = make_clean_dataset(20, read_len=100, gap=100)
    rid = "wide"
    pairs.append(qc.PairInfo(rid, 100, 100))
    candidates += [cand(rid, 1, "chr1", 1000, "+", 9),
                   cand(rid, 2, "chr1", 1000 + 100 + 300, "-", 9)]  # span 500
    report = qc.analyze(pairs, qc.parse_candidates(json.dumps(
        {"candidates": candidates}).encode()))
    check("MAD floor: tolerance stays 50 when MAD is 0",
          report["calibration"]["span_tolerance"] == 50.0)
    check("MAD floor: wide pair flagged span",
          report["counts"][qc.SPAN] == 1)


def test_calibration_needs_20_pairs() -> None:
    pairs, candidates = make_clean_dataset(19)
    try:
        qc.analyze(pairs, qc.parse_candidates(json.dumps(
            {"candidates": candidates}).encode()))
        check("calibration: 19 clean pairs rejected", False)
    except qc.CalibrationError as exc:
        check("calibration: 19 clean pairs rejected", "at least 20" in str(exc),
              str(exc))

    # 19 clean inward pairs plus one cross-reference pair do not rescue it:
    # only single-candidate inward pairs count toward calibration.
    rid = "xref"
    pairs.append(qc.PairInfo(rid, 50, 50))
    candidates += [cand(rid, 1, "chr1", 10, "+", 1),
                   cand(rid, 2, "chr2", 10, "-", 1)]
    try:
        qc.analyze(pairs, qc.parse_candidates(json.dumps(
            {"candidates": candidates}).encode()))
        check("calibration: non-inward pair excluded from calibration", False)
    except qc.CalibrationError:
        check("calibration: non-inward pair excluded from calibration", True)


def test_cross_reference_and_orientation() -> None:
    pairs, candidates = make_clean_dataset(22)
    xref, orient = "xref", "orient"
    pairs += [qc.PairInfo(xref, 50, 50), qc.PairInfo(orient, 50, 50)]
    candidates += [
        cand(xref, 1, "chr1", 100, "+", 30),
        cand(xref, 2, "chr2", 100, "-", 30),
        # Same reference but outward: R1 reverse on the right of R2 forward.
        cand(orient, 1, "chr1", 500, "-", 30),
        cand(orient, 2, "chr1", 100, "+", 30),
    ]
    report = qc.analyze(pairs, qc.parse_candidates(json.dumps(
        {"candidates": candidates}).encode()))
    check("classification: cross_reference counted",
          report["counts"][qc.CROSS_REFERENCE] == 1)
    check("classification: orientation counted",
          report["counts"][qc.ORIENTATION] == 1)
    cats = {a["read_id"]: a["category"] for a in report["anomalies"]}
    check("classification: anomaly categories attached to read ids",
          cats.get(xref) == qc.CROSS_REFERENCE
          and cats.get(orient) == qc.ORIENTATION, str(cats))


def test_best_combined_score_wins() -> None:
    pairs, candidates = make_clean_dataset(20)
    rid = "multi"
    pairs.append(qc.PairInfo(rid, 50, 50))
    # R1: chr1 forward (score 60, normal placement) vs chr2 forward (score 5).
    candidates += [
        cand(rid, 1, "chr1", 1000, "+", 60),
        cand(rid, 1, "chr2", 1000, "+", 5),
        cand(rid, 2, "chr1", 1000 + 50 + 200, "-", 40),
    ]
    report = qc.analyze(pairs, qc.parse_candidates(json.dumps(
        {"candidates": candidates}).encode()))
    # Best combination = 60+40 on chr1 (normal); cross-reference combos lose.
    check("scoring: highest combined score chosen, pair normal",
          report["counts"][qc.NORMAL] == 21 and report["anomalies"] == [])

    # Flip the scores: the cross-reference placement now wins.
    for c in candidates:
        if c["read_id"] == rid and c["side"] == 1 and c["reference"] == "chr1":
            c["score"] = 5
        if c["read_id"] == rid and c["side"] == 1 and c["reference"] == "chr2":
            c["score"] = 60
    report2 = qc.analyze(pairs, qc.parse_candidates(json.dumps(
        {"candidates": candidates}).encode()))
    check("scoring: lower-scoring normal combo loses to cross-reference",
          report2["counts"][qc.CROSS_REFERENCE] == 1)


def test_ties_consistent_verdict_adopted() -> None:
    pairs, candidates = make_clean_dataset(20)
    rid = "tied-same"
    pairs.append(qc.PairInfo(rid, 50, 50))
    # Two tied R1 choices on different references; both combinations are
    # cross-reference because R2 sits on chr9 -> verdict consistent.
    candidates += [
        cand(rid, 1, "chr1", 100, "+", 20),
        cand(rid, 1, "chr2", 100, "+", 20),
        cand(rid, 2, "chr9", 100, "-", 20),
    ]
    report = qc.analyze(pairs, qc.parse_candidates(json.dumps(
        {"candidates": candidates}).encode()))
    check("tie: consistent cross_reference verdicts adopted",
          report["counts"][qc.CROSS_REFERENCE] == 1
          and report["counts"][qc.UNDETERMINED] == 0)


def test_ties_conflicting_verdict_undetermined() -> None:
    pairs, candidates = make_clean_dataset(20)
    rid = "tied-conflict"
    pairs.append(qc.PairInfo(rid, 50, 50))
    # Tied best score: R1 chr1 (inward normal with R2 chr1) vs R1 chr2
    # (cross-reference with R2 chr1).
    candidates += [
        cand(rid, 1, "chr1", 1000, "+", 25),
        cand(rid, 1, "chr2", 1000, "+", 25),
        cand(rid, 2, "chr1", 1000 + 50 + 200, "-", 25),
    ]
    report = qc.analyze(pairs, qc.parse_candidates(json.dumps(
        {"candidates": candidates}).encode()))
    check("tie: conflicting verdicts -> undetermined",
          report["counts"][qc.UNDETERMINED] == 1
          and report["anomalies"][-1]["category"] == qc.UNDETERMINED)


def test_tie_normal_vs_span_is_undetermined() -> None:
    # 20 calibration pairs with MAD large enough for a wide tolerance, then a
    # pair whose tied combos disagree between normal and span.
    pairs, candidates = make_clean_dataset(20, read_len=50, gap=100)  # span 200
    # Widen MAD: perturb calibration starts so MAD = 25 -> tolerance 100.
    for i, c in enumerate(candidates):
        if c["side"] == 2 and (i // 2) % 2 == 1:
            c["start"] += 50
    # median stays 200; MAD = 25 -> tolerance = max(50, 100) = 100.
    rid = "edge-tie"
    pairs.append(qc.PairInfo(rid, 50, 50))
    candidates += [
        cand(rid, 1, "chr1", 1000, "+", 10),
        cand(rid, 2, "chr1", 1000 + 50 + 100, "-", 10),   # span 200 normal
        cand(rid, 2, "chr1", 1000 + 50 + 350, "-", 10),   # span 450, delta 250 -> span
    ]
    report = qc.analyze(pairs, qc.parse_candidates(json.dumps(
        {"candidates": candidates}).encode()))
    check("tie: normal vs span conflict -> undetermined",
          report["counts"][qc.UNDETERMINED] == 1, str(report["counts"]))
    check("tie: tolerance derived from 4xMAD (100 > floor 50)",
          report["calibration"]["span_tolerance"] == 100.0,
          str(report["calibration"]))


def test_input_validation_errors() -> None:
    def expect_error(name: str, pairs: list[qc.PairInfo], raw_cands: bytes):
        try:
            qc.analyze(pairs, qc.parse_candidates(raw_cands))
            check(name, False, "no error raised")
        except qc.QcInputError as exc:
            check(name, True, str(exc)[:80])

    pairs, candidates = make_clean_dataset(20)
    good = json.dumps({"candidates": candidates}).encode()

    # Missing candidate for one end.
    broken = [c for c in candidates
              if not (c["read_id"] == "r0000" and c["side"] == 2)]
    expect_error("validation: missing R2 candidate -> whole error",
                 pairs, json.dumps({"candidates": broken}).encode())

    # Unknown candidate id.
    extra = candidates + [cand("ghost", 1, "chr1", 5, "+")]
    expect_error("validation: unknown candidate id -> whole error",
                 pairs, json.dumps({"candidates": extra}).encode())

    # Non-positive start.
    pairs2 = pairs + [qc.PairInfo("fresh", 50, 50)]
    bad2 = candidates + [
        cand("fresh", 1, "chr1", 0, "+", 1),
        cand("fresh", 2, "chr1", 60, "-", 1),
    ]
    expect_error("validation: start < 1 rejected", pairs2,
                 json.dumps({"candidates": bad2}).encode())

    # Duplicate candidate identity.
    dup = candidates + [dict(candidates[0])]
    expect_error("validation: duplicate candidate rejected", pairs,
                 json.dumps({"candidates": dup}).encode())

    # Four candidates for one end (distinct identities so the per-end cap,
    # not the duplicate rule, is what fires).
    four = candidates + [
        cand("r0000", 1, "chr1", 900, "+", 1),
        cand("r0000", 1, "chr1", 901, "+", 1),
        cand("r0000", 1, "chr2", 100, "+", 1),
    ]
    try:
        qc.analyze(pairs, qc.parse_candidates(json.dumps(
            {"candidates": four}).encode()))
        check("validation: >3 candidates per end rejected", False)
    except qc.QcInputError as exc:
        check("validation: >3 candidates per end rejected",
              "at most 3" in str(exc), str(exc)[:100])

    # Malformed JSON and wrong top-level shape.
    expect_error("validation: malformed JSON rejected", pairs, b"{not json")
    expect_error("validation: wrong document shape rejected", pairs,
                 json.dumps(candidates).encode())

    # Bad enum / type cases.
    base = candidates[:]
    c = dict(candidates[0]); c["strand"] = "?"
    expect_error("validation: illegal strand rejected", pairs,
                 json.dumps({"candidates": base + [c]}).encode())
    c = dict(candidates[0]); c["score"] = 3.5
    expect_error("validation: non-integer score rejected", pairs,
                 json.dumps({"candidates": base + [c]}).encode())
    c = dict(candidates[0]); c["side"] = 3
    expect_error("validation: side outside {1,2} rejected", pairs,
                 json.dumps({"candidates": base + [c]}).encode())
    expect_error("validation: boolean disguised as integer rejected", pairs,
                 json.dumps({"candidates":
                             base + [cand("r0001", 1, "chr1", True, "+", 1)]}).encode())
    # good dataset itself must pass
    check("validation: clean dataset accepted",
          qc.analyze(pairs, qc.parse_candidates(good))["pair_count"] == 20)


def test_pair_and_length_limits() -> None:
    # 2001 pairs -> rejected.
    ids = [f"q{i:05d}" for i in range(2001)]
    r1 = make_fastq([(rid, "A" * 20) for rid in ids], 1)
    r2 = make_fastq([(rid, "C" * 20) for rid in ids], 2)
    pairs = pairs_from(r1, r2)
    cands: list[dict] = []
    for rid in ids:
        cands += [cand(rid, 1, "chr1", 10, "+", 1),
                  cand(rid, 2, "chr1", 200, "-", 1)]
    try:
        qc.analyze(pairs, qc.parse_candidates(json.dumps(
            {"candidates": cands}).encode()))
        check("limit: 2001 pairs rejected", False)
    except qc.QcInputError as exc:
        check("limit: 2001 pairs rejected", "2000" in str(exc), str(exc))

    # Exactly 2000 pairs accepted at the validation layer (calibration may
    # still fail in this synthetic layout, which is a different error path).
    pairs_ok = pairs[:2000]
    cands_ok = [c for c in cands if c["read_id"] in {p.read_id for p in pairs_ok}]
    report = qc.analyze(pairs_ok, qc.parse_candidates(json.dumps(
        {"candidates": cands_ok}).encode()))
    check("limit: 2000 pairs accepted", report["pair_count"] == 2000)

    # 501-base read rejected.
    r1 = make_fastq([("long", "A" * 501)], 1)
    r2 = make_fastq([("long", "C" * 20)], 2)
    long_pairs = pairs_from(r1, r2)
    try:
        qc.analyze(long_pairs, qc.parse_candidates(json.dumps({"candidates": [
            cand("long", 1, "chr1", 1, "+", 1),
            cand("long", 2, "chr1", 10, "-", 1)]}).encode()))
        check("limit: 501-base read rejected", False)
    except qc.QcInputError as exc:
        check("limit: 501-base read rejected", "500" in str(exc), str(exc))


def test_fastq_parsing() -> None:
    good = make_fastq([("a", "ACGT"), ("b", "ACGT")], 1)
    recs = parse_fastq(good, 1)
    check("fastq: ids parsed and /1 suffix stripped",
          [r.read_id for r in recs] == ["a", "b"])
    check("fastq: sequences retained", recs[0].sequence == b"ACGT")

    for name, blob in [
        ("empty file", b""),
        ("truncated tail", b"@x/1\nACGT\n+\n~~~~\n@y/1\nACGT\n"),
        ("bad header", b"x/1\nACGT\n+\n~~~~\n"),
        ("bad plus line", b"@x/1\nACGT\n-\n~~~~\n"),
        ("seq/qual length mismatch", b"@x/1\nACGT\n+\n~~\n"),
        ("empty sequence", b"@x/1\n\n+\n~\n"),
    ]:
        try:
            parse_fastq(blob, 1)
            check(f"fastq: {name} rejected", False)
        except Exception as exc:
            from app.fastq import FastqFormatError
            check(f"fastq: {name} rejected", isinstance(exc, FastqFormatError),
                  str(exc)[:60])

    # Pair-file mismatches.
    r1 = make_fastq([("a", "ACGT"), ("b", "ACGT")], 1)
    r2_missing = make_fastq([("a", "ACGT")], 2)
    try:
        validate_paired_files(parse_fastq(r1, 1), parse_fastq(r2_missing, 2))
        check("fastq: different record counts rejected", False)
    except Exception as exc:
        check("fastq: different record counts rejected", True, str(exc)[:60])

    r2_order = make_fastq([("b", "ACGT"), ("a", "ACGT")], 2)
    try:
        validate_paired_files(parse_fastq(r1, 1), parse_fastq(r2_order, 2))
        check("fastq: id order mismatch rejected", False)
    except Exception:
        check("fastq: id order mismatch rejected", True)

    r2_dup = make_fastq([("a", "ACGT"), ("a", "ACGT")], 2)
    try:
        validate_paired_files(parse_fastq(r1, 1), parse_fastq(r2_dup, 2))
        check("fastq: duplicate id within side rejected", False)
    except Exception:
        check("fastq: duplicate id within side rejected", True)


# --------------------------------------------------------------------------- #
# 2. HTTP end-to-end tests
# --------------------------------------------------------------------------- #

def post_analyze(r1: bytes, r2: bytes, cands: bytes,
                 task_id: str = "task-123"):
    return client.post(
        "/qc/analyze",
        data={"task_id": task_id},
        files={
            "r1": ("R1.fastq", r1, "application/octet-stream"),
            "r2": ("R2.fastq", r2, "application/octet-stream"),
            "candidates": ("candidates.json", cands, "application/json"),
        },
    )


def test_http_happy_path() -> None:
    pairs, candidates = make_clean_dataset(25, read_len=100, gap=100)
    ids = [p.read_id for p in pairs]
    r1 = make_fastq([(rid, "A" * 100) for rid in ids], 1)
    r2 = make_fastq([(rid, "C" * 100) for rid in ids], 2)
    resp = post_analyze(r1, r2, json.dumps({"candidates": candidates}).encode())
    check("http: happy path returns 200", resp.status_code == 200,
          str(resp.status_code) + resp.text[:200])
    body = resp.json()
    check("http: report shape and counts",
          body["pair_count"] == 25 and body["counts"]["normal"] == 25
          and body["calibration"]["median_span"] == 300.0)
    check("http: anomalies list empty and all five count keys present",
          body["anomalies"] == []
          and set(body["counts"]) == set(qc.ALL_CATEGORIES))


def test_http_full_failure_no_partial_results() -> None:
    pairs, candidates = make_clean_dataset(20)
    ids = [p.read_id for p in pairs]
    r1 = make_fastq([(rid, "A" * 50) for rid in ids], 1)
    r2 = make_fastq([(rid, "C" * 50) for rid in ids], 2)

    # Corrupted pair file -> 400 and an error envelope, never a report.
    bad_r1 = r1 + b"@tail/1\nACGT\n"
    resp = post_analyze(bad_r1, r2, json.dumps(
        {"candidates": candidates}).encode())
    check("http: corrupted pair file -> 400", resp.status_code == 400,
          str(resp.status_code))
    check("http: error envelope without result fields",
          "error" in resp.json() and "counts" not in resp.json())

    # Candidate with non-positive start -> 400.
    broken = [dict(c) for c in candidates]
    broken[0]["start"] = -3
    resp = post_analyze(r1, r2, json.dumps({"candidates": broken}).encode())
    check("http: invalid start -> 400 whole-request error",
          resp.status_code == 400 and "counts" not in resp.json())

    # Calibration failure -> 400.
    few_r1 = make_fastq([(rid, "A" * 50) for rid in ids[:19]], 1)
    few_r2 = make_fastq([(rid, "C" * 50) for rid in ids[:19]], 2)
    few_c = [c for c in candidates if c["read_id"] in set(ids[:19])]
    resp = post_analyze(few_r1, few_r2, json.dumps(
        {"candidates": few_c}).encode())
    check("http: calibration failure -> 400, no partial result",
          resp.status_code == 400
          and "calibration failed" in resp.json()["detail"])


def test_http_missing_field_and_empty_file() -> None:
    r1 = make_fastq([("a", "ACGT")], 1)
    r2 = make_fastq([("a", "ACGT")], 2)
    resp = client.post(
        "/qc/analyze",
        data={"task_id": "t"},
        files={"r1": ("R1.fastq", r1, "application/octet-stream"),
               "r2": ("R2.fastq", r2, "application/octet-stream")},
    )
    check("http: missing candidates field -> 422", resp.status_code == 422)

    resp = post_analyze(r1, r2, b"")
    check("http: empty candidates upload -> 400", resp.status_code == 400,
          str(resp.status_code))
    check("http: empty file detail message",
          "empty" in resp.json()["detail"]["detail"])


def test_http_health_and_docs() -> None:
    check("http: /health 200", client.get("/health").status_code == 200)
    check("http: /docs available", client.get("/docs").status_code == 200)
    openapi = client.get("/openapi.json").json()
    check("http: endpoint registered in openapi schema",
          "/qc/analyze" in openapi["paths"])


def main() -> int:
    tests = [
        test_calibration_median_mad_tolerance,
        test_span_anomaly_inside_and_outside_tolerance,
        test_mad_based_tolerance,
        test_calibration_needs_20_pairs,
        test_cross_reference_and_orientation,
        test_best_combined_score_wins,
        test_ties_consistent_verdict_adopted,
        test_ties_conflicting_verdict_undetermined,
        test_tie_normal_vs_span_is_undetermined,
        test_input_validation_errors,
        test_pair_and_length_limits,
        test_fastq_parsing,
        test_http_happy_path,
        test_http_full_failure_no_partial_results,
        test_http_missing_field_and_empty_file,
        test_http_health_and_docs,
    ]
    for test in tests:
        try:
            test()
        except Exception as exc:  # noqa: BLE001 - report every failure
            record(test.__name__, False, f"unexpected exception: {exc!r}")

    passed = sum(1 for _, ok, _ in RESULTS if ok)
    print(f"\n{passed}/{len(RESULTS)} checks passed")
    if passed != len(RESULTS):
        print("FAILED:", [name for name, ok, _ in RESULTS if not ok])
        return 1
    print("ALL SELF-TESTS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
