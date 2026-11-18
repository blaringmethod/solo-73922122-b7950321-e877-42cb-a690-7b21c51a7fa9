"""Core paired-end fragment QC logic (pure functions, no I/O framework).

Pipeline:

1. Validate candidate coverage against the paired FASTQ files.
2. Calibrate the expected fragment span from the "clean" subset: pairs with a
   single candidate on each end, both on the same reference, R1 forward on the
   left and R2 reverse on the right (inward facing). At least 20 such pairs are
   required, otherwise calibration fails and the whole request errors out.
3. For every pair, enumerate every R1 x R2 candidate combination, keep the
   combination(s) with the highest combined score and classify them:
   different references -> cross-reference; same reference but not inward ->
   orientation; inward -> normal/span by the calibrated tolerance. Tied
   combinations only count when they all reach the same conclusion, otherwise
   the pair is undetermined.
"""

from __future__ import annotations

import json
import statistics
from dataclasses import dataclass
from typing import Any

from . import config

CROSS_REFERENCE = "cross_reference"
ORIENTATION = "orientation"
SPAN = "span"
NORMAL = "normal"
UNDETERMINED = "undetermined"

ABNORMAL_CATEGORIES = (CROSS_REFERENCE, ORIENTATION, SPAN)
ALL_CATEGORIES = (NORMAL, CROSS_REFERENCE, ORIENTATION, SPAN, UNDETERMINED)

FORWARD = "+"
REVERSE = "-"


class QcInputError(ValueError):
    """A request-level input error; no partial result may be returned."""


class CalibrationError(QcInputError):
    """Raised when fewer than 20 clean pairs are available for calibration."""


@dataclass(frozen=True)
class Candidate:
    """One full-length, gap-free alignment candidate for one end of a read."""

    read_id: str
    side: int
    reference: str
    start: int  # 1-based, inclusive
    strand: str  # "+" forward or "-" reverse
    score: int


@dataclass(frozen=True)
class PairInfo:
    """One read pair as derived from the two paired FASTQ files."""

    read_id: str
    r1_length: int
    r2_length: int


def _is_int(value: Any) -> bool:
    # ``bool`` is a subclass of int; JSON true/false must not pass as numbers.
    return isinstance(value, int) and not isinstance(value, bool)


def parse_candidates(raw: bytes) -> list[Candidate]:
    """Parse and structurally validate the candidate JSON document.

    Accepted shape: ``{"candidates": [ {read_id, side, reference, start,
    strand, score}, ... ]}``. Any structural problem raises QcInputError.
    """
    try:
        document = json.loads(raw.decode("utf-8"))
    except UnicodeDecodeError as exc:
        raise QcInputError(f"candidate file is not valid UTF-8: {exc}") from exc
    except json.JSONDecodeError as exc:
        raise QcInputError(f"candidate file is not valid JSON: {exc}") from exc

    if not isinstance(document, dict) or "candidates" not in document:
        raise QcInputError(
            "candidate document must be a JSON object with a 'candidates' array"
        )
    entries = document["candidates"]
    if not isinstance(entries, list):
        raise QcInputError("'candidates' must be an array")

    parsed: list[Candidate] = []
    seen: set[tuple[Any, ...]] = set()
    for index, entry in enumerate(entries):
        where = f"candidate #{index + 1}"
        if not isinstance(entry, dict):
            raise QcInputError(f"{where}: entry must be a JSON object")
        missing = [
            key
            for key in ("read_id", "side", "reference", "start", "strand", "score")
            if key not in entry
        ]
        if missing:
            raise QcInputError(f"{where}: missing field(s) {', '.join(missing)}")

        read_id = entry["read_id"]
        side = entry["side"]
        reference = entry["reference"]
        start = entry["start"]
        strand = entry["strand"]
        score = entry["score"]

        if not isinstance(read_id, str) or not read_id:
            raise QcInputError(f"{where}: 'read_id' must be a non-empty string")
        if not _is_int(side) or side not in (1, 2):
            raise QcInputError(f"{where}: 'side' must be the integer 1 or 2")
        if not isinstance(reference, str) or not reference:
            raise QcInputError(
                f"{where}: 'reference' must be a non-empty reference name"
            )
        if not _is_int(start):
            raise QcInputError(f"{where}: 'start' must be an integer")
        if start < 1:
            raise QcInputError(
                f"{where}: start must be a positive 1-based coordinate, got {start}"
            )
        if strand not in (FORWARD, REVERSE):
            raise QcInputError(f"{where}: 'strand' must be '+' or '-'")
        if not _is_int(score):
            raise QcInputError(f"{where}: 'score' must be an integer")

        identity = (read_id, side, reference, start, strand)
        if identity in seen:
            raise QcInputError(
                f"{where}: duplicate candidate for read '{read_id}' R{side} "
                f"({reference}:{start} {strand})"
            )
        seen.add(identity)
        parsed.append(
            Candidate(
                read_id=read_id,
                side=side,
                reference=reference,
                start=start,
                strand=strand,
                score=score,
            )
        )
    if not parsed:
        raise QcInputError("candidate document contains no candidates")
    return parsed


def _index_candidates(
    candidates: list[Candidate],
) -> dict[tuple[str, int], list[Candidate]]:
    index: dict[tuple[str, int], list[Candidate]] = {}
    for cand in candidates:
        index.setdefault((cand.read_id, cand.side), []).append(cand)
    return index


def validate_request(
    pairs: list[PairInfo], candidates: list[Candidate]
) -> dict[tuple[str, int], list[Candidate]]:
    """Cross-validate pair files against the candidate set.

    Enforces the 2000-pair / 500-base / 3-candidate caps, complete coverage
    (no missing or unknown candidates) and unique candidate identities.
    """
    if len(pairs) > config.MAX_PAIRS:
        raise QcInputError(
            f"{len(pairs)} read pairs exceed the limit of {config.MAX_PAIRS}"
        )

    index = _index_candidates(candidates)
    pair_ids = {pair.read_id for pair in pairs}

    # Candidates referencing reads absent from the pair files are also a
    # coverage failure (unknown / dangling candidate id).
    for (read_id, side) in index:
        if read_id not in pair_ids:
            raise QcInputError(
                f"candidate references unknown read id '{read_id}' (R{side}); "
                "it is absent from the paired files"
            )

    for pair in pairs:
        if pair.r1_length > config.MAX_READ_LENGTH:
            raise QcInputError(
                f"read '{pair.read_id}' R1 length {pair.r1_length} exceeds the "
                f"{config.MAX_READ_LENGTH}-base limit"
            )
        if pair.r2_length > config.MAX_READ_LENGTH:
            raise QcInputError(
                f"read '{pair.read_id}' R2 length {pair.r2_length} exceeds the "
                f"{config.MAX_READ_LENGTH}-base limit"
            )
        for side in (1, 2):
            end_candidates = index.get((pair.read_id, side))
            if not end_candidates:
                raise QcInputError(
                    f"missing alignment candidate(s) for read '{pair.read_id}' "
                    f"R{side}"
                )
            if len(end_candidates) > config.MAX_CANDIDATES_PER_END:
                raise QcInputError(
                    f"read '{pair.read_id}' R{side} has {len(end_candidates)} "
                    f"candidates; at most {config.MAX_CANDIDATES_PER_END} are allowed"
                )
    return index


def _is_clean_inward(
    r1: list[Candidate],
    r2: list[Candidate],
    pair: PairInfo,
) -> int | None:
    """Return the fragment span for a calibration pair, else None.

    Calibration pairs carry exactly one candidate per end, on the same
    reference, with R1 forward on the left and R2 reverse on the right.
    """
    if len(r1) != 1 or len(r2) != 1:
        return None
    c1, c2 = r1[0], r2[0]
    if c1.reference != c2.reference:
        return None
    if not (c1.strand == FORWARD and c2.strand == REVERSE and c1.start < c2.start):
        return None
    return _fragment_span(c1, c2, pair)


def _fragment_span(c1: Candidate, c2: Candidate, pair: PairInfo) -> int:
    """Span of an inward full-length alignment, in reference bases.

    R1 starts at ``c1.start`` and extends right by its read length; R2 starts
    (its 5' end, 1-based) at ``c2.start`` on the reverse strand and extends
    left, so its rightmost aligned base is ``c2.start + r2_length - 1``.
    """
    rightmost = c2.start + pair.r2_length - 1
    return rightmost - c1.start + 1


def calibrate(
    pairs: list[PairInfo], index: dict[tuple[str, int], list[Candidate]]
) -> tuple[float, float, float, int]:
    """Compute median span, MAD and tolerance from the clean inward subset.

    Returns ``(median_span, mad_span, tolerance, calibration_pair_count)``.
    Raises CalibrationError when fewer than ``MIN_CALIBRATION_PAIRS`` clean
    pairs exist.
    """
    spans: list[int] = []
    for pair in pairs:
        span = _is_clean_inward(
            index[(pair.read_id, 1)], index[(pair.read_id, 2)], pair
        )
        if span is not None:
            spans.append(span)

    if len(spans) < config.MIN_CALIBRATION_PAIRS:
        raise CalibrationError(
            "calibration failed: only "
            f"{len(spans)} single-candidate inward pair(s) on one reference; "
            f"at least {config.MIN_CALIBRATION_PAIRS} are required"
        )

    median_span = statistics.median(spans)
    mad_span = statistics.median(abs(span - median_span) for span in spans)
    tolerance = max(
        float(config.SPAN_ABSOLUTE_TOLERANCE),
        config.SPAN_MAD_MULTIPLIER * mad_span,
    )
    return median_span, mad_span, tolerance, len(spans)


def _classify_combination(
    c1: Candidate,
    c2: Candidate,
    pair: PairInfo,
    median_span: float,
    tolerance: float,
) -> str:
    if c1.reference != c2.reference:
        return CROSS_REFERENCE
    inward = (
        c1.strand == FORWARD
        and c2.strand == REVERSE
        and c1.start < c2.start
    )
    if not inward:
        return ORIENTATION
    span = _fragment_span(c1, c2, pair)
    if abs(span - median_span) <= tolerance:
        return NORMAL
    return SPAN


def analyze(
    pairs: list[PairInfo], candidates: list[Candidate]
) -> dict[str, Any]:
    """Run the full QC analysis and return a JSON-serializable report dict."""
    index = validate_request(pairs, candidates)
    median_span, mad_span, tolerance, calibration_count = calibrate(pairs, index)

    counts = {category: 0 for category in ALL_CATEGORIES}
    anomalies: list[dict[str, Any]] = []

    for position, pair in enumerate(pairs, start=1):
        r1_candidates = index[(pair.read_id, 1)]
        r2_candidates = index[(pair.read_id, 2)]
        combinations = [
            (c1, c2, c1.score + c2.score)
            for c1 in r1_candidates
            for c2 in r2_candidates
        ]
        best_score = max(combined[2] for combined in combinations)
        tied = [combined for combined in combinations if combined[2] == best_score]
        verdicts = {
            _classify_combination(c1, c2, pair, median_span, tolerance)
            for c1, c2, _ in tied
        }
        category = next(iter(verdicts)) if len(verdicts) == 1 else UNDETERMINED
        counts[category] += 1
        if category in ABNORMAL_CATEGORIES or category == UNDETERMINED:
            anomalies.append(
                {"position": position, "read_id": pair.read_id, "category": category}
            )

    return {
        "pair_count": len(pairs),
        "calibration": {
            "pair_count": calibration_count,
            "median_span": median_span,
            "mad_span": mad_span,
            "span_tolerance": tolerance,
        },
        "counts": counts,
        "anomalies": anomalies,
    }
