"""FastAPI entry point: paired-end alignment candidate QC review.

The endpoint is synchronous and stateless. It consumes the two paired FASTQ
files and one candidate JSON document, validates everything up front, and
either returns the complete QC report or fails the whole request — partial
results are never published and the uploaded originals are never modified.
"""

from __future__ import annotations

from fastapi import FastAPI, File, Form, HTTPException, UploadFile, status
from fastapi.responses import JSONResponse
from pydantic import BaseModel, Field

from . import __version__, config, qc
from .fastq import FastqFormatError, parse_fastq, validate_paired_files
from .qc import QcInputError

#: Hard ceiling for one whole multipart upload (all four fields combined).
#: At 2000 pairs x 500 bases per mate plus <= 12000 candidate entries the
#: real payload is far smaller; this only guards against abusive uploads.
MAX_UPLOAD_BYTES = 64 * 1024 * 1024


class CalibrationOut(BaseModel):
    pair_count: int = Field(description="Number of clean inward pairs used")
    median_span: float = Field(description="Median fragment span, reference bases")
    mad_span: float = Field(
        description="Median absolute deviation of the fragment span"
    )
    span_tolerance: float = Field(
        description="max(50, 4 x MAD); span deviations beyond this are abnormal"
    )


class AnomalyOut(BaseModel):
    position: int = Field(description="1-based position of the pair in R1 order")
    read_id: str
    category: str = Field(
        description="cross_reference | orientation | span | undetermined"
    )


class QcReport(BaseModel):
    task_id: str = Field(description="Task number echoed from the request")
    pair_count: int
    calibration: CalibrationOut
    counts: dict[str, int] = Field(
        description="Per-category counts: normal, cross_reference, orientation, "
        "span, undetermined"
    )
    anomalies: list[AnomalyOut]


app = FastAPI(
    title="Paired-end Fragment QC Review Service",
    version=__version__,
    description=(
        "Review anomalous paired-end fragments from a successfully paired task. "
        "The request carries the two paired FASTQ files and alignment candidates "
        f"for every read pair (at most {config.MAX_CANDIDATES_PER_END} per end, "
        f"each with reference name, 1-based start, strand and integer score; "
        "full-length gap-free alignments). The service calibrates the expected "
        "fragment span on single-candidate inward pairs and classifies every "
        "pair as normal, cross_reference, orientation, span or undetermined. "
        "Any input defect fails the whole request with no partial result."
    ),
)


@app.exception_handler(QcInputError)
async def qc_input_error_handler(_, exc: QcInputError) -> JSONResponse:
    return JSONResponse(
        status_code=status.HTTP_400_BAD_REQUEST,
        content={"error": "invalid_request", "detail": str(exc)},
    )


@app.exception_handler(FastqFormatError)
async def fastq_error_handler(_, exc: FastqFormatError) -> JSONResponse:
    return JSONResponse(
        status_code=status.HTTP_400_BAD_REQUEST,
        content={"error": "invalid_pair_files", "detail": str(exc)},
    )


async def _read_bounded(upload: UploadFile, label: str) -> bytes:
    data = await upload.read()
    if not data:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail={"error": "missing_file", "detail": f"{label} is empty"},
        )
    if len(data) > MAX_UPLOAD_BYTES:
        raise HTTPException(
            status_code=status.HTTP_413_REQUEST_ENTITY_TOO_LARGE,
            detail={
                "error": "upload_too_large",
                "detail": f"{label} exceeds the {MAX_UPLOAD_BYTES}-byte upload cap",
            },
        )
    return data


@app.post(
    "/qc/analyze",
    response_model=QcReport,
    responses={
        400: {"description": "Invalid pair files or candidate document (whole "
                             "request rejected, no partial result)"},
        413: {"description": "Upload exceeds the hard size cap"},
        422: {"description": "Missing required multipart field"},
    },
    summary="Classify paired-end fragments for structural-variant review",
)
async def analyze_pairs(
    task_id: str = Form(
        ..., description="Task number of the successfully paired upstream task"
    ),
    r1: UploadFile = File(..., description="Paired R1 FASTQ (uncompressed)"),
    r2: UploadFile = File(..., description="Paired R2 FASTQ (uncompressed)"),
    candidates: UploadFile = File(
        ...,
        description=(
            "JSON document: {\"candidates\": [{read_id, side (1|2), reference, "
            "start (1-based, >=1), strand ('+'|'-'), score (int)}, ...]} "
            "covering every read pair, at most 3 candidates per end"
        ),
    ),
) -> QcReport:
    r1_data = await _read_bounded(r1, "r1")
    r2_data = await _read_bounded(r2, "r2")
    candidate_data = await _read_bounded(candidates, "candidates")
    # The uploads are only read; nothing is ever written back to them.

    r1_records = parse_fastq(r1_data, side=1)
    r2_records = parse_fastq(r2_data, side=2)
    pair_tuples = validate_paired_files(r1_records, r2_records)
    pairs = [
        qc.PairInfo(read_id=read_id, r1_length=len1, r2_length=len2)
        for read_id, len1, len2 in pair_tuples
    ]

    candidate_list = qc.parse_candidates(candidate_data)
    report = qc.analyze(pairs, candidate_list)
    return QcReport(task_id=task_id, **report)


@app.get("/health", summary="Liveness probe")
def health() -> dict[str, str]:
    return {"status": "ok"}
