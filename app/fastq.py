"""FASTQ parsing and pair-file validation.

A request is capped at 2000 pairs with at most 500 bases per mate, so the
whole payload is small; records are parsed into memory. The parser only ever
reads the uploaded bytes, never writes to or rewrites the original files.
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass


class FastqFormatError(ValueError):
    """Raised when a FASTQ file violates the four-line record contract."""


@dataclass(frozen=True)
class FastqRecord:
    """One parsed FASTQ record: canonical read id plus the raw sequence."""

    read_id: str
    sequence: bytes


def extract_read_id(header_line: bytes) -> str:
    """Return the canonical read id from a FASTQ header line.

    The id is the field before the first whitespace, without the leading
    ``@`` and without a trailing ``/1`` or ``/2``.
    """
    if not header_line.startswith(b"@"):
        raise FastqFormatError("header line does not start with '@'")
    first_field = header_line[1:].split(None, 1)[0] if len(header_line) > 1 else b""
    if not first_field:
        raise FastqFormatError("empty read id in header line")
    if len(first_field) >= 2 and first_field[-2:-1] == b"/" and first_field[-1:] in (
        b"1",
        b"2",
    ):
        first_field = first_field[:-2]
    if not first_field:
        raise FastqFormatError("empty read id after stripping '/1' or '/2'")
    # Latin-1 maps every byte 1:1 to a unicode code point, so no input byte
    # sequence can fail to decode and ids round-trip losslessly.
    return first_field.decode("latin-1")


def _iter_records(data: bytes, side_label: str) -> Iterator[FastqRecord]:
    if not data:
        raise FastqFormatError(f"{side_label}: file is empty")
    # splitlines() treats a final newline correctly and rejects a stray
    # trailing segment because the total line count must be divisible by 4.
    lines = data.splitlines()
    if len(lines) % 4 != 0:
        raise FastqFormatError(
            f"{side_label}: truncated trailing record (expected four lines per "
            f"record, found {len(lines) % 4} trailing line(s))"
        )
    record_no = 0
    for i in range(0, len(lines), 4):
        record_no += 1
        header, sequence, plus, quality = lines[i : i + 4]
        if not header.startswith(b"@"):
            raise FastqFormatError(
                f"{side_label}: record {record_no} header does not start with '@'"
            )
        if not plus.startswith(b"+"):
            raise FastqFormatError(
                f"{side_label}: record {record_no} third line does not start with '+'"
            )
        if not sequence:
            raise FastqFormatError(
                f"{side_label}: record {record_no} has an empty sequence"
            )
        if not quality:
            raise FastqFormatError(
                f"{side_label}: record {record_no} has an empty quality string"
            )
        if len(sequence) != len(quality):
            raise FastqFormatError(
                f"{side_label}: record {record_no} sequence length {len(sequence)} "
                f"!= quality length {len(quality)}"
            )
        yield FastqRecord(extract_read_id(header), sequence)


def parse_fastq(data: bytes, side: int) -> list[FastqRecord]:
    """Parse one FASTQ byte string into ordered records."""
    return list(_iter_records(data, f"R{side}"))


def validate_paired_files(
    r1_records: list[FastqRecord], r2_records: list[FastqRecord]
) -> list[tuple[str, int, int]]:
    """Validate the two pair files and return ``(id, len_R1, len_R2)`` pairs.

    Both files must contain the same read ids, each exactly once per side, in
    the same order (the re-pairing task upstream already aligns them). Any
    violation fails the whole request.
    """
    r1_ids = [r.read_id for r in r1_records]
    r2_ids = [r.read_id for r in r2_records]

    seen: set[str] = set()
    for rid in r1_ids:
        if rid in seen:
            raise FastqFormatError(f"R1: duplicate read id '{rid}'")
        seen.add(rid)
    seen.clear()
    for rid in r2_ids:
        if rid in seen:
            raise FastqFormatError(f"R2: duplicate read id '{rid}'")
        seen.add(rid)

    if len(r1_ids) != len(r2_ids):
        raise FastqFormatError(
            f"pair files hold different numbers of records: R1={len(r1_ids)}, "
            f"R2={len(r2_ids)}"
        )
    pairs: list[tuple[str, int, int]] = []
    for idx, (rec1, rec2) in enumerate(zip(r1_records, r2_records), start=1):
        if rec1.read_id != rec2.read_id:
            raise FastqFormatError(
                f"pair position {idx}: R1 id '{rec1.read_id}' does not match "
                f"R2 id '{rec2.read_id}'"
            )
        pairs.append((rec1.read_id, len(rec1.sequence), len(rec2.sequence)))
    if not pairs:
        raise FastqFormatError("pair files contain no records")
    return pairs
