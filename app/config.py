"""Runtime configuration for the paired-alignment QC review service.

All limits come straight from the delivery requirement and may be overridden
through environment variables (the defaults already match the specification).
"""

from __future__ import annotations

import os

#: Maximum number of read pairs accepted in one request.
MAX_PAIRS = int(os.environ.get("QC_MAX_PAIRS", 2000))

#: Maximum sequence length (bases) of either mate of a pair.
MAX_READ_LENGTH = int(os.environ.get("QC_MAX_READ_LENGTH", 500))

#: Maximum number of alignment candidates supplied for one end of one read.
MAX_CANDIDATES_PER_END = int(os.environ.get("QC_MAX_CANDIDATES_PER_END", 3))

#: Minimum number of inward, single-candidate calibration pairs required.
MIN_CALIBRATION_PAIRS = int(os.environ.get("QC_MIN_CALIBRATION_PAIRS", 20))

#: Absolute slack (bases) used together with the MAD-derived span tolerance.
SPAN_ABSOLUTE_TOLERANCE = int(os.environ.get("QC_SPAN_ABSOLUTE_TOLERANCE", 50))

#: Multiplier applied to the MAD when building the span tolerance.
SPAN_MAD_MULTIPLIER = int(os.environ.get("QC_SPAN_MAD_MULTIPLIER", 4))
