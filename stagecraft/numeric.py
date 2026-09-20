"""Numeric guards for variance and rank checks. Not scientific thresholds."""

from __future__ import annotations

import math


def nonpositive(value: float, atol: float = 1e-15) -> bool:
    return (not math.isfinite(value)) or value <= atol
