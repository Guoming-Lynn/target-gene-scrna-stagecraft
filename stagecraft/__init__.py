"""Shared helpers for target-gene-scrna-stagecraft CLI scripts."""

from __future__ import annotations

__version__ = "2.2.0"

# CLI exit codes. 0 success; 1 usage/validation; 2 gate/environment failure;
# 3 scientific non-estimability or empty KO-eligible set.
EXIT_OK = 0
EXIT_USAGE = 1
EXIT_GATE = 2
EXIT_NOT_ESTIMABLE = 3
