"""Shared helpers for target-gene-scrna-stagecraft CLI scripts."""

from __future__ import annotations

import sys

__version__ = "2.3.0"

# CLI exit codes. 0 success; 1 usage/validation; 2 gate/environment failure;
# 3 scientific non-estimability or empty KO-eligible set.
EXIT_OK = 0
EXIT_USAGE = 1
EXIT_GATE = 2
EXIT_NOT_ESTIMABLE = 3


def stop(message: str, code: int = EXIT_USAGE) -> None:
    """Print ``message`` and exit. Gates pass ``EXIT_GATE``; bad input stays 1."""
    print(message, file=sys.stderr)
    raise SystemExit(code)
