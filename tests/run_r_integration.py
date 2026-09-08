"""Execute the Part 5 audit boundaries and real-model integration tests."""
from __future__ import annotations

import argparse
import os
import pathlib
import shutil
import subprocess


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rscript", default=os.environ.get("STAGECRAFT_RSCRIPT") or shutil.which("Rscript"))
    args = parser.parse_args()
    root = pathlib.Path(__file__).parents[1]
    if not args.rscript:
        parser.error("Set STAGECRAFT_RSCRIPT, pass --rscript, or add Rscript to PATH")
    for test in ("audit_boundaries.R", "scientific_regression.R"):
        completed = subprocess.run([args.rscript, "--vanilla", str(root / "tests" / test)],
                                   cwd=root, check=False, timeout=600)
        if completed.returncode:
            return completed.returncode
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

