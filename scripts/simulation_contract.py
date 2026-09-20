#!/usr/bin/env python3
"""Create a frozen calibration manifest; this command does not run simulations."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from stagecraft.io import require_new  # noqa: E402

SCENARIOS = (
    "NULL",
    "DONOR_EFFECT",
    "SOURCE_CONFOUNDED",
    "DEPTH_COLLINEAR",
    "ONE_DONOR_DRIVEN",
    "QC_ATTRITION",
    "FAILED_HOLDOUT",
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument(
        "--declared-replicates",
        "--replicates",
        dest="replicates",
        type=int,
        default=1000,
        help="Declared replicate count for a later stage run. This CLI writes a manifest only.",
    )
    parser.add_argument("--scenario", choices=SCENARIOS, action="append")
    args = parser.parse_args(argv)
    if args.replicates < 100:
        raise SystemExit("Declare at least 100 replicates")
    payload = {
        "protocol": "target-gene-scrna-stagecraft simulation calibration",
        "seed": args.seed,
        "replicates": args.replicates,
        "scenarios": args.scenario or list(SCENARIOS),
        "required_outputs": [
            "false_positive_rate",
            "FDR",
            "power",
            "CI_coverage",
            "NOT_ESTIMABLE_rate",
            "LODO_failure_rate",
        ],
        "status": "MANIFEST_ONLY",
        "note": (
            "This command does not sample data or fit models. "
            "Run scripts/calibrate_part5_null.R for executable null draws. "
            "Generic OLS is not a substitute."
        ),
    }
    path = require_new(args.out)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"wrote simulation manifest (MANIFEST_ONLY): {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
