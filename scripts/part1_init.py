#!/usr/bin/env python3
"""Create a Part 1 stage directory. Does not run QC, Harmony, or Leiden.

Part 1 clustering is executed from references/part1-qc-and-global-atlas.md
with Scanpy. After a human-locked clustering object exists, use
cluster_review_tables.py for marker worksheets.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from stagecraft import EXIT_OK, EXIT_USAGE  # noqa: E402
from stagecraft.io import require_new  # noqa: E402

STAGE_DIRS = (
    "00_protocol_manifest",
    "01_code",
    "02_tables",
    "03_figures",
    "05_logs",
    "06_reports",
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--gene", required=True, help="Pre-specified TARGET_GENE symbol")
    args = parser.parse_args(argv)
    root = args.out
    if root.exists():
        raise SystemExit(f"Refusing to overwrite: {root}")
    for name in STAGE_DIRS:
        (root / name).mkdir(parents=True)
    template = (_ROOT / "references" / "protocol-template.md").read_text(encoding="utf-8")
    protocol = template.replace("`TARGET_GENE`: SYMBOL", f"`TARGET_GENE`: {args.gene}")
    protocol_path = require_new(root / "00_protocol_manifest" / "PROTOCOL.md")
    protocol_path.write_text(protocol, encoding="utf-8")
    readme = require_new(root / "README.md")
    readme.write_text(
        "\n".join(
            [
                f"# Part 1 stage for `{args.gene}`",
                "",
                "This directory was created by `scripts/part1_init.py`.",
                "It is not a completed atlas.",
                "",
                "Execute QC, strict-common merge, Harmony and Leiden from",
                "`references/part1-qc-and-global-atlas.md`. Human gates remain",
                "required for resolution choice and KEEP/DELETE naming.",
                "",
                "After clustering exists, export review tables with",
                "`python scripts/cluster_review_tables.py`.",
                "",
            ]
        ),
        encoding="utf-8",
    )
    print(f"wrote Part 1 stage scaffold: {root}")
    print("Next: freeze PROTOCOL.md, then execute the Part 1 spec. No Scanpy runner is bundled.")
    return EXIT_OK if args.gene.strip() else EXIT_USAGE


if __name__ == "__main__":
    raise SystemExit(main())
