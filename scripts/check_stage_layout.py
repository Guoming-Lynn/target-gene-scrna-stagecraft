#!/usr/bin/env python3
"""Check a stage directory against the layout contract.

Usage:
    python scripts/check_stage_layout.py path/to/analysis/NN_name
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

REQUIRED_DIRS = [
    "00_protocol_manifest",
    "01_code",
    "05_logs",
    "06_reports",
]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage_dir", type=Path)
    args = parser.parse_args(argv)
    root = args.stage_dir
    if not root.is_dir():
        print(f"missing dir: {root}", file=sys.stderr)
        return 1
    failures = []
    for name in REQUIRED_DIRS:
        if not (root / name).is_dir():
            failures.append(f"missing {name}/")
    proto = root / "00_protocol_manifest"
    if proto.is_dir():
        prots = list(proto.glob("PROTOCOL*.md")) + list(proto.glob("FROZEN_PROTOCOL*.md"))
        if not prots:
            failures.append("00_protocol_manifest has no PROTOCOL*.md")
        else:
            for protocol in prots:
                text = protocol.read_text(encoding="utf-8", errors="replace")
                if len(text.strip()) < 1200 or "NOT_ESTIMABLE" not in text:
                    failures.append(f"{protocol.name} is not a substantive validated protocol")
    logs = root / "05_logs"
    if logs.is_dir() and not (logs / "verdict.json").exists():
        reports = root / "06_reports"
        if reports.is_dir() and any(reports.glob("*.md")):
            failures.append("06_reports exists but 05_logs/verdict.json is missing")
    tables = (
        (root / "02_tables").is_dir()
        or (root / "03_tables").is_dir()
        or (root / "tables").is_dir()
    )
    if not tables:
        failures.append("missing 02_tables/ or 03_tables/ or tables/")
    report_files = list((root / "06_reports").glob("*.md")) if (root / "06_reports").is_dir() else []
    if report_files and not (logs / "verdict.json").exists():
        failures.append("reports cannot exist before 05_logs/verdict.json")
    if failures:
        print("FAIL", root)
        for item in failures:
            print(" -", item)
        return 1
    print("OK", root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
