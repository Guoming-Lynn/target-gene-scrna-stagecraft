#!/usr/bin/env python3
"""Explicitly authorized cell-level exploratory fallback."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
from scipy.stats import mannwhitneyu

AUTH = "I_ACCEPT_CELL_LEVEL_FALSE_POSITIVE_RISK"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("table", type=Path)
    parser.add_argument("--group", required=True)
    parser.add_argument("--value", required=True)
    parser.add_argument("--group-a", required=True)
    parser.add_argument("--group-b", required=True)
    parser.add_argument("--authorization", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.authorization != AUTH:
        raise SystemExit(f"Exact authorization required: {AUTH}")
    frame = pd.read_csv(args.table)
    required = {args.group, args.value}
    if not required.issubset(frame.columns):
        raise SystemExit(f"missing columns: {sorted(required - set(frame.columns))}")
    group = frame[args.group].astype(str)
    values = pd.to_numeric(frame[args.value], errors="coerce")
    left = values.loc[group.eq(args.group_a)].dropna()
    right = values.loc[group.eq(args.group_b)].dropna()
    if len(left) < 2 or len(right) < 2:
        raise SystemExit("both groups need at least two finite cells")
    out = {
        "analysis": "CELL_LEVEL_EXPLORATORY_FALLBACK",
        "status": "EXPLORATORY_ONLY",
        "false_positive_risk_accepted": True,
        "authorization": AUTH,
        "group": args.group,
        "value": args.value,
        "group_a": args.group_a,
        "group_b": args.group_b,
        "n_cells_a": len(left),
        "n_cells_b": len(right),
        "test": "Mann-Whitney U; cells are dependent subsamples",
        "p_value": float(mannwhitneyu(left, right, alternative="two-sided").pvalue),
        "warning": "Not a donor-level finding; not replication; not a Part 5 verdict.",
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    if args.out.exists():
        raise SystemExit(f"Refusing to overwrite: {args.out}")
    args.out.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(f"wrote exploratory cell-level result: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
