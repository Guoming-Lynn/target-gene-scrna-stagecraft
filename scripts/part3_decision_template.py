#!/usr/bin/env python3
"""Blank KEEP/DELETE worksheet for a Part 3 STOP POINT 2 round.

Does not name clusters. Empty decision cells are intentional: the human
must fill KEEP or DELETE for every cluster before prepare-removal.

Usage:
    python scripts/part3_decision_template.py cluster_qc.csv \\
        --round round_00_initial --leiden-key leiden_r0_5 \\
        --out tables/round_00_initial/manual_decision_template.csv
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def _cluster_column(frame: pd.DataFrame) -> str:
    for name in ("cluster_id", "cluster", "leiden"):
        if name in frame.columns:
            return name
    raise SystemExit("QC table needs a cluster / cluster_id column.")


def build_template(
    qc: pd.DataFrame,
    *,
    round_name: str,
    leiden_key: str,
) -> pd.DataFrame:
    cluster_col = _cluster_column(qc)
    frame = qc.copy()
    frame[cluster_col] = frame[cluster_col].astype(str)
    n_cells = frame["n_cells"] if "n_cells" in frame.columns else pd.NA
    out = pd.DataFrame(
        {
            "round": round_name,
            "selected_leiden_column": leiden_key,
            "cluster_id": frame[cluster_col].to_numpy(),
            "n_cells": n_cells,
            "decision": "",
            "reason": "",
            "reviewer": "",
            "review_date": "",
        }
    )
    return out.sort_values(
        "cluster_id",
        key=lambda s: s.map(lambda x: int(x) if str(x).isdigit() else x),
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("qc_csv", type=Path)
    parser.add_argument("--round", required=True)
    parser.add_argument("--leiden-key", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)

    qc = pd.read_csv(args.qc_csv, dtype=str)
    if "n_cells" in qc.columns:
        qc["n_cells"] = pd.to_numeric(qc["n_cells"], errors="coerce")
    template = build_template(qc, round_name=args.round, leiden_key=args.leiden_key)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    if args.out.exists():
        raise SystemExit(f"Refusing to overwrite existing decision template: {args.out}")
    template.to_csv(args.out, index=False, encoding="utf-8-sig")
    print(f"wrote blank KEEP/DELETE template ({len(template)} clusters): {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

