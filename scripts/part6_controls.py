#!/usr/bin/env python3
"""Match ≤10 negative-control genes to TARGET_GENE on raw detection + mean.

Controls are descriptive comparators. They are not a permutation null.
One window expansion is allowed. A second expansion is not.

Required columns on --genes:
    gene_symbol, detection_fraction, mean_raw_counts_per_cell, in_model_vocabulary

Usage:
    python scripts/part6_controls.py gene_stats.csv \\
        --target TARGET_GENE \\
        --endpoint-union endpoint_members.txt \\
        --out 05_controls/frozen_controls.csv
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

CYCLE = {
    "MKI67", "TOP2A", "PCNA", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7",
    "CDK1", "CCNA2", "CCNB1", "CCNB2", "UBE2C", "BIRC5", "CENPF", "TYMS",
}
TECHNICAL = {
    "ACTB", "GAPDH", "B2M", "MALAT1", "XIST", "TMSB10", "TMSB4X",
    "EEF1A1",
}
FAMILY = re.compile(
    r"^(MT-|RPS|RPL|MRPS|MRPL|JUN|FOS|STAT|IRF|NFKB|CEBP|KLF|ZNF|GATA|"
    r"SPI|RUNX|MYC|EGR|ATF|CREB|FOXO|HIF|REL)"
)


def exclusion_reason(symbol: str, mean: float, low: float, high: float, target: str, endpoint_union: set[str]) -> str:
    if symbol == target:
        return "target"
    if symbol in endpoint_union:
        return "frozen_endpoint_union"
    if symbol in CYCLE:
        return "cell_cycle_core"
    if symbol in TECHNICAL:
        return "technical_extreme"
    if FAMILY.match(symbol):
        return "excluded_family"
    if mean <= low or mean >= high:
        return "expression_extreme"
    return ""


def logit(p: np.ndarray) -> np.ndarray:
    clipped = np.clip(p.astype(float), 1e-6, 1 - 1e-6)
    return np.log(clipped / (1 - clipped))


def select_controls(
    genes: pd.DataFrame,
    *,
    target: str,
    endpoint_union: set[str],
    excluded_genes: set[str] | None = None,
    max_n: int = 10,
    window_detection: float = 0.05,
    window_log2_mean: float = 0.5,
    expand_detection: float = 0.10,
    expand_log2_mean: float = 1.0,
) -> tuple[pd.DataFrame, pd.DataFrame, str]:
    required = {"gene_symbol", "detection_fraction", "mean_raw_counts_per_cell", "in_model_vocabulary"}
    missing = required.difference(genes.columns)
    if missing:
        raise SystemExit(f"gene stats missing columns: {sorted(missing)}")
    frame = genes.loc[genes["in_model_vocabulary"].astype(bool)].copy()
    frame = frame.loc[~frame["gene_symbol"].duplicated(keep=False)].copy()
    if frame["gene_symbol"].eq(target).sum() != 1:
        raise SystemExit(f"{target} must have exactly one model-visible row")
    tgt = frame.loc[frame["gene_symbol"].eq(target)].iloc[0]
    low, high = frame["mean_raw_counts_per_cell"].quantile([0.01, 0.99])
    frame["exclusion_reason"] = [
        exclusion_reason(str(sym), float(mean), float(low), float(high), target, endpoint_union)
        for sym, mean in zip(frame["gene_symbol"], frame["mean_raw_counts_per_cell"])
    ]
    frame.loc[frame["gene_symbol"].isin(excluded_genes or set()), "exclusion_reason"] = "frozen_project_exclusion"
    eligible = frame.loc[frame["exclusion_reason"].eq("")].copy()

    def in_window(pool: pd.DataFrame, det: float, log2m: float) -> pd.DataFrame:
        mean_ratio = np.log2(
            (pool["mean_raw_counts_per_cell"].to_numpy(float) + 1e-8)
            / (float(tgt.mean_raw_counts_per_cell) + 1e-8)
        )
        det_delta = (pool["detection_fraction"].to_numpy(float) - float(tgt.detection_fraction))
        mask = np.abs(det_delta) <= det
        mask &= np.abs(mean_ratio) <= log2m
        return pool.loc[mask].copy()

    window = "primary"
    selected_pool = in_window(eligible, window_detection, window_log2_mean)
    if len(selected_pool) < max_n:
        selected_pool = in_window(eligible, expand_detection, expand_log2_mean)
        window = "expanded_once"

    x1 = logit(eligible["detection_fraction"].to_numpy(float))
    x2 = np.log1p(eligible["mean_raw_counts_per_cell"].to_numpy(float))
    z1 = (x1 - x1.mean()) / (x1.std() if x1.std() else 1.0)
    z2 = (x2 - x2.mean()) / (x2.std() if x2.std() else 1.0)
    t1 = logit(np.array([float(tgt.detection_fraction)]))[0]
    t2 = np.log1p(float(tgt.mean_raw_counts_per_cell))
    tz1 = (t1 - x1.mean()) / (x1.std() if x1.std() else 1.0)
    tz2 = (t2 - x2.mean()) / (x2.std() if x2.std() else 1.0)
    eligible = eligible.copy()
    eligible["distance"] = np.sqrt((z1 - tz1) ** 2 + (z2 - tz2) ** 2)
    selected_pool = selected_pool.merge(
        eligible[["gene_symbol", "distance"]], on="gene_symbol", how="left", validate="one_to_one"
    )
    selected = selected_pool.sort_values(["distance", "gene_symbol"], kind="stable").head(max_n).copy()
    selected["control_rank"] = np.arange(1, len(selected) + 1)
    selected["matching_window"] = window
    selected["target_gene"] = target
    return selected, frame.merge(eligible[["gene_symbol", "distance"]], on="gene_symbol", how="left"), window


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("genes", type=Path)
    parser.add_argument("--target", required=True)
    parser.add_argument("--endpoint-union", type=Path, required=True)
    parser.add_argument("--exclude-genes", type=Path, help="Frozen project-specific exclusions, one symbol per line")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--max-n", type=int, default=10)
    args = parser.parse_args(argv)
    union = {line.strip() for line in args.endpoint_union.read_text(encoding="utf-8").splitlines() if line.strip()}
    excluded = {line.strip() for line in args.exclude_genes.read_text(encoding="utf-8").splitlines() if line.strip()} if args.exclude_genes else set()
    selected, candidates, window = select_controls(
        pd.read_csv(args.genes), target=args.target, endpoint_union=union, excluded_genes=excluded, max_n=args.max_n
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    if args.out.exists():
        raise SystemExit(f"Refusing to overwrite: {args.out}")
    selected.to_csv(args.out, index=False)
    candidates.to_csv(args.out.with_name("control_candidates.csv"), index=False)
    print(f"n_controls={len(selected)} window={window}")
    if len(selected) < 5:
        print("control ranks are NOT_ESTIMABLE (<5)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

