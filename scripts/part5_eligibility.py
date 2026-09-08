#!/usr/bin/env python3
"""Eligibility of the design that will actually be fit.

Computes rank and residual df in Python before R. Does not fit limma.
Does not override a Part 4 LOW_N / NO_EXPOSURE_RANGE arm into formal.

Usage:
    python scripts/part5_eligibility.py 02_tables/metadata_with_source_block.csv \\
        --categorical subtype,dataset \\
        --numeric z_log1p_n_cells,z_log1p_mean_umi,jeffreys_per_10pct \\
        --exposure jeffreys_per_10pct --out 02_tables/eligibility.csv
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def _require_new(path: Path) -> None:
    if path.exists():
        raise SystemExit(f"Refusing to overwrite: {path}")


def _split(raw: str) -> list[str]:
    return [part.strip() for part in raw.split(",") if part.strip()]


def drop_single_level(frame: pd.DataFrame, columns: list[str]) -> tuple[list[str], list[str]]:
    kept = []
    dropped = []
    for column in columns:
        if column not in frame.columns:
            raise SystemExit(f"missing design column {column}")
        if frame[column].isna().any():
            raise SystemExit(f"missing categorical values: {column}")
        n_level = frame[column].astype(str).nunique(dropna=True)
        if n_level < 2:
            dropped.append(column)
        else:
            kept.append(column)
    return kept, dropped


def design_matrix(
    frame: pd.DataFrame,
    categorical: list[str],
    numeric: list[str],
) -> tuple[np.ndarray, list[str]]:
    pieces = [np.ones((len(frame), 1), dtype=float)]
    names = ["intercept"]
    for column in categorical:
        if frame[column].isna().any():
            raise SystemExit(f"missing categorical values: {column}")
        dummies = pd.get_dummies(frame[column].astype(str), drop_first=True, prefix=column)
        if dummies.shape[1] == 0:
            continue
        pieces.append(dummies.to_numpy(float))
        names.extend(list(dummies.columns))
    for column in numeric:
        if column not in frame.columns:
            raise SystemExit(f"missing numeric column {column}")
        values = frame[column].to_numpy(float)
        if not np.isfinite(values).all():
            raise SystemExit(f"{column} contains non-finite values; no implicit imputation")
        pieces.append(values.reshape(-1, 1))
        names.append(column)
    matrix = np.concatenate(pieces, axis=1)
    return matrix, names


def residual_df(matrix: np.ndarray) -> tuple[int, int, int]:
    rank = int(np.linalg.matrix_rank(matrix))
    n = int(matrix.shape[0])
    return n, rank, n - rank


def flag_arm(
    *,
    n_units: int,
    n_datasets: int,
    n_blocks: int,
    rdf: int,
    rank: int,
    n_cols: int,
    exposure_sd: float,
    n_formal: int,
    n_exploratory: int,
    min_datasets_formal: int,
    min_rdf_formal: int,
    min_rdf_fit: int,
    sd_floor: float,
    part4_flag: str,
) -> tuple[str, str]:
    if part4_flag in {"LOW_N", "NO_EXPOSURE_RANGE"}:
        return "NOT_ESTIMABLE", f"Part 4 forecast {part4_flag} was not overridden"
    if rank < n_cols:
        return "NOT_ESTIMABLE", f"rank-deficient design (rank {rank} < {n_cols} columns)"
    if not np.isfinite(exposure_sd) or exposure_sd < sd_floor:
        return "NOT_ESTIMABLE", "NO_EXPOSURE_RANGE"
    if n_units < n_exploratory or rdf < min_rdf_fit:
        return "NOT_ESTIMABLE", f"n_units={n_units} rdf={rdf} below fit floor"
    if n_units >= n_formal and n_blocks >= min_datasets_formal and rdf >= min_rdf_formal:
        return "formal", ""
    return "exploratory", "below formal n/dataset/rdf gate"


def evaluate(
    meta: pd.DataFrame,
    *,
    arm_key: str,
    categorical: list[str],
    numeric: list[str],
    exposure: str,
    dataset_key: str,
    source_key: str,
    n_formal: int,
    n_exploratory: int,
    min_datasets_formal: int,
    min_rdf_formal: int,
    min_rdf_fit: int,
    sd_floor: float,
    drop_single: bool,
    part4: pd.DataFrame | None,
) -> pd.DataFrame:
    work = meta.loc[meta["eligible"].astype(bool)].copy() if "eligible" in meta.columns else meta.copy()
    if work.empty:
        raise SystemExit("no eligible units")
    rows = []
    groups = [(arm_key, "ALL", work)] if arm_key not in work.columns else [
        (arm_key, str(name), sub) for name, sub in work.groupby(arm_key, observed=True, sort=False)
    ]
    # Always also evaluate the pooled eligible table as ALL if grouping.
    if arm_key in work.columns:
        groups.append((arm_key, "ALL", work))

    for _, arm, sub in groups:
        cat = list(categorical)
        dropped: list[str] = []
        if drop_single:
            cat, dropped = drop_single_level(sub, cat)
        try:
            matrix, names = design_matrix(sub, cat, numeric)
        except SystemExit as exc:
            rows.append(
                {
                    "arm": arm,
                    "n_units": int(len(sub)),
                    "status": "NOT_ESTIMABLE",
                    "reason": str(exc),
                }
            )
            continue
        n, rank, rdf = residual_df(matrix)
        values = sub[exposure].to_numpy(float) if exposure in sub.columns else np.array([])
        sd = float(np.nanstd(values, ddof=1)) if len(values) >= 2 else float("nan")
        part4_flag = ""
        if part4 is not None and not part4.empty:
            first = part4.columns[0]
            mask = part4[first].astype(str).eq(arm)
            if "source" in part4.columns:
                all_mask = mask & part4["source"].astype(str).eq("ALL")
                hit = part4.loc[all_mask if all_mask.any() else mask]
            else:
                hit = part4.loc[mask]
            if not hit.empty and "flag" in hit.columns:
                part4_flag = str(hit["flag"].iloc[0])
        status, reason = flag_arm(
            n_units=int(sub["dataset_donor_id"].nunique()) if "dataset_donor_id" in sub else n,
            n_datasets=int(sub[dataset_key].nunique()) if dataset_key in sub.columns else 0,
            n_blocks=int(sub[source_key].nunique()) if source_key in sub.columns else 0,
            rdf=rdf,
            rank=rank,
            n_cols=int(matrix.shape[1]),
            exposure_sd=sd,
            n_formal=n_formal,
            n_exploratory=n_exploratory,
            min_datasets_formal=min_datasets_formal,
            min_rdf_formal=min_rdf_formal,
            min_rdf_fit=min_rdf_fit,
            sd_floor=sd_floor,
            part4_flag=part4_flag,
        )
        rows.append(
            {
                "arm": arm,
                "n_units": int(sub["dataset_donor_id"].nunique()) if "dataset_donor_id" in sub.columns else n,
                "n_donors": int(sub["dataset_donor_id"].nunique()) if "dataset_donor_id" in sub.columns else n,
                "n_datasets": int(sub[dataset_key].nunique()) if dataset_key in sub.columns else pd.NA,
                "n_source_blocks": int(sub[source_key].nunique()) if source_key in sub.columns else pd.NA,
                "n_cells": int(sub["n_cells"].sum()) if "n_cells" in sub.columns else pd.NA,
                "design_rank": rank,
                "design_columns": int(matrix.shape[1]),
                "rdf": rdf,
                "exposure": exposure,
                "exposure_sd": sd,
                "dropped_single_level": ",".join(dropped),
                "part4_flag": part4_flag,
                "status": status,
                "reason": reason,
                "terms": ",".join(names),
            }
        )
    result = pd.DataFrame(rows)
    result["calibration_status"] = "NOT_ESTABLISHED_FOR_THIS_ANALYSIS"
    result["power_status"] = "NOT_ASSESSED_REQUIRES_DESIGN_SIMULATION"
    result["mde"] = np.nan
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("metadata", type=Path)
    parser.add_argument("--arm-key", default="subtype")
    parser.add_argument("--categorical", default="subtype,dataset")
    parser.add_argument("--numeric", default="z_log1p_n_cells,z_log1p_mean_umi,jeffreys_per_10pct")
    parser.add_argument("--exposure", default="jeffreys_per_10pct")
    parser.add_argument("--dataset-key", default="dataset")
    parser.add_argument("--source-key", default="source_block")
    parser.add_argument("--n-formal", type=int, default=12)
    parser.add_argument("--n-exploratory", type=int, default=8)
    parser.add_argument("--min-datasets-formal", type=int, default=3)
    parser.add_argument("--min-rdf-formal", type=int, default=6)
    parser.add_argument("--min-rdf-fit", type=int, default=4)
    parser.add_argument("--sd-floor", type=float, default=0.5)
    parser.add_argument("--keep-single-level", action="store_true")
    parser.add_argument("--part4", type=Path, default=None)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)

    meta = pd.read_csv(args.metadata)
    part4 = pd.read_csv(args.part4) if args.part4 else None
    table = evaluate(
        meta,
        arm_key=args.arm_key,
        categorical=_split(args.categorical),
        numeric=_split(args.numeric),
        exposure=args.exposure,
        dataset_key=args.dataset_key,
        source_key=args.source_key if args.source_key in meta.columns else args.dataset_key,
        n_formal=args.n_formal,
        n_exploratory=args.n_exploratory,
        min_datasets_formal=args.min_datasets_formal,
        min_rdf_formal=args.min_rdf_formal,
        min_rdf_fit=args.min_rdf_fit,
        sd_floor=args.sd_floor,
        drop_single=not args.keep_single_level,
        part4=part4,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    _require_new(args.out)
    table.to_csv(args.out, index=False, encoding="utf-8-sig")
    audit = {
        "n_arms": int(len(table)),
        "n_formal": int((table["status"] == "formal").sum()),
        "n_exploratory": int((table["status"] == "exploratory").sum()),
        "n_not_estimable": int((table["status"] == "NOT_ESTIMABLE").sum()),
    }
    audit_path = args.out.with_suffix(".audit.json")
    _require_new(audit_path)
    audit_path.write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(f"wrote eligibility ({len(table)} arms): {args.out}")
    print(
        f"formal={audit['n_formal']} exploratory={audit['n_exploratory']} "
        f"NOT_ESTIMABLE={audit['n_not_estimable']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


