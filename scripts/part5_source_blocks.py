#!/usr/bin/env python3
"""Map datasets to source_block, audit exposure range, relabel donor LOO.

Does not merge GSE accessions automatically. Consecutive IDs with
interleaved donor numbering are a suggestion for the protocol.

Usage:
    python scripts/part5_source_blocks.py 03_pseudobulk/metadata.csv \\
        --map 00_protocol_manifest/source_block_map.yaml \\
        --exposure jeffreys_per_10pct --out 02_tables/
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import yaml
except ImportError:  # pragma: no cover
    yaml = None


GSE_RE = re.compile(r"GSE(\d+)", re.I)
DONOR_NUM_RE = re.compile(r"(\d+)$")


def _require_new(path: Path) -> None:
    if path.exists():
        raise SystemExit(f"Refusing to overwrite: {path}")


def load_map(path: Path | None) -> dict[str, str]:
    if path is None:
        return {}
    if yaml is None:
        raise SystemExit("PyYAML required: pip install pyyaml")
    if not path.is_file():
        raise SystemExit(f"missing source_block map: {path}")
    raw = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if not isinstance(raw, dict):
        raise SystemExit("source_block map must be dataset → block")
    if any(not isinstance(v, str) or not v.strip() for v in raw.values()):
        raise SystemExit("Each source block must be a nonempty string.")
    return {str(k): str(v) for k, v in raw.items() if not str(k).startswith("#")}


def gse_number(name: str) -> int | None:
    match = GSE_RE.search(str(name))
    if not match:
        return None
    return int(match.group(1))


def donor_tail(name: str) -> int | None:
    match = DONOR_NUM_RE.search(str(name).strip())
    if not match:
        return None
    return int(match.group(1))


def suggest_siblings(meta: pd.DataFrame, dataset_key: str, donor_key: str) -> pd.DataFrame:
    datasets = list(dict.fromkeys(meta[dataset_key].astype(str)))
    numbered = [(name, gse_number(name)) for name in datasets]
    numbered = [(name, num) for name, num in numbered if num is not None]
    rows = []
    for i, (a, na) in enumerate(numbered):
        for b, nb in numbered[i + 1 :]:
            if abs(na - nb) != 1:
                continue
            donors_a = {donor_tail(x) for x in meta.loc[meta[dataset_key].astype(str).eq(a), donor_key].astype(str)}
            donors_b = {donor_tail(x) for x in meta.loc[meta[dataset_key].astype(str).eq(b), donor_key].astype(str)}
            donors_a.discard(None)
            donors_b.discard(None)
            if not donors_a or not donors_b:
                interleaved = False
            else:
                interleaved = bool(not donors_a & donors_b and
                                   max(min(donors_a), min(donors_b)) < min(max(donors_a), max(donors_b)))
            rows.append(
                {
                    "dataset_a": a,
                    "dataset_b": b,
                    "consecutive_gse": True,
                    "interleaved_donor_numbers": interleaved,
                    "suggestion": "SAME_STUDY_LIKELY" if interleaved else "UNRESOLVED",
                    "note": "Do not merge unless the protocol freezes this map.",
                }
            )
    return pd.DataFrame(rows)


def between_within_fraction(values: np.ndarray, groups: np.ndarray) -> dict[str, float]:
    frame = pd.DataFrame({"y": values, "g": groups}).dropna()
    if len(frame) < 2:
        return {"between_fraction": float("nan"), "within_fraction": float("nan")}
    y = frame["y"].to_numpy(float)
    grand = float(np.mean(y))
    ss_total = float(np.sum((y - grand) ** 2))
    if ss_total == 0:
        return {"between_fraction": float("nan"), "within_fraction": float("nan")}
    ss_between = 0.0
    for _, sub in frame.groupby("g"):
        n = len(sub)
        ss_between += n * (float(sub["y"].mean()) - grand) ** 2
    between = ss_between / ss_total
    return {"between_fraction": between, "within_fraction": max(0.0, 1.0 - between)}


def exposure_table(
    meta: pd.DataFrame,
    *,
    group_key: str,
    exposure: str,
    min_sd: float,
) -> pd.DataFrame:
    rows = []
    for (dataset, block), sub in meta.groupby([group_key, "source_block"], observed=True, sort=False):
        values = sub[exposure].to_numpy(float)
        finite = values[np.isfinite(values)]
        sd = float(np.std(finite, ddof=1)) if finite.size >= 2 else float("nan")
        q75, q25 = (np.percentile(finite, [75, 25]) if finite.size else (np.nan, np.nan))
        rows.append(
            {
                group_key: str(dataset),
                "source_block": str(block),
                "n_units": int(len(sub)),
                "n_donors": int(sub["dataset_donor_id"].nunique()) if "dataset_donor_id" in sub else int(len(sub)),
                "n_cells": int(sub["n_cells"].sum()) if "n_cells" in sub else pd.NA,
                "mean": float(np.mean(finite)) if finite.size else float("nan"),
                "sd": sd,
                "min": float(np.min(finite)) if finite.size else float("nan"),
                "max": float(np.max(finite)) if finite.size else float("nan"),
                "iqr": float(q75 - q25) if finite.size else float("nan"),
                "flag": "NO_EXPOSURE_RANGE" if (not np.isfinite(sd) or sd < min_sd) else "RANGE_OK",
                "exposure": exposure,
            }
        )
    return pd.DataFrame(rows)


def relabel_loo(loo: pd.DataFrame, meta: pd.DataFrame, dataset_key: str) -> pd.DataFrame:
    identities = meta[["dataset_donor_id", "source_block", dataset_key]].drop_duplicates()
    if identities["dataset_donor_id"].duplicated().any():
        raise SystemExit("A donor maps to conflicting dataset/source identities.")
    lookup = (
        meta.drop_duplicates("dataset_donor_id")
        .set_index("dataset_donor_id")[["source_block", dataset_key]]
        .to_dict("index")
    )
    out = loo.copy()
    donor_col = "dataset_donor_id" if "dataset_donor_id" in out.columns else "omitted_donor"
    if donor_col not in out.columns:
        raise SystemExit("LOO table needs dataset_donor_id or omitted_donor")
    if out[donor_col].isna().any() or not out[donor_col].astype(str).isin(lookup).all():
        raise SystemExit("LOO donor not found in metadata.")
    out["source_block"] = out[donor_col].astype(str).map(lambda x: lookup.get(x, {}).get("source_block", ""))
    out[dataset_key] = out[donor_col].astype(str).map(lambda x: lookup.get(x, {}).get(dataset_key, ""))
    return out


def loo_unbalanced(loo: pd.DataFrame, threshold: float = 0.75) -> dict[str, float | str | bool]:
    if "source_block" not in loo.columns or loo.empty:
        return {"loo_source_unbalanced": False, "dominant_block": "", "fraction": float("nan")}
    counts = loo["source_block"].astype(str).value_counts(normalize=True)
    dominant = str(counts.index[0])
    fraction = float(counts.iloc[0])
    return {
        "loo_source_unbalanced": fraction >= threshold,
        "dominant_block": dominant,
        "fraction": fraction,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("metadata", type=Path)
    parser.add_argument("--map", type=Path, default=None)
    parser.add_argument("--dataset-key", default="dataset")
    parser.add_argument("--donor-key", default="donor_id")
    parser.add_argument("--exposure", default="jeffreys_per_10pct")
    parser.add_argument("--min-sd", type=float, default=0.5)
    parser.add_argument("--all-units", action="store_true", help="include ineligible units")
    parser.add_argument("--loo", type=Path, default=None)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)

    meta = pd.read_csv(args.metadata)
    if args.dataset_key not in meta.columns:
        raise SystemExit(f"metadata missing {args.dataset_key}")
    if args.exposure not in meta.columns:
        raise SystemExit(f"metadata missing {args.exposure}")
    for key in (args.dataset_key, args.donor_key, "dataset_donor_id"):
        if key in meta and (meta[key].isna().any() or meta[key].astype(str).str.strip().eq("").any()):
            raise SystemExit(f"Missing identity in {key}.")
    if "dataset_donor_id" not in meta.columns:
        if args.donor_key in meta.columns:
            meta["dataset_donor_id"] = meta[args.dataset_key].astype(str) + "_" + meta[args.donor_key].astype(str)
        else:
            raise SystemExit("metadata needs dataset_donor_id or donor_id")
    if "eligible" in meta.columns and not args.all_units:
        eligible = meta["eligible"].astype(str).str.lower().map({"true": True, "false": False, "1": True, "0": False})
        if eligible.isna().any():
            raise SystemExit("eligible must be true/false or 1/0.")
        work = meta.loc[eligible].copy()
    else:
        work = meta.copy()
    for key in (args.dataset_key, "dataset_donor_id"):
        if work[key].isna().any() or work[key].astype(str).str.strip().eq("").any():
            raise SystemExit(f"Missing identity in {key}.")
    if not np.isfinite(work[args.exposure].to_numpy(float)).all():
        raise SystemExit("Eligible exposure must be finite.")

    mapping = load_map(args.map)
    work["source_block"] = work[args.dataset_key].astype(str).map(lambda x: mapping.get(x, x))
    identities = work[["dataset_donor_id", args.dataset_key, "source_block"]].drop_duplicates()
    if identities["dataset_donor_id"].duplicated().any():
        raise SystemExit("A donor maps to conflicting dataset/source identities.")
    if mapping:
        unknown = sorted(set(work[args.dataset_key].astype(str)) - set(mapping))
        if unknown:
            raise ValueError(f"Datasets missing from source-block map: {unknown}. Freeze the map before LODO.")

    args.out.mkdir(parents=True, exist_ok=True)
    donor_col = "dataset_donor_id" if "dataset_donor_id" in work.columns else args.dataset_key
    grouped = work.groupby([args.dataset_key, "source_block"], observed=True)
    block_map = grouped.agg(n_units=(donor_col, "size"), n_donors=(donor_col, "nunique")).reset_index()
    if "n_cells" in work.columns:
        block_map["n_cells"] = grouped["n_cells"].sum().to_numpy()
    exposure = exposure_table(work, group_key=args.dataset_key, exposure=args.exposure, min_sd=args.min_sd)
    suggestions = suggest_siblings(work, args.dataset_key, args.donor_key)
    var_share = between_within_fraction(work[args.exposure].to_numpy(float), work[args.dataset_key].astype(str).to_numpy())

    for name, table in {
        "source_block_map.csv": block_map,
        "exposure_range_by_dataset.csv": exposure,
        "source_block_suggestions.csv": suggestions,
    }.items():
        path = args.out / name
        _require_new(path)
        table.to_csv(path, index=False, encoding="utf-8-sig")

    annotated_path = args.out / "metadata_with_source_block.csv"
    _require_new(annotated_path)
    work.to_csv(annotated_path, index=False, encoding="utf-8-sig")

    loo_audit = {}
    if args.loo is not None:
        loo = relabel_loo(pd.read_csv(args.loo), work, args.dataset_key)
        loo_path = args.out / "loo_with_source.csv"
        _require_new(loo_path)
        loo.to_csv(loo_path, index=False, encoding="utf-8-sig")
        loo_audit = loo_unbalanced(loo)

    audit = {
        "exposure": args.exposure,
        "n_datasets": int(work[args.dataset_key].nunique()),
        "n_source_blocks": int(work["source_block"].nunique()),
        "mapped": bool(mapping),
        "variance": var_share,
        "loo": loo_audit,
        "note": "Suggestions are not a merge. Freeze the map in the protocol.",
    }
    audit_path = args.out.parent.joinpath("05_logs") / "exposure_variance.json"
    if "02_tables" in args.out.as_posix() or args.out.name == "02_tables":
        logs = args.out.parent / "05_logs"
        logs.mkdir(parents=True, exist_ok=True)
        audit_path = logs / "exposure_variance.json"
    else:
        audit_path = args.out / "exposure_variance.json"
    _require_new(audit_path)
    audit_path.write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(f"source blocks: {audit['n_source_blocks']} from {audit['n_datasets']} datasets")
    if not suggestions.empty:
        print(f"sibling suggestions: {len(suggestions)} pair(s) — protocol must freeze the merge")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

