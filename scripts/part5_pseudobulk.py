#!/usr/bin/env python3
"""Build target-excluded donor-unit pseudobulk from a locked h5ad.

Counts only. TARGET_GENE is the exposure sidecar, not an outcome gene.
Does not fit limma. Does not reopen labels.

Usage:
    python scripts/part5_pseudobulk.py locked_subtype.h5ad \\
        --gene SYMBOL --group subtype --out 03_pseudobulk/
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

try:
    import anndata as sc
except ImportError:  # pragma: no cover
    sc = None

UNLIKELY_PATTERN = re.compile(
    r"unresolved|stressed|doublet|debris|contaminant|low[\s_-]?qc",
    re.I,
)


def _require_new(path: Path) -> None:
    if path.exists():
        raise SystemExit(f"Refusing to overwrite: {path}")


def counts_csr(adata):
    matrix = adata.layers["counts"]
    if sparse.issparse(matrix):
        return matrix.tocsr()
    return sparse.csr_matrix(np.asarray(matrix))


def gene_counts(matrix, var_names: pd.Index, gene: str) -> np.ndarray:
    if gene not in var_names:
        raise SystemExit(f"TARGET_GENE '{gene}' is not in var_names. Stop.")
    idx = int(var_names.get_loc(gene))
    column = matrix[:, idx]
    if sparse.issparse(column):
        return np.asarray(column.toarray()).ravel()
    return np.asarray(column).ravel()


def unit_key(obs: pd.DataFrame, dataset_key: str, donor_key: str, unit_col: str) -> pd.Series:
    if unit_col in obs.columns:
        return obs[unit_col].astype(str)
    if dataset_key not in obs.columns or donor_key not in obs.columns:
        raise SystemExit(f"Need {unit_col} or both {dataset_key} and {donor_key}.")
    return obs[dataset_key].astype(str) + "_" + obs[donor_key].astype(str)


def jeffreys_per_10pct(n_detected: np.ndarray, n_cells: np.ndarray) -> np.ndarray:
    return 10.0 * (n_detected.astype(float) + 0.5) / (n_cells.astype(float) + 1.0)


def zscore_full(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    sd = np.nanstd(values, ddof=1)
    if not np.isfinite(sd) or sd == 0:
        return np.zeros_like(values)
    return (values - np.nanmean(values)) / sd


def build_pseudobulk(
    adata,
    *,
    gene: str,
    group_key: str,
    dataset_key: str,
    donor_key: str,
    unit_col: str,
    min_cells: int,
    include_group_in_unit: bool,
) -> tuple[sparse.csr_matrix, pd.DataFrame, pd.DataFrame, np.ndarray]:
    if "counts" not in adata.layers:
        raise SystemExit("Need layers['counts'].")
    if not adata.obs_names.is_unique or adata.obs_names.isna().any():
        raise SystemExit("Cell barcodes must be unique and nonmissing before pseudobulk.")
    if min_cells < 1:
        raise SystemExit("min_cells must be positive.")
    if group_key not in adata.obs:
        raise SystemExit(f"missing obs column {group_key}")
    if adata.obs[group_key].isna().any():
        raise SystemExit(f"{group_key} has NA labels. Part 3 is not locked.")

    matrix = counts_csr(adata)
    if matrix.shape[1] != adata.n_vars:
        raise SystemExit("counts layer shape does not match n_vars.")
    if not adata.var_names.is_unique:
        raise SystemExit("Gene names must be unique before pseudobulk.")
    if not np.isfinite(matrix.data).all() or (matrix.data < 0).any() or not np.equal(matrix.data, np.floor(matrix.data)).all():
        raise SystemExit("counts must be finite non-negative integers.")

    obs = adata.obs.copy()
    for key in (dataset_key, donor_key, unit_col):
        if key in obs and (obs[key].isna().any() or obs[key].astype(str).str.strip().eq("").any()):
            raise SystemExit(f"Missing identity in {key}.")
    obs["_donor_unit"] = unit_key(obs, dataset_key, donor_key, unit_col)
    # A string key must encode exactly one dataset/donor pair, in both directions.
    if dataset_key in obs and donor_key in obs:
        identities = obs[[dataset_key, donor_key, "_donor_unit"]].drop_duplicates()
        if identities["_donor_unit"].duplicated().any() or identities.duplicated([dataset_key, donor_key]).any():
            raise SystemExit("Donor-unit keys collide or split a dataset/donor identity.")
    obs["_group"] = obs[group_key].astype(str)
    if include_group_in_unit:
        obs["_row_unit"] = obs["_donor_unit"] + "|" + obs["_group"]
    else:
        obs["_row_unit"] = obs["_donor_unit"]
    components = ["_donor_unit", "_group"] if include_group_in_unit else ["_donor_unit"]
    if obs[["_row_unit", *components]].drop_duplicates()["_row_unit"].duplicated().any():
        raise SystemExit("Row-unit keys collide across donor/group identities.")

    gene_vec = gene_counts(matrix, adata.var_names, gene)
    library = np.asarray(matrix.sum(axis=1)).ravel()
    outcome_mask = np.asarray(adata.var_names != gene)
    outcome_matrix = matrix[:, outcome_mask]
    outcome_library = library - gene_vec
    detected = (gene_vec > 0).astype(int)

    rows = []
    unit_sums = []
    for unit_id, sub in obs.groupby("_row_unit", sort=False):
        iloc = adata.obs.index.get_indexer(sub.index)
        if (iloc < 0).any():
            raise SystemExit("unit barcodes do not match adata.obs.index")
        n = int(len(iloc))
        n_det = int(detected[iloc].sum())
        gene_umi = float(gene_vec[iloc].sum())
        total_umi = float(library[iloc].sum())
        mean_umi = float(np.mean(outcome_library[iloc])) if n else float("nan")
        dataset = str(sub[dataset_key].iloc[0]) if dataset_key in sub else ""
        donor = str(sub[donor_key].iloc[0]) if donor_key in sub else ""
        group = str(sub["_group"].iloc[0]) if sub["_group"].nunique() == 1 else "POOLED"
        rows.append(
            {
                "unit_id": str(unit_id),
                "dataset_donor_id": str(sub["_donor_unit"].iloc[0]),
                dataset_key: dataset,
                donor_key: donor,
                group_key: group,
                "n_cells": n,
                "n_detected": n_det,
                "detection_fraction": n_det / n if n else float("nan"),
                "jeffreys_per_10pct": float(jeffreys_per_10pct(np.array([n_det]), np.array([n]))[0]),
                "gene_umi": gene_umi,
                "total_umi": total_umi,
                "mean_umi_per_cell": mean_umi,
                "log2cpm": float(np.log2(1e6 * gene_umi / total_umi + 1.0)) if total_umi > 0 else float("nan"),
                "eligible": n >= min_cells,
                "unlikely_part5_arm": any(bool(UNLIKELY_PATTERN.search(label)) for label in sub["_group"].unique()),
            }
        )
        unit_sums.append(outcome_matrix[iloc].sum(axis=0))

    meta = pd.DataFrame(rows)
    if meta.empty:
        raise SystemExit("No units constructed.")
    stacked = sparse.vstack([sparse.csr_matrix(block) for block in unit_sums])
    counts = stacked.T.tocsr()

    eligible = meta["eligible"].to_numpy()
    meta["z_log1p_n_cells"] = np.nan
    meta["z_log1p_mean_umi"] = np.nan
    meta["z_log2cpm"] = np.nan
    if eligible.any():
        meta.loc[eligible, "z_log1p_n_cells"] = zscore_full(np.log1p(meta.loc[eligible, "n_cells"].to_numpy(float)))
        meta.loc[eligible, "z_log1p_mean_umi"] = zscore_full(np.log1p(meta.loc[eligible, "mean_umi_per_cell"].to_numpy(float)))
        meta.loc[eligible, "z_log2cpm"] = zscore_full(meta.loc[eligible, "log2cpm"].to_numpy(float))

    genes = pd.DataFrame(
        {
            "gene": list(adata.var_names[outcome_mask].astype(str)),
            "is_target": False,
        }
    )
    return counts, meta, genes, gene_vec


def write_pseudobulk(
    counts: sparse.spmatrix,
    meta: pd.DataFrame,
    genes: pd.DataFrame,
    out_dir: Path,
    *,
    gene: str,
) -> None:
    if counts.shape != (len(genes), len(meta)) or genes["gene"].eq(gene).any():
        raise SystemExit("Outcome matrix must align with genes/units and exclude TARGET_GENE.")
    out_dir.mkdir(parents=True, exist_ok=True)
    mtx_path = out_dir / "counts.mtx"
    meta_path = out_dir / "metadata.csv"
    genes_path = out_dir / "genes.csv"
    audit_path = out_dir / "pseudobulk_audit.json"
    for path in (mtx_path, meta_path, genes_path, audit_path):
        _require_new(path)

    mmwrite(mtx_path, counts)
    meta.to_csv(meta_path, index=False, encoding="utf-8-sig")
    genes.to_csv(genes_path, index=False, encoding="utf-8-sig")
    audit = {
        "target_gene": gene,
        "n_genes": int(counts.shape[0]),
        "n_units": int(counts.shape[1]),
        "n_eligible_units": int(meta["eligible"].sum()),
        "target_in_outcome_matrix": False,
        "outcome_gene_manifest": "genes.csv",
        "layer": "counts",
    }
    audit_path.write_text(json.dumps(audit, indent=2), encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("h5ad")
    parser.add_argument("--gene", required=True)
    parser.add_argument("--group", default="subtype")
    parser.add_argument("--dataset-key", default="dataset")
    parser.add_argument("--donor-key", default="donor_id")
    parser.add_argument("--unit-key", default="dataset_donor_id")
    parser.add_argument("--min-cells", type=int, default=20)
    parser.add_argument("--pool-groups", action="store_true", help="one row per donor, pooling subtypes")
    parser.add_argument("--out", required=True)
    args = parser.parse_args(argv)
    if sc is None:
        raise SystemExit("anndata required")

    adata = sc.read_h5ad(args.h5ad)
    counts, meta, genes, _ = build_pseudobulk(
        adata,
        gene=args.gene,
        group_key=args.group,
        dataset_key=args.dataset_key,
        donor_key=args.donor_key,
        unit_col=args.unit_key,
        min_cells=args.min_cells,
        include_group_in_unit=not args.pool_groups,
    )
    write_pseudobulk(counts, meta, genes, Path(args.out), gene=args.gene)
    print(f"wrote {counts.shape[0]} genes × {counts.shape[1]} units under {args.out}")
    print(f"eligible units: {int(meta['eligible'].sum())} / {len(meta)}")
    print("TARGET_GENE is the exposure sidecar, not an outcome gene.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

