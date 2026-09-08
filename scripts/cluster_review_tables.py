#!/usr/bin/env python3
"""Export Leiden top-20 markers and per-cluster QC for human review.

Part 1: Wilcoxon worksheet (proposed labels stay blank).
Part 3: add --strict-positive and --exclude-genes TARGET_GENE; still no names.

Does not name clusters and does not delete cells.

Usage:
    python scripts/cluster_review_tables.py clustered.h5ad \\
        --leiden-key leiden_r0_5 \\
        --out tables/Leiden_markers \\
        --lineage-genes ACTA2,PECAM1,CD3D,CD79A,LYZ
    python scripts/cluster_review_tables.py round_leiden.h5ad \\
        --leiden-key leiden_r0_5 --out tables/round_00_initial/markers \\
        --strict-positive --exclude-genes TARGET_GENE \\
        --lineage-genes PTPRC,CD3D,PECAM1,ACTA2,CD79A,LYZ
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse as sp

STRESS_SYMBOLS = {"MALAT1", "NEAT1", "KCNQ1OT1", "XIST"}


def _median(series: pd.Series) -> float:
    return float(series.median()) if len(series) else float("nan")


def markers(adata, leiden_key: str, n_genes: int) -> pd.DataFrame:
    import scanpy as sc
    rank_key = f"rank_genes_{leiden_key}"
    sc.tl.rank_genes_groups(
        adata,
        groupby=leiden_key,
        method="wilcoxon",
        use_raw=False,
        layer="normalized",
        n_genes=n_genes,
        pts=True,
        tie_correct=True,
        key_added=rank_key,
    )
    df = sc.get.rank_genes_groups_df(adata, group=None, key=rank_key)
    df = df.rename(
        columns={
            "group": "cluster",
            "names": "gene",
            "scores": "score",
            "pvals": "pval",
            "pvals_adj": "pval_adj",
            "logfoldchanges": "logfoldchange",
        }
    )
    df["cluster"] = df["cluster"].astype(str)
    df["rank"] = df.groupby("cluster", observed=True).cumcount() + 1
    keep = [
        c
        for c in [
            "cluster",
            "gene",
            "score",
            "pval",
            "pval_adj",
            "logfoldchange",
            "rank",
            "pct_nz_group",
            "pct_nz_reference",
        ]
        if c in df.columns
    ]
    return df[keep]


def strict_positive(df: pd.DataFrame) -> pd.DataFrame:
    """Cluster-vs-rest display filter used in Part 3 marker review."""
    out = df.copy()
    if "score" in out:
        out = out.loc[out["score"] > 0]
    if "logfoldchange" in out:
        out = out.loc[out["logfoldchange"] > 0]
    if "pval_adj" in out:
        out = out.loc[out["pval_adj"] <= 0.05]
    if {"pct_nz_group", "pct_nz_reference"}.issubset(out.columns):
        out = out.loc[out["pct_nz_group"] > out["pct_nz_reference"]]
        out = out.loc[out["pct_nz_group"] >= 0.10]
    out = out.copy()
    out["rank"] = out.groupby("cluster", observed=True).cumcount() + 1
    return out


def cluster_qc(
    adata,
    leiden_key: str,
    top20: pd.DataFrame,
    lineage: set[str],
    target_gene: str | None,
) -> pd.DataFrame:
    obs = adata.obs.copy()
    obs[leiden_key] = obs[leiden_key].astype(str)
    rows = []
    for cluster, sub in obs.groupby(leiden_key, observed=True):
        genes20 = (
            top20.loc[top20["cluster"].eq(cluster) & top20["rank"].le(20), "gene"]
            .astype(str)
            .str.upper()
            .tolist()
        )
        gene_set = set(genes20)
        n_ds = sub["dataset"].nunique() if "dataset" in sub else 1
        if "dataset" in sub:
            max_frac = float(sub["dataset"].value_counts(normalize=True).iloc[0])
        else:
            max_frac = 1.0
        det = float("nan")
        if target_gene and target_gene in adata.var_names and "counts" in adata.layers:
            x = adata[sub.index, target_gene].layers["counts"]
            if sp.issparse(x):
                vals = np.asarray(x.toarray()).ravel()
            else:
                vals = np.asarray(x).ravel()
            det = float((vals > 0).mean())
        rows.append(
            {
                "cluster": cluster,
                "n_cells": int(len(sub)),
                "n_donors": (int(sub["dataset_donor_id"].nunique()) if "dataset_donor_id" in sub
                             else int(sub[["dataset", "donor_id"]].drop_duplicates().shape[0])
                             if {"dataset", "donor_id"}.issubset(sub.columns)
                             else int(sub["donor_id"].nunique()) if "donor_id" in sub else pd.NA),
                "n_datasets": int(n_ds),
                "n_libraries": int(sub["library_id"].nunique())
                if "library_id" in sub
                else pd.NA,
                "median_total_counts": _median(sub["total_counts"])
                if "total_counts" in sub
                else pd.NA,
                "median_n_genes": _median(sub["n_genes_by_counts"])
                if "n_genes_by_counts" in sub
                else pd.NA,
                "median_pct_mt": _median(sub["pct_counts_mt"])
                if "pct_counts_mt" in sub
                else pd.NA,
                "median_pct_hb": _median(sub["pct_counts_hb"])
                if "pct_counts_hb" in sub
                else pd.NA,
                "median_pct_ribo": _median(sub["pct_counts_ribo"])
                if "pct_counts_ribo" in sub
                else pd.NA,
                "median_pct_top50": _median(sub["pct_counts_in_top_50_genes"])
                if "pct_counts_in_top_50_genes" in sub
                else pd.NA,
                "median_log10_genes_per_umi": _median(sub["log10_genes_per_umi"])
                if "log10_genes_per_umi" in sub
                else pd.NA,
                "median_scrublet_score": _median(sub["scrublet_doublet_score"])
                if "scrublet_doublet_score" in sub
                else pd.NA,
                "max_dataset_fraction": max_frac,
                "top20_contains_malat1_neat1": bool(gene_set & STRESS_SYMBOLS),
                "lineage_panel_hits": int(len(gene_set & {g.upper() for g in lineage})),
                "target_gene_detection_fraction": det,
                "proposed_label": "",
            }
        )
    return pd.DataFrame(rows).sort_values(
        "cluster", key=lambda s: s.map(lambda x: int(x) if str(x).isdigit() else x)
    )


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("h5ad")
    p.add_argument("--leiden-key", required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--n-genes", type=int, default=50)
    p.add_argument("--lineage-genes", default="", help="comma-separated symbols")
    p.add_argument("--target-gene", default="")
    p.add_argument(
        "--strict-positive",
        action="store_true",
        help="Part 3 filter: score>0, logFC>0, adjP<=0.05, pct_in>pct_out, pct_in>=0.10",
    )
    p.add_argument(
        "--exclude-genes",
        default="",
        help="comma-separated symbols dropped from annotation worksheets (e.g. TARGET_GENE)",
    )
    args = p.parse_args()

    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(args.h5ad)
    if args.leiden_key not in adata.obs:
        raise SystemExit(f"missing obs column {args.leiden_key}")
    if "normalized" not in adata.layers:
        raise SystemExit("missing layers['normalized']")

    lineage = {g.strip() for g in args.lineage_genes.split(",") if g.strip()}
    target = args.target_gene.strip() or None

    mark = markers(adata, args.leiden_key, args.n_genes)
    mark.insert(0, "leiden_key", args.leiden_key)
    mark.to_csv(out / f"{args.leiden_key}_top50.csv", index=False, encoding="utf-8-sig")
    display = mark.copy()
    if args.strict_positive:
        display = strict_positive(display)
        display.to_csv(
            out / f"{args.leiden_key}_top50_strict_positive.csv",
            index=False,
            encoding="utf-8-sig",
        )
    exclude = {g.strip().upper() for g in args.exclude_genes.split(",") if g.strip()}
    if exclude:
        display = display.loc[~display["gene"].astype(str).str.upper().isin(exclude)]
        display["rank"] = display.groupby("cluster", observed=True).cumcount() + 1
    top20 = display.loc[display["rank"] <= 20].copy()
    top20.to_csv(out / f"{args.leiden_key}_top20.csv", index=False, encoding="utf-8-sig")

    sizes = (
        adata.obs[args.leiden_key]
        .astype(str)
        .value_counts()
        .rename_axis("cluster")
        .rename("n_cells")
        .reset_index()
    )
    compact = (
        top20.sort_values(["cluster", "rank"])
        .groupby("cluster", observed=True)["gene"]
        .apply(lambda s: ", ".join(s.astype(str)))
        .rename("top20_markers")
        .reset_index()
    )
    compact = sizes.merge(compact, on="cluster", how="left")
    compact["top20_markers"] = compact["top20_markers"].fillna("")
    compact["proposed_label"] = ""
    compact[["cluster", "n_cells", "top20_markers", "proposed_label"]].to_csv(
        out / f"{args.leiden_key}_top20_by_cluster.csv",
        index=False,
        encoding="utf-8-sig",
    )

    qc = cluster_qc(adata, args.leiden_key, top20, lineage, target)
    qc.to_csv(out / f"{args.leiden_key}_cluster_qc.csv", index=False, encoding="utf-8-sig")

    lines = [
        f"# Marker annotation worksheet: {args.leiden_key}",
        "",
        "- Marker test: wilcoxon on layers['normalized']",
        "- Proposed labels are blank for manual review.",
    ]
    if args.strict_positive:
        lines.append(
            "- Part 3: do not name clusters here. Fill the KEEP/DELETE decision CSV."
        )
    lines.extend(
        [
            "",
            "| Cluster | Cells | Top 20 markers | Proposed label |",
            "|---:|---:|---|---|",
        ]
    )
    for row in compact.sort_values(
        "cluster", key=lambda s: s.map(lambda x: int(x) if str(x).isdigit() else x)
    ).itertuples(index=False):
        lines.append(
            f"| {row.cluster} | {int(row.n_cells)} | {row.top20_markers} |  |"
        )
    (out / f"{args.leiden_key}_marker_annotation.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

