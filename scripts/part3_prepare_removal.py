#!/usr/bin/env python3
"""Archive DELETE clusters and write a child raw-count checkpoint.

Does not recompute HVG, PCA, Harmony, neighbors, UMAP, or Leiden.
Does not name clusters. Overwrite policy: forbid.

Usage:
    python scripts/part3_prepare_removal.py \\
        --parent-raw objects/round_00_initial/raw_counts.h5ad \\
        --parent-clustered objects/round_00_initial/leiden.h5ad \\
        --leiden-key leiden_r0_5 \\
        --decision tables/round_00_initial/manual_decision.csv \\
        --removed-out objects/round_00_initial/removed_clusters.h5ad \\
        --child-raw objects/round_01_after_removal/raw_counts.h5ad \\
        --tables-out tables/round_00_initial
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import anndata as ad
except ImportError:  # pragma: no cover
    ad = None

KEEP = {"keep", "retain"}
DELETE = {"delete", "remove"}
EMBEDDING_UNS = {"neighbors", "pca", "umap", "harmony", "paga"}
DROP_OBS_PREFIXES = ("leiden",)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(2**20):
            digest.update(block)
    return digest.hexdigest()


def _require_new(path: Path) -> Path:
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite existing artifact: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def normalize_decision(value: object) -> str:
    text = str(value).strip().upper()
    key = text.lower()
    if key in KEEP:
        return "KEEP"
    if key in DELETE:
        return "DELETE"
    if text in {"KEEP", "DELETE"}:
        return text
    raise ValueError(
        f"Illegal decision {value!r}. Use KEEP or DELETE "
        "(aliases: keep/retain, delete/remove). Names do not belong here."
    )


def load_decisions(path: Path, leiden_key: str, cluster_ids: set[str]) -> pd.DataFrame:
    required = {
        "selected_leiden_column",
        "cluster_id",
        "decision",
        "reason",
        "reviewer",
        "review_date",
    }
    frame = pd.read_csv(path, dtype=str).fillna("")
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"Decision file missing columns: {sorted(missing)}")
    frame["cluster_id"] = frame["cluster_id"].astype(str).str.strip()
    frame["decision"] = frame["decision"].map(normalize_decision)
    if frame["cluster_id"].eq("").any():
        raise ValueError("Empty cluster_id in decision file.")
    if frame["cluster_id"].duplicated().any():
        raise ValueError("Duplicate cluster_id in decision file.")
    if not frame["selected_leiden_column"].eq(leiden_key).all():
        raise ValueError("Decision file selected_leiden_column does not match --leiden-key.")
    listed = set(frame["cluster_id"])
    extra = listed - cluster_ids
    missing_ids = cluster_ids - listed
    if extra or missing_ids:
        raise ValueError(
            f"Decision file must cover every cluster exactly once. "
            f"extra={sorted(extra)} missing={sorted(missing_ids, key=lambda x: int(x) if x.isdigit() else x)}"
        )
    deletes = frame.loc[frame["decision"].eq("DELETE")]
    if deletes["reason"].str.strip().eq("").any():
        raise ValueError("Every DELETE row needs a non-empty human reason.")
    for field in ("reviewer", "review_date"):
        if frame[field].str.strip().eq("").any():
            raise ValueError(f"Every decision row needs {field}")
    if deletes.empty:
        raise ValueError(
            "Decision file has no DELETE rows. Do not run prepare-removal; "
            "go to STOP POINT 3 naming instead."
        )
    return frame


def clean_analysis_state(adata) -> None:
    """Drop embeddings so the child round cannot reuse the parent manifold."""
    for key in list(adata.obsm.keys()):
        del adata.obsm[key]
    for key in list(adata.obsp.keys()):
        del adata.obsp[key]
    for key in list(adata.varm.keys()):
        del adata.varm[key]
    adata.raw = None
    stale = [key for key in adata.var if key.startswith("highly_variable") or key in
             {"means", "dispersions", "dispersions_norm", "variances", "variances_norm"}]
    adata.var = adata.var.drop(columns=stale)
    adata.uns.pop("hvg", None)
    for key in list(adata.uns.keys()):
        lowered = str(key).lower()
        if key in EMBEDDING_UNS or lowered.startswith("leiden") or lowered.startswith("rank_genes"):
            del adata.uns[key]
    drop_obs = [
        col
        for col in adata.obs.columns
        if col.lower().startswith(DROP_OBS_PREFIXES) or col.lower().startswith("rank_genes")
    ]
    if drop_obs:
        adata.obs = adata.obs.drop(columns=drop_obs)
    if "normalized" in adata.layers:
        del adata.layers["normalized"]
    if "counts" not in adata.layers:
        raise SystemExit("Parent raw checkpoint is missing layers['counts'].")
    adata.X = adata.layers["counts"].copy()


def partition(raw, clustered, leiden_key: str, delete_ids: set[str]):
    if leiden_key not in clustered.obs:
        raise SystemExit(f"Clustered object missing {leiden_key}")
    if not raw.obs_names.is_unique or not clustered.obs_names.is_unique:
        raise SystemExit("Parent barcodes must be unique")
    if set(raw.obs_names) != set(clustered.obs_names):
        raise SystemExit("Raw and clustered parent must contain exactly the same barcodes")
    if clustered.obs[leiden_key].isna().any():
        raise SystemExit("Missing parent cluster labels")
    labels = clustered.obs.loc[raw.obs_names, leiden_key].astype(str)
    if not delete_ids or not delete_ids.issubset(set(labels)):
        raise SystemExit("DELETE IDs must be nonempty and present in the parent")
    mask = labels.isin(delete_ids).to_numpy()
    if mask.all():
        raise SystemExit("Removal would leave no cells; stop the lineage")
    removed = raw[mask].copy()
    retained = raw[~mask].copy()
    if removed.n_obs + retained.n_obs != raw.n_obs:
        raise AssertionError("Removal partition integrity failed: counts do not sum to parent.")
    if set(removed.obs_names) & set(retained.obs_names):
        raise AssertionError("Removal partition integrity failed: barcode overlap.")
    return removed, retained, labels, mask


def _audit_columns(obs: pd.DataFrame) -> list[str]:
    preferred = [
        "dataset",
        "donor_id",
        "library_id",
        "sample_id",
        "total_counts",
        "n_genes_by_counts",
        "pct_counts_mt",
        "scrublet_doublet_score",
        "previous_round_cluster",
        "exclusion_reason",
    ]
    return [c for c in preferred if c in obs.columns]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parent-raw", type=Path, required=True)
    parser.add_argument("--parent-clustered", type=Path, required=True)
    parser.add_argument("--leiden-key", required=True)
    parser.add_argument("--decision", type=Path, required=True)
    parser.add_argument("--removed-out", type=Path, required=True)
    parser.add_argument("--child-raw", type=Path, required=True)
    parser.add_argument("--tables-out", type=Path, required=True)
    parser.add_argument("--child-tables", type=Path, default=None)
    parser.add_argument("--child-round-name", default="")
    args = parser.parse_args(argv)
    if ad is None:
        raise SystemExit("anndata required")

    clustered = ad.read_h5ad(args.parent_clustered)
    cluster_ids = set(clustered.obs[args.leiden_key].astype(str))
    decisions = load_decisions(args.decision, args.leiden_key, cluster_ids)
    delete_ids = set(decisions.loc[decisions["decision"].eq("DELETE"), "cluster_id"])
    raw = ad.read_h5ad(args.parent_raw)
    if "counts" not in raw.layers:
        raise SystemExit("Parent raw checkpoint is missing layers['counts'].")

    removed, retained, labels, mask = partition(raw, clustered, args.leiden_key, delete_ids)
    reason_map = dict(
        zip(
            decisions.loc[decisions["decision"].eq("DELETE"), "cluster_id"],
            decisions.loc[decisions["decision"].eq("DELETE"), "reason"],
        )
    )
    removed.obs["previous_round_cluster"] = labels.loc[removed.obs_names].to_numpy()
    removed.obs["exclusion_reason"] = removed.obs["previous_round_cluster"].map(reason_map)
    retained.obs["previous_round_cluster"] = labels.loc[retained.obs_names].to_numpy()
    history = "; ".join(
        f"removed {row.cluster_id}: {row.reason}"
        for row in decisions.loc[decisions["decision"].eq("DELETE")].itertuples()
    )
    retained.obs["removal_history"] = history
    child_name = args.child_round_name.strip() or _infer_child_name(args.child_raw)
    retained.obs["analysis_round"] = child_name

    clean_analysis_state(retained)

    removed_path = _require_new(args.removed_out)
    child_path = _require_new(args.child_raw)
    tables = args.tables_out
    tables.mkdir(parents=True, exist_ok=True)

    removed.write_h5ad(removed_path, compression="gzip")
    retained.write_h5ad(child_path, compression="gzip")

    audit = removed.obs.loc[:, _audit_columns(removed.obs)].copy()
    audit.insert(0, "cell_id", removed.obs_names.astype(str))
    _write_new(audit, tables / "removed_clusters_audit.csv")
    for field, name in (
        ("dataset", "removed_clusters_by_dataset.csv"),
        ("donor_id", "removed_clusters_by_donor.csv"),
        ("library_id", "removed_clusters_by_library.csv"),
        ("sample_id", "removed_clusters_by_sample.csv"),
    ):
        if field in removed.obs:
            summary = (
                removed.obs.groupby(field, observed=True).size().rename("n_cells").reset_index()
            )
            _write_new(summary, tables / name)

    manifest = pd.DataFrame(
        {
            "parent_cell_id": raw.obs_names.astype(str),
            "status": np.where(mask, "removed", "retained"),
            "parent_cluster": labels.to_numpy(),
        }
    )
    _write_new(manifest, tables / "parent_round_partition_manifest.csv")
    if args.child_tables is not None and args.child_tables.resolve() != tables.resolve():
        args.child_tables.mkdir(parents=True, exist_ok=True)
        _write_new(manifest, args.child_tables / "parent_round_partition_manifest.csv")

    payload = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "parent_raw": str(args.parent_raw.resolve()),
        "parent_clustered": str(args.parent_clustered.resolve()),
        "leiden_key": args.leiden_key,
        "decision_file": str(args.decision.resolve()),
        "decision_file_sha256": _sha256(args.decision),
        "removed_h5ad": str(removed_path.resolve()),
        "child_raw": str(child_path.resolve()),
        "n_parent": int(raw.n_obs),
        "n_removed": int(removed.n_obs),
        "n_retained": int(retained.n_obs),
        "removed_clusters": decisions.loc[decisions["decision"].eq("DELETE")].to_dict(
            orient="records"
        ),
        "kept_clusters": decisions.loc[decisions["decision"].eq("KEEP"), "cluster_id"].tolist(),
    }
    decision_json = tables / "parent_removal_decision.json"
    if decision_json.exists():
        raise FileExistsError(f"Refusing to overwrite: {decision_json}")
    decision_json.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    ids = "_".join(sorted(delete_ids, key=lambda x: int(x) if x.isdigit() else x))
    print(
        f"Prepared child raw checkpoint with {retained.n_obs:,} retained cells "
        f"(removed clusters {ids}). Run run-leiden on the child; "
        "no old embedding is reused."
    )
    return 0


def _write_new(frame: pd.DataFrame, path: Path) -> None:
    _require_new(path)
    frame.to_csv(path, index=False, encoding="utf-8-sig")


def _infer_child_name(child_raw: Path) -> str:
    for part in child_raw.parts:
        if re.match(r"round_\d+", part):
            return part
    return child_raw.parent.name


if __name__ == "__main__":
    raise SystemExit(main())

