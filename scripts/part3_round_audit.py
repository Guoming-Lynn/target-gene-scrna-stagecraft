#!/usr/bin/env python3
"""Parent vs child Leiden membership after a Part 3 HVG recompute.

New IDs are not the old IDs. This table is the source for the membership
heatmap; it is not a mapping YAML.

Usage:
    python scripts/part3_round_audit.py \\
        --parent parent_leiden.h5ad --parent-key leiden_r0_5 \\
        --child child_leiden.h5ad --child-key leiden_r0_5 \\
        --out tables/round_01_after_removal/old_to_new_membership.csv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
from stagecraft.io import ensure_repo_on_path as _ensure_repo_on_path  # noqa: E402

_ensure_repo_on_path(__file__)

from stagecraft.io import CSV_EXCEL, publish_new_files  # noqa: E402

try:
    import scanpy as sc
except ImportError:  # pragma: no cover
    sc = None


def membership_counts(
    parent_labels: pd.Series,
    child_labels: pd.Series,
    *,
    min_shared_fraction: float = 0.5,
) -> tuple[pd.DataFrame, list[str], list[str]]:
    parent_labels = parent_labels.astype(str)
    child_labels = child_labels.astype(str)
    shared = parent_labels.index.intersection(child_labels.index)
    if len(shared) == 0:
        raise SystemExit("No shared barcodes between parent and child objects.")
    fraction = len(shared) / len(child_labels)
    if fraction < min_shared_fraction:
        raise SystemExit(
            f"Shared barcodes {len(shared)} / {len(child_labels)} child barcodes "
            f"is below {min_shared_fraction}"
        )
    frame = pd.DataFrame(
        {
            "parent_cluster": parent_labels.loc[shared].to_numpy(),
            "child_cluster": child_labels.loc[shared].to_numpy(),
        }
    )
    counts = (
        frame.groupby(["parent_cluster", "child_cluster"], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    parent_n = parent_labels.loc[shared].value_counts().rename("parent_n_shared")
    child_n = child_labels.loc[shared].value_counts().rename("child_n_shared")
    counts = counts.merge(parent_n, left_on="parent_cluster", right_index=True, how="left")
    counts = counts.merge(child_n, left_on="child_cluster", right_index=True, how="left")
    counts["fraction_of_parent"] = counts["n_cells"] / counts["parent_n_shared"]
    counts["fraction_of_child"] = counts["n_cells"] / counts["child_n_shared"]
    only_parent = sorted(set(parent_labels.index) - set(child_labels.index))
    only_child = sorted(set(child_labels.index) - set(parent_labels.index))
    return counts, only_parent, only_child


def _sort_id(series: pd.Series) -> pd.Series:
    return series.map(lambda x: (0, int(x)) if str(x).isdigit() else (1, str(x)))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parent", type=Path, required=True)
    parser.add_argument("--parent-key", required=True)
    parser.add_argument("--child", type=Path, required=True)
    parser.add_argument("--child-key", required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--min-shared-fraction", type=float, default=0.5)
    args = parser.parse_args(argv)
    if sc is None:
        raise SystemExit("scanpy required")

    parent = sc.read_h5ad(args.parent, backed="r")
    child = sc.read_h5ad(args.child, backed="r")
    if args.parent_key not in parent.obs:
        raise SystemExit(f"parent missing {args.parent_key}")
    if args.child_key not in child.obs:
        raise SystemExit(f"child missing {args.child_key}")

    counts, only_parent, only_child = membership_counts(
        parent.obs[args.parent_key],
        child.obs[args.child_key],
        min_shared_fraction=args.min_shared_fraction,
    )
    counts = counts.sort_values(["parent_cluster", "child_cluster"], key=_sort_id)
    orphan = args.out.with_name(args.out.stem + "_barcode_orphans.csv")
    orphans = pd.DataFrame(
        {
            "barcode": only_parent + only_child,
            "where": (["parent_only"] * len(only_parent)) + (["child_only"] * len(only_child)),
        }
    )

    def write(partials: list[Path]) -> None:
        counts.to_csv(partials[0], index=False, encoding=CSV_EXCEL)
        orphans.to_csv(partials[1], index=False, encoding=CSV_EXCEL)

    publish_new_files([args.out, orphan], write)
    print(
        f"wrote {args.out} ({len(counts)} parent×child cells; "
        f"{len(only_parent)} parent-only barcodes, {len(only_child)} child-only)"
    )
    print("New Leiden IDs are not the old IDs.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

