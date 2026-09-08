#!/usr/bin/env python3
"""Leave-one-donor-out embedding axes and Δaxis.

The held-out donor must not enter its own reference axis. Perturbed CLS
must not enter axis construction.

Usage (after original CLS is cached):
    python scripts/part6_axes.py \\
        --cls original_cls.npy --cells cell_order.csv \\
        --scores baseline_scores.csv --endpoint ENDPOINT_A \\
        --out 03_geneformer/axes/ENDPOINT_A.npz
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def l2_normalize(vectors: np.ndarray, axis: int = -1, eps: float = 1e-12) -> np.ndarray:
    norms = np.linalg.norm(vectors, axis=axis, keepdims=True)
    return vectors / np.maximum(norms, eps)


def top_bottom_indices(
    scores: np.ndarray,
    cell_ids: np.ndarray,
    *,
    quartile: float = 0.25,
    min_side: int = 3,
) -> tuple[np.ndarray, np.ndarray]:
    n = int(len(scores))
    if not np.isfinite(scores).all() or len(np.unique(scores)) < 2:
        raise ValueError("Finite, nonconstant baseline scores required")
    k = max(min_side, int(np.ceil(quartile * n)))
    if n < 2 * k:
        raise ValueError(f"need ≥{2 * k} cells to form top/bottom, got {n}")
    order = np.lexsort((np.asarray(cell_ids, dtype=str), np.asarray(scores, dtype=float)))
    return order[-k:], order[:k]


def donor_axis(cls: np.ndarray, top_idx: np.ndarray, bottom_idx: np.ndarray) -> np.ndarray:
    if not np.isfinite(cls).all():
        raise ValueError("Non-finite CLS")
    top = l2_normalize(np.asarray(cls)[top_idx], axis=1).mean(axis=0)
    bottom = l2_normalize(np.asarray(cls)[bottom_idx], axis=1).mean(axis=0)
    if np.linalg.norm(top-bottom) <= 1e-12:
        raise ValueError("Zero donor axis")
    return l2_normalize(top - bottom, axis=0).ravel()


def loo_axis(donor_axes: dict[str, np.ndarray], held_out: str) -> np.ndarray:
    others = [vec for key, vec in donor_axes.items() if key != held_out]
    if not others:
        raise ValueError(f"no training axes left after holding out {held_out}")
    stacked = np.stack(others, axis=0)
    if np.linalg.norm(stacked.mean(axis=0)) <= 1e-12:
        raise ValueError("Training donor axes cancel")
    return l2_normalize(stacked.mean(axis=0), axis=0).ravel()


def delta_axis(cls_perturbed: np.ndarray, cls_original: np.ndarray, axis: np.ndarray) -> np.ndarray:
    shift = l2_normalize(np.asarray(cls_perturbed), axis=-1) - l2_normalize(np.asarray(cls_original), axis=-1)
    axis = np.asarray(axis).ravel()
    if shift.ndim == 1:
        return np.array([float(np.dot(shift, axis))])
    return shift @ axis


def assert_no_self_leakage(donor_axes: dict[str, np.ndarray], held_out: str, loo: np.ndarray) -> None:
    """Check the normalized mean of training donors; coincident directions are valid."""
    if held_out not in donor_axes:
        raise AssertionError(f"{held_out} has no donor axis")
    if not np.allclose(loo, loo_axis(donor_axes, held_out)):
        raise AssertionError("LOO axis differs from the training-only reference")


def build_axes(
    cls: np.ndarray,
    cells: pd.DataFrame,
    scores: pd.Series,
    *,
    quartile: float = 0.25,
    min_side: int = 3,
    min_cells: int = 10,
    min_training: int = 10,
) -> dict[str, np.ndarray]:
    required = {"cell_id", "dataset_donor_id"}
    missing = required.difference(cells.columns)
    if missing:
        raise SystemExit(f"cell table missing: {sorted(missing)}")
    if len(cls) != len(cells):
        raise SystemExit("CLS rows must match cell table rows, same order")
    if np.asarray(cls).ndim != 2 or not np.isfinite(cls).all():
        raise SystemExit("CLS must be a finite two-dimensional matrix.")
    if cells[list(required)].isna().any().any() or cells["cell_id"].duplicated().any():
        raise SystemExit("Cell identities must be nonmissing and cell_id unique.")
    cells = cells.reset_index(drop=True)
    scores = pd.Series(np.asarray(scores), index=cells.index)
    donor_axes: dict[str, np.ndarray] = {}
    donors = cells["dataset_donor_id"].astype(str).to_numpy()
    cell_ids = cells["cell_id"].to_numpy()
    score_values = scores.to_numpy()
    for donor in sorted(set(donors)):
        index = np.flatnonzero(donors == donor)
        if len(index) < min_cells:
            continue
        try:
            top, bottom = top_bottom_indices(
                score_values[index],
                cell_ids[index],
                quartile=quartile,
                min_side=min_side,
            )
        except ValueError:
            continue
        try:
            donor_axes[str(donor)] = donor_axis(cls[index], top, bottom)
        except ValueError:
            continue
    if len(donor_axes) < min_training + 1:
        raise SystemExit(f"need ≥{min_training + 1} valid donor axes, got {len(donor_axes)}")
    loo = {donor: loo_axis(donor_axes, donor) for donor in donor_axes}
    for donor, axis in loo.items():
        assert_no_self_leakage(donor_axes, donor, axis)
    return loo


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cls", type=Path, required=True)
    parser.add_argument("--cells", type=Path, required=True)
    parser.add_argument("--scores", type=Path, required=True)
    parser.add_argument("--endpoint", required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--score-col", default=None)
    args = parser.parse_args(argv)
    cls = np.load(args.cls)
    cells = pd.read_csv(args.cells)
    score_table = pd.read_csv(args.scores)
    col = args.score_col or args.endpoint
    if col not in score_table.columns:
        raise SystemExit(f"score column {col} missing")
    if "cell_id" not in score_table.columns:
        raise SystemExit("Scores require cell_id for verified alignment.")
    if set(cells["cell_id"]) != set(score_table["cell_id"]):
        raise SystemExit("Score and CLS cell identities must match exactly.")
    score_table = cells[["cell_id"]].merge(score_table, on="cell_id", how="left", validate="one_to_one")
    loo = build_axes(cls, cells, score_table[col])
    args.out.parent.mkdir(parents=True, exist_ok=True)
    if args.out.exists():
        raise SystemExit(f"Refusing to overwrite: {args.out}")
    np.savez(args.out, **{key: value for key, value in loo.items()})
    print(f"n_loo_axes={len(loo)} wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

