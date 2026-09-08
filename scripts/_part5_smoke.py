#!/usr/bin/env python3
"""Smoke tests for Part 5 helpers that do not need scanpy or R."""
from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from part5_eligibility import design_matrix, evaluate, residual_df  # noqa: E402
from part5_pseudobulk import jeffreys_per_10pct, zscore_full  # noqa: E402
from part5_source_blocks import between_within_fraction, suggest_siblings  # noqa: E402
from part5_verdict import apply_table  # noqa: E402


def test_jeffreys() -> None:
    # 5 / 10 detected -> (5.5 / 11) * 10 ≈ 5.0
    got = float(jeffreys_per_10pct(np.array([5]), np.array([10]))[0])
    assert abs(got - 5.0) < 1e-12, got


def test_design_rdf() -> None:
    frame = pd.DataFrame(
        {
            "dataset": ["A"] * 6 + ["B"] * 6,
            "z": np.array([0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1], dtype=float),
            "exposure": np.linspace(0, 2, 12),
        }
    )
    matrix, names = design_matrix(frame, ["dataset"], ["z", "exposure"])
    n, rank, rdf = residual_df(matrix)
    assert n == 12
    assert rank == 4  # intercept + 1 dummy + z + exposure
    assert rdf == 8
    assert "exposure" in names


def test_verdict_first_hit() -> None:
    rows = [
        {"token": "NOT_ESTIMABLE", "when": {"full_status": "NOT_ESTIMABLE"}},
        {"token": "SINGLE_SOURCE_DEPENDENT", "when": {"drop_dominant_status": "NOT_ESTIMABLE"}},
        {"token": "FROZEN_PASS", "when": {"full_status": "SUCCESS"}},
        {"token": "INCONCLUSIVE", "when": {}},
    ]
    hit = apply_table(rows, {"full_status": "SUCCESS", "drop_dominant_status": "NOT_ESTIMABLE"})
    assert hit["verdict"] == "SINGLE_SOURCE_DEPENDENT"
    assert hit["matched_row"] == 2


def test_sibling_suggestion() -> None:
    meta = pd.DataFrame(
        {
            "dataset": ["GSE000001"] * 3 + ["GSE000002"] * 3,
            "donor_id": ["DONOR_5", "DONOR_7", "DONOR_8", "DONOR_4", "DONOR_6", "DONOR_9"],
        }
    )
    sug = suggest_siblings(meta, "dataset", "donor_id")
    assert len(sug) == 1
    assert bool(sug.iloc[0]["interleaved_donor_numbers"])
    assert sug.iloc[0]["suggestion"] == "SAME_STUDY_LIKELY"


def test_evaluate_dead_arm() -> None:
    meta = pd.DataFrame(
        {
            "subtype": ["A"] * 6,
            "dataset": ["X"] * 6,
            "dataset_donor_id": [f"d{i}" for i in range(6)],
            "jeffreys_per_10pct": np.ones(6),
            "z_log1p_n_cells": np.zeros(6),
            "eligible": True,
            "n_cells": 30,
            "source_block": ["X"] * 6,
        }
    )
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "el.csv"
        table = evaluate(
            meta,
            arm_key="subtype",
            categorical=["subtype", "dataset"],
            numeric=["z_log1p_n_cells", "jeffreys_per_10pct"],
            exposure="jeffreys_per_10pct",
            dataset_key="dataset",
            source_key="source_block",
            n_formal=12,
            n_exploratory=8,
            min_datasets_formal=3,
            min_rdf_formal=6,
            min_rdf_fit=4,
            sd_floor=0.5,
            drop_single=True,
            part4=None,
        )
        assert (table["status"] == "NOT_ESTIMABLE").any(), table
        path.write_text(table.to_csv(index=False), encoding="utf-8")
        json.dumps(table.to_dict(orient="records"))


def test_variance_share() -> None:
    y = np.array([0.0, 0.1, 0.0, 0.1, 10.0, 10.1, 10.0, 10.1])
    g = np.array(["a", "a", "a", "a", "b", "b", "b", "b"])
    share = between_within_fraction(y, g)
    assert share["between_fraction"] > 0.9, share


if __name__ == "__main__":
    test_jeffreys()
    test_design_rdf()
    test_verdict_first_hit()
    test_sibling_suggestion()
    test_evaluate_dead_arm()
    test_variance_share()
    z = zscore_full(np.array([1.0, 2.0, 3.0]))
    assert abs(z.mean()) < 1e-12
    print("part5 smoke ok")


