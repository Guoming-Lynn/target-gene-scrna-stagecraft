#!/usr/bin/env python3
"""Smoke tests for Part 6 helpers that do not need Geneformer or GPU."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from part6_axes import build_axes, delta_axis, loo_axis, l2_normalize  # noqa: E402
from part6_controls import select_controls  # noqa: E402
from part6_eligibility import donor_eligibility  # noqa: E402
from part6_endpoints import coverage_row, primary_blocked  # noqa: E402
from part6_sign_tests import bh_adjust_fixed_family, exact_two_sided_sign_p, sign_tests  # noqa: E402
from part6_smoke_gate import judge_smoke  # noqa: E402
from part6_token_audit import ILLEGAL, KO_OK, TRUNCATION, audit_ledger, classify_row  # noqa: E402
from part5_verdict import apply_table  # noqa: E402


def test_token_truncation_vs_illegal() -> None:
    assert classify_row(1, True, 4096, 10) == KO_OK
    assert classify_row(1, False, 4096, 4211) == TRUNCATION
    assert classify_row(1, False, 3000, 10) == ILLEGAL
    assert classify_row(0, True, 100, 1) == ILLEGAL
    ledger = pd.DataFrame(
        {
            "cell_id": ["a", "b", "c", "d"],
            "raw_count": [1, 1, 0, 0],
            "final_token_present": [True, False, False, False],
            "sequence_length": [100, 4096, 100, 100],
            "pretruncation_rank": [10, 5000, 10, 10],
        }
    )
    out = audit_ledger(ledger)
    assert int(out["ko_eligible"].sum()) == 1
    assert int(out["illegal"].sum()) == 0
    assert int(out["token_class"].eq(TRUNCATION).sum()) == 1


def test_target_dropped_from_endpoint() -> None:
    row = coverage_row(
        endpoint_id="E",
        role="primary",
        members=["A", "B", "TARGET_GENE", "C"],
        model_visible={"A", "B", "C", "TARGET_GENE"},
        target="TARGET_GENE",
        min_frac=0.70,
        min_genes=2,
        path="x",
        sha256="y",
    )
    assert row["n_target_dropped"] == 1
    assert row["n_frozen"] == 3
    assert row["n_model_visible"] == 3
    assert row["status"] == "PASS"
    failed = pd.DataFrame([dict(row, status="STOPPED", role="primary")])
    assert primary_blocked(failed) is True


def test_sign_p_all_positive() -> None:
    p = exact_two_sided_sign_p(5, 5)
    assert abs(p - 2 / 32) < 1e-12, p


def test_bh_does_not_shrink_family() -> None:
    # Three finite tests, declared family 4. Rank-1 p=0.01 → 0.01*4/1 = 0.04
    adjusted = bh_adjust_fixed_family([0.01, 0.02, 0.03, float("nan")], 4)
    assert abs(adjusted[0] - 0.04) < 1e-12, adjusted
    shrunk = bh_adjust_fixed_family([0.01, 0.02, 0.03], 3)
    assert shrunk[0] < adjusted[0]


def test_sign_family_size_in_output() -> None:
    donors = pd.DataFrame(
        {
            "dataset_donor_id": [f"d{i}" for i in range(8)] * 2,
            "target_symbol": ["G"] * 16,
            "perturbation": ["KO"] * 16,
            "endpoint": ["E1"] * 8 + ["E2"] * 8,
            "median_delta_axis": [0.1] * 7 + [-0.1] + [0.1] * 8,
            "eligible": True,
        }
    )
    table = sign_tests(donors, target="G", family_size=2, perturbations=["KO"], endpoints=["E1", "E2"])
    assert len(table) == 2
    assert int(table["family_size"].iloc[0]) == 2
    assert table["raw_p_value"].notna().all()


def test_eligibility_before_summary() -> None:
    cells = pd.DataFrame(
        {
            "cell_id": [f"c{i}" for i in range(10)],
            "dataset_donor_id": ["d1"] * 4 + ["d2"] * 6,
            "target_symbol": ["G"] * 10,
            "perturbation": ["KO"] * 10,
            "analysis_population": ["KO_token_present"] * 10,
            "endpoint": ["E"] * 10,
            "delta_axis": np.linspace(-0.1, 0.1, 10),
            "run_status": ["RUN"] * 10,
        }
    )
    ledger, eligible = donor_eligibility(cells, gates={"KO": 5, "OE": 10, "OE_SYMMETRY": 5})
    assert bool(ledger.loc[ledger["dataset_donor_id"].eq("d1"), "eligible"].iloc[0]) is False
    assert bool(ledger.loc[ledger["dataset_donor_id"].eq("d2"), "eligible"].iloc[0]) is True
    assert set(eligible["dataset_donor_id"]) == {"d2"}


def test_loo_excludes_held_out() -> None:
    rng = np.random.default_rng(0)
    n_donors, n_per, dim = 12, 12, 8
    cls = rng.normal(size=(n_donors * n_per, dim))
    # Make donor 0's high-score cells uniquely aligned to e0
    for i in range(n_per):
        cls[i] = 0
        if i >= 9:
            cls[i, 0] = 5
        elif i < 3:
            cls[i, 0] = -5
    cells = pd.DataFrame(
        {
            "cell_id": [f"c{i:03d}" for i in range(n_donors * n_per)],
            "dataset_donor_id": [f"d{i // n_per}" for i in range(n_donors * n_per)],
        }
    )
    scores = pd.Series(np.tile(np.linspace(0, 1, n_per), n_donors))
    loo = build_axes(cls, cells, scores, min_cells=10, min_side=3, min_training=10)
    own = {}
    donors = cells["dataset_donor_id"].to_numpy()
    for donor in sorted(set(donors)):
        idx = np.flatnonzero(donors == donor)
        # rebuild own axis from the same helper by isolating
        from part6_axes import donor_axis, top_bottom_indices

        top, bottom = top_bottom_indices(scores.to_numpy()[idx], cells["cell_id"].to_numpy()[idx])
        own[donor] = donor_axis(cls[idx], top, bottom)
    held = "d0"
    reconstructed = loo_axis(own, held)
    assert np.allclose(loo[held], reconstructed)
    assert not np.allclose(loo[held], own[held])


def test_delta_axis_sign() -> None:
    axis = l2_normalize(np.array([1.0, 0.0, 0.0]))
    orig = np.array([1.0, 0.0, 0.0])
    pert = np.array([0.0, 1.0, 0.0])
    # moving away from axis-aligned original toward y → negative or smaller projection along x
    value = float(delta_axis(pert, orig, axis)[0])
    assert value < 0


def test_control_window_expands_once() -> None:
    n = 32
    symbols = ["TARGET"] + [f"G{i:03d}" for i in range(1, n)]
    det = np.full(n, 0.10, dtype=float)
    mean = np.full(n, 1.0, dtype=float)
    det[0] = 0.50
    mean[0] = 2.0
    det[1:13] = 0.58
    mean[1:13] = 2.5
    mean[-1] = 100.0
    mean[-2] = 1e-4
    genes = pd.DataFrame(
        {
            "gene_symbol": symbols,
            "detection_fraction": det,
            "mean_raw_counts_per_cell": mean,
            "in_model_vocabulary": True,
        }
    )
    selected, _, window = select_controls(genes, target="TARGET", endpoint_union=set(), max_n=10)
    assert window == "expanded_once"
    assert len(selected) == 10
    assert "TARGET" not in set(selected["gene_symbol"])


def test_smoke_requires_declared_oe() -> None:
    table = pd.DataFrame(
        {
            "perturbation": ["KO"],
            "max_abs_cosine_diff": [1e-7],
            "max_abs_embedding_rerun_diff": [0.0],
        }
    )
    ko_only = judge_smoke(table, perturbations=["KO"])
    assert ko_only["pass"] is True
    both = judge_smoke(table, perturbations=["KO", "OE"])
    assert both["pass"] is False


def test_verdict_first_hit() -> None:
    rows = [
        {"token": "STOPPED", "when": {"blocker": True}},
        {"token": "SHAM_DRIFT", "when": {"sham_pass": False}},
        {"token": "EMBEDDING_SHIFT_CONSISTENT", "when": {"ko_primary_bh_pass": True}},
        {"token": "INCONCLUSIVE", "when": {}},
    ]
    hit = apply_table(rows, {"blocker": False, "sham_pass": False, "ko_primary_bh_pass": True})
    assert hit["verdict"] == "SHAM_DRIFT"


def main() -> int:
    tests = [
        test_token_truncation_vs_illegal,
        test_target_dropped_from_endpoint,
        test_sign_p_all_positive,
        test_bh_does_not_shrink_family,
        test_sign_family_size_in_output,
        test_eligibility_before_summary,
        test_loo_excludes_held_out,
        test_delta_axis_sign,
        test_control_window_expands_once,
        test_smoke_requires_declared_oe,
        test_verdict_first_hit,
    ]
    for fn in tests:
        fn()
        print("ok", fn.__name__)
    print("part6 smoke passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

