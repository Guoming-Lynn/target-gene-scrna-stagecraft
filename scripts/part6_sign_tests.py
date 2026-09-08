#!/usr/bin/env python3
"""Closed-family exact sign tests on donor medians.

BH uses the declared family size even when some cells are NA.
Do not shrink the family. Do not test cells.

Usage:
    python scripts/part6_sign_tests.py 02_tables/donor_effects_eligible.csv \\
        --target TARGET_GENE --family-size 2 --out 02_tables/sign_tests.csv
"""
from __future__ import annotations

import argparse
from math import comb
from pathlib import Path

import numpy as np
import pandas as pd


def median_order_statistic_interval(values, confidence=0.95):
    """Central order-statistic interval for a population median.

    Coverage is at least nominal for independent identically distributed donors;
    ties are conservative. Shared fitted axes can violate this assumption.
    Keep zeros: the median estimand includes all eligible donors. This is not
    a Hodges-Lehmann/Walsh-average interval for a pseudomedian.
    """
    values = np.sort(np.asarray(values, dtype=float))
    if values.ndim != 1 or not np.isfinite(values).all():
        raise ValueError("Expected finite one-dimensional donor effects")
    if not 0 < confidence < 1:
        raise ValueError("confidence must lie between zero and one")
    n = len(values)
    if not n:
        return float("nan"), float("nan"), float("nan"), "NO_DONORS"
    k, coverage = 0, 1.0
    tail = 0
    for candidate in range(1, (n + 1) // 2 + 1):
        tail += comb(n, candidate - 1)
        candidate_coverage = 1.0 - 2.0 * tail / (2 ** n)
        if candidate_coverage < confidence:
            break
        k, coverage = candidate, candidate_coverage
    if not k:
        return -float("inf"), float("inf"), 1.0, "UNBOUNDED_SMALL_N"
    return float(values[k - 1]), float(values[n - k]), coverage, "FINITE"


def exact_two_sided_sign_p(n_positive: int, n_nonzero: int) -> float:
    if n_nonzero <= 0:
        return float("nan")
    if n_positive < 0 or n_positive > n_nonzero:
        raise ValueError("n_positive out of range")
    denom = 2 ** n_nonzero
    lower = sum(comb(n_nonzero, k) for k in range(0, n_positive + 1)) / denom
    upper = sum(comb(n_nonzero, k) for k in range(n_positive, n_nonzero + 1)) / denom
    return float(min(1.0, 2.0 * min(lower, upper)))


def bh_adjust_fixed_family(pvalues: list[float], family_size: int) -> list[float]:
    """BH on the finite tests, multiplier = declared family size (not n_tested)."""
    if family_size < 1:
        raise ValueError("family_size must be ≥1")
    n = len(pvalues)
    adjusted = [float("nan")] * n
    finite = [(i, p) for i, p in enumerate(pvalues) if p == p]
    if not finite:
        return adjusted
    order = sorted(finite, key=lambda item: item[1])
    running = 1.0
    ranked: list[tuple[int, float]] = []
    for rank in range(len(order), 0, -1):
        index, pval = order[rank - 1]
        running = min(running, min(1.0, pval * family_size / rank))
        ranked.append((index, running))
    for index, value in ranked:
        adjusted[index] = value
    # restore original-index order already stored in `adjusted`
    return adjusted


def sign_tests(
    donor: pd.DataFrame,
    *,
    target: str,
    family_size: int,
    min_nonzero: int = 5,
    near_zero: float = 1e-6,
    perturbations: list[str] | None = None,
    endpoints: list[str] | None = None,
) -> pd.DataFrame:
    frame = donor.loc[donor["target_symbol"].astype(str).eq(target)].copy()
    if "eligible" in frame.columns:
        frame = frame.loc[frame["eligible"].astype(bool)].copy()
    endpoints = endpoints or sorted(frame["endpoint"].astype(str).unique())
    perturbations = perturbations or sorted(frame["perturbation"].astype(str).unique())
    records = []
    for endpoint in endpoints:
        for pert in perturbations:
            subset = frame.loc[frame["endpoint"].astype(str).eq(endpoint) & frame["perturbation"].astype(str).eq(pert)]
            values = subset["median_delta_axis"].astype(float)
            if not np.isfinite(values).all():
                raise SystemExit("Eligible donor effects must be finite")
            if "dataset_donor_id" not in subset or subset["dataset_donor_id"].duplicated().any():
                raise SystemExit("Require one biological donor row per endpoint/perturbation; separate populations")
            nonzero = values.loc[values.abs() > near_zero]
            n_pos = int((nonzero > 0).sum())
            n_neg = int((nonzero < 0).sum())
            n_nz = int(len(nonzero))
            ci_l, ci_r, ci_coverage, ci_status = median_order_statistic_interval(values)
            if n_nz >= min_nonzero:
                raw_p, reason = exact_two_sided_sign_p(n_pos, n_nz), ""
            else:
                raw_p, reason = float("nan"), "FEWER_THAN_MIN_NONZERO_DONOR_MEDIANS"
            records.append(
                {
                    "target_symbol": target,
                    "endpoint": endpoint,
                    "perturbation": pert,
                    "n_evaluable_donors": int(len(values)),
                    "n_nonzero_donor_medians": n_nz,
                    "n_positive_donors": n_pos,
                    "n_negative_donors": n_neg,
                    "across_donor_median": float(values.median()) if len(values) else float("nan"),
                    "median_ci_lower": ci_l,
                    "median_ci_upper": ci_r,
                    "median_ci_confidence": 0.95,
                    "median_ci_coverage_under_iid": ci_coverage,
                    "median_ci_status": ci_status,
                    "median_ci_method": "central_order_statistics_all_eligible_donors",
                    "median_ci_inference_status": "NOMINAL_SHARED_AXIS_DEPENDENCE_NOT_CALIBRATED",
                    "raw_p_value": raw_p,
                    "bh_adjusted_p_value": float("nan"),
                    "test": "exact_two_sided_binomial_sign_test",
                    "inference_status": "NOMINAL_SHARED_AXIS_DEPENDENCE_NOT_CALIBRATED",
                    "evidence_class": "exploratory_embedding_only",
                    "causal_inference": False,
                    "expression_prediction": False,
                    "reason_not_tested": reason,
                    "near_zero_threshold": near_zero,
                    "family_size": family_size,
                }
            )
    result = pd.DataFrame(records)
    declared = int(family_size)
    if len(result) > declared:
        raise SystemExit(
            f"constructed {len(result)} family cells but family_size={declared}. "
            "Do not silently enlarge the family."
        )
    if len(result) < declared:
        # Pad is not added here; BH still uses declared size (conservative).
        pass
    result["bh_adjusted_p_value"] = bh_adjust_fixed_family(
        result["raw_p_value"].astype(float).tolist(), declared
    )
    return result


def across_donor_ranks(
    eligible: pd.DataFrame,
    *,
    endpoint: str,
    perturbation: str,
    target: str,
    ko_more_positive: bool = True,
) -> pd.DataFrame:
    """Directional + absolute ranks among TARGET_GENE and matched controls.

    Not a permutation p-value. Requires one row per donor × gene.
    """
    block = eligible.loc[
        eligible["endpoint"].astype(str).eq(endpoint)
        & eligible["perturbation"].astype(str).eq(perturbation)
    ]
    if block.empty:
        return pd.DataFrame()
    across = (
        block.groupby("target_symbol", sort=True)["median_delta_axis"]
        .median()
        .rename("across_donor_median")
        .reset_index()
    )
    signed = across["across_donor_median"].astype(float)
    if str(perturbation).upper() == "KO" and ko_more_positive:
        direction = signed
    else:
        direction = -signed
    across["directional_rank"] = direction.rank(ascending=False, method="min").astype(int)
    across["absolute_rank"] = signed.abs().rank(ascending=False, method="min").astype(int)
    across["n_genes"] = len(across)
    across["is_target"] = across["target_symbol"].astype(str).eq(target)
    across["endpoint"] = endpoint
    across["perturbation"] = perturbation
    across["note"] = "ranks among matched genes only; not a permutation p-value"
    return across


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("donors", type=Path)
    parser.add_argument("--target", required=True)
    parser.add_argument("--family-size", type=int, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--endpoints", default="")
    parser.add_argument("--perturbations", default="KO")
    args = parser.parse_args(argv)
    endpoints = [item.strip() for item in args.endpoints.split(",") if item.strip()] or None
    perts = [item.strip() for item in args.perturbations.split(",") if item.strip()]
    table = sign_tests(
        pd.read_csv(args.donors),
        target=args.target,
        family_size=args.family_size,
        perturbations=perts,
        endpoints=endpoints,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    if args.out.exists():
        raise SystemExit(f"Refusing to overwrite: {args.out}")
    table.to_csv(args.out, index=False)
    n_run = int(table["raw_p_value"].notna().sum())
    print(f"n_family={args.family_size} n_tested={n_run} (family size not shrunk)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

