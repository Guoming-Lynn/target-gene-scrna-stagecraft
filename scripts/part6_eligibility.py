#!/usr/bin/env python3
"""Donor eligibility BEFORE across-donor summaries.

KO ≥5 token-present successes; OE ≥10 tokenized successes.
Ineligible donors stay in the ledger and do not enter medians or tests.

Usage:
    python scripts/part6_eligibility.py cell_effects.csv \\
        --out 02_tables/donor_eligibility.csv
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

DEFAULT_GATES = {
    "KO": 5,
    "OE": 10,
    "OE_SYMMETRY": 5,
}


def _min_for(perturbation: str, population: str, gates: dict[str, int]) -> int:
    pert = str(perturbation).upper()
    pop = str(population).upper()
    if "SYMMETRY" in pop:
        return int(gates.get("OE_SYMMETRY", 5))
    if pert == "KO":
        return int(gates.get("KO", 5))
    if pert == "OE":
        return int(gates.get("OE", 10))
    return int(gates.get(pert, 5))


def donor_eligibility(
    cells: pd.DataFrame,
    *,
    gates: dict[str, int] | None = None,
    near_zero: float = 1e-6,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    required = {
        "cell_id",
        "dataset_donor_id",
        "target_symbol",
        "perturbation",
        "analysis_population",
        "endpoint",
        "delta_axis",
        "run_status",
    }
    missing = required.difference(cells.columns)
    if missing:
        raise SystemExit(f"cell effects missing columns: {sorted(missing)}")
    gates = gates or DEFAULT_GATES
    run = cells.copy()
    run["success"] = run["run_status"].astype(str).eq("RUN") & np.isfinite(run["delta_axis"].astype(float))
    group_cols = ["dataset_donor_id", "target_symbol", "perturbation", "analysis_population", "endpoint"]
    if run[["cell_id", *group_cols]].isna().any().any():
        raise SystemExit("Cell effects require nonmissing identities.")
    if run.duplicated(["cell_id", *group_cols]).any():
        raise SystemExit("Duplicate cell effect within a donor/endpoint/perturbation population.")
    if not run["perturbation"].astype(str).str.upper().isin(["KO", "OE"]).all():
        raise SystemExit("Unknown perturbation; expected KO or OE.")
    rows = []
    eligible_effects = []
    for keys, block in run.groupby(group_cols, sort=True, dropna=False):
        donor, target, pert, pop, endpoint = keys
        n_success = int(block["success"].sum())
        minimum = _min_for(pert, pop, gates)
        ok = n_success >= minimum
        values = block.loc[block["success"], "delta_axis"].astype(float)
        median = float(values.median()) if len(values) else float("nan")
        q1 = float(values.quantile(0.25)) if len(values) else float("nan")
        q3 = float(values.quantile(0.75)) if len(values) else float("nan")
        reason = "" if ok else f"n_success={n_success} < {minimum}"
        rows.append(
            {
                "dataset_donor_id": donor,
                "target_symbol": target,
                "perturbation": pert,
                "analysis_population": pop,
                "endpoint": endpoint,
                "n_cells": int(len(block)),
                "n_success": n_success,
                "min_required": minimum,
                "eligible": ok,
                "reason": reason,
                "median_delta_axis": median,
                "q1_delta_axis": q1,
                "q3_delta_axis": q3,
                "iqr_delta_axis": (q3 - q1) if ok else float("nan"),
                "n_cell_near_zero": int((values.abs() <= near_zero).sum()) if len(values) else 0,
            }
        )
        if ok:
            eligible_effects.append(rows[-1])
    ledger = pd.DataFrame(rows)
    eligible = pd.DataFrame(eligible_effects, columns=ledger.columns)
    return ledger, eligible


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cells", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--ko-min", type=int, default=5)
    parser.add_argument("--oe-min", type=int, default=10)
    args = parser.parse_args(argv)
    ledger, eligible = donor_eligibility(
        pd.read_csv(args.cells),
        gates={"KO": args.ko_min, "OE": args.oe_min, "OE_SYMMETRY": args.ko_min},
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    eligible_out = args.out.with_name("donor_effects_eligible.csv")
    if args.out.resolve() == eligible_out.resolve():
        raise SystemExit("Ledger and eligible output paths must differ")
    for path in (args.out, eligible_out):
        if path.exists():
            raise SystemExit(f"Refusing to overwrite: {path}")
    ledger.to_csv(args.out, index=False)
    eligible.to_csv(eligible_out, index=False)
    print(f"n_ledger={len(ledger)} n_eligible_rows={len(eligible)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

