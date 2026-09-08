#!/usr/bin/env python3
"""Judge Geneformer wrapper-vs-official smoke parity. Does not run the model.

Required columns on --parity:
    perturbation, max_abs_cosine_diff, max_abs_embedding_rerun_diff

OE must be present if the protocol declared OE. A KO-only pass with a
missing OE row is a fail when OE is in --perturbations.

Usage:
    python scripts/part6_smoke_gate.py parity.csv --perturbations KO,OE
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import pandas as pd

DEFAULT_COSINE = 1e-5
DEFAULT_DETERMINISM = 1e-6


def judge_smoke(
    table: pd.DataFrame,
    *,
    perturbations: list[str],
    cosine_max: float = DEFAULT_COSINE,
    determinism_max: float = DEFAULT_DETERMINISM,
) -> dict[str, object]:
    required = {"perturbation", "max_abs_cosine_diff", "max_abs_embedding_rerun_diff"}
    missing = required.difference(table.columns)
    if missing:
        raise SystemExit(f"parity table missing: {sorted(missing)}")
    present = set(table["perturbation"].astype(str).str.upper())
    needed = {p.upper() for p in perturbations}
    if not needed or not needed.issubset({"KO", "OE"}):
        raise SystemExit("Declare a non-empty set of KO/OE perturbations.")
    if not all(math.isfinite(x) and x >= 0 for x in (cosine_max, determinism_max)):
        raise SystemExit("Smoke thresholds must be finite and non-negative.")
    missing_pert = sorted(needed - present)
    rows = []
    ok = True
    reasons = []
    if missing_pert:
        ok = False
        reasons.append("missing_perturbation:" + ",".join(missing_pert))
    for pert in needed:
        block = table.loc[table["perturbation"].astype(str).str.upper().eq(pert)]
        if block.empty:
            continue
        metrics = block[["max_abs_cosine_diff", "max_abs_embedding_rerun_diff"]].apply(pd.to_numeric, errors="coerce")
        if not all(math.isfinite(x) and x >= 0 for x in metrics.to_numpy().ravel()):
            ok = False
            reasons.append(f"{pert}_invalid_or_missing_metrics")
            continue
        cosine = float(metrics["max_abs_cosine_diff"].max())
        rerun = float(metrics["max_abs_embedding_rerun_diff"].max())
        cosine_ok = cosine <= cosine_max
        det_ok = rerun <= determinism_max
        if not cosine_ok:
            ok = False
            reasons.append(f"{pert}_cosine={cosine}>{cosine_max}")
        if not det_ok:
            ok = False
            reasons.append(f"{pert}_determinism={rerun}>{determinism_max}")
        rows.append(
            {
                "perturbation": pert,
                "max_abs_cosine_diff": cosine,
                "max_abs_embedding_rerun_diff": rerun,
                "cosine_pass": cosine_ok,
                "determinism_pass": det_ok,
            }
        )
    return {
        "pass": ok,
        "cosine_max": cosine_max,
        "determinism_max": determinism_max,
        "required_perturbations": sorted(needed),
        "reasons": reasons,
        "rows": rows,
        "note": "Do not relax the cosine gate. Do not drop OE to make smoke pass. Do not switch checkpoint.",
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("parity", type=Path)
    parser.add_argument("--perturbations", default="KO")
    parser.add_argument("--cosine-max", type=float, default=DEFAULT_COSINE)
    parser.add_argument("--determinism-max", type=float, default=DEFAULT_DETERMINISM)
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args(argv)
    perts = [p.strip() for p in args.perturbations.split(",") if p.strip()]
    verdict = judge_smoke(
        pd.read_csv(args.parity),
        perturbations=perts,
        cosine_max=args.cosine_max,
        determinism_max=args.determinism_max,
    )
    text = json.dumps(verdict, indent=2)
    print(text)
    if args.out:
        if args.out.exists():
            raise SystemExit(f"Refusing to overwrite: {args.out}")
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text + "\n", encoding="utf-8")
    return 0 if verdict["pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())

