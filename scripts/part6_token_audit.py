#!/usr/bin/env python3
"""Classify KO / OE eligibility from a token ledger.

Illegal mismatches stop the chapter. Cap truncation is technical
provenance, not a biological subgroup.

Required columns:
    cell_id, raw_count, final_token_present, sequence_length, pretruncation_rank

Usage:
    python scripts/part6_token_audit.py ledger.csv --out 02_tables/token_audit.csv
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import numpy as np

ILLEGAL = "ILLEGAL_MISMATCH"
TRUNCATION = "CAP_TRUNCATION"
KO_OK = "KO_ELIGIBLE"
OE_ONLY = "OE_PRIMARY_ONLY"


def classify_row(
    raw_count: float,
    token_present: bool,
    sequence_length: int,
    pretruncation_rank: float | None,
    *,
    max_length: int = 4096,
    max_gene_tokens: int = 4094,
) -> str:
    raw_pos = float(raw_count) > 0
    if raw_pos and token_present:
        return KO_OK
    if (not raw_pos) and token_present:
        return ILLEGAL
    if raw_pos and not token_present:
        rank = float("nan") if pretruncation_rank is None else float(pretruncation_rank)
        at_cap = int(sequence_length) == int(max_length)
        rank_beyond = rank == rank and rank > max_gene_tokens
        if at_cap and rank_beyond:
            return TRUNCATION
        return ILLEGAL
    return OE_ONLY


def audit_ledger(
    frame: pd.DataFrame,
    *,
    max_length: int = 4096,
    max_gene_tokens: int = 4094,
) -> pd.DataFrame:
    required = {"cell_id", "raw_count", "final_token_present", "sequence_length", "pretruncation_rank"}
    missing = required.difference(frame.columns)
    if missing:
        raise SystemExit(f"token ledger missing columns: {sorted(missing)}")
    out = frame.copy()
    if out["cell_id"].isna().any() or out["cell_id"].duplicated().any():
        raise SystemExit("Token ledger needs unique nonmissing cell_id.")
    present = out["final_token_present"].astype(str).str.lower().map({"true": True, "false": False, "1": True, "0": False})
    if present.isna().any():
        raise SystemExit("final_token_present must be true/false or 1/0.")
    for col in ("raw_count", "sequence_length"):
        values = pd.to_numeric(out[col], errors="coerce").to_numpy(float)
        if not np.isfinite(values).all() or (values < 0).any() or (values != np.floor(values)).any():
            raise SystemExit(f"{col} must contain finite nonnegative integers.")
    if (out["sequence_length"].astype(float) > max_length).any():
        raise SystemExit("sequence_length exceeds the token cap.")
    ranks = pd.to_numeric(out["pretruncation_rank"], errors="coerce")
    known = ranks.dropna().to_numpy(float)
    if (out["pretruncation_rank"].notna() & ranks.isna()).any() or not np.isfinite(known).all() or (known < 1).any() or (known != np.floor(known)).any():
        raise SystemExit("Known pretruncation ranks must be finite positive integers.")
    out["token_class"] = [
        classify_row(
            raw,
            bool(tok),
            int(length),
            None if pd.isna(rank) else float(rank),
            max_length=max_length,
            max_gene_tokens=max_gene_tokens,
        )
        for raw, tok, length, rank in zip(
            out["raw_count"], present, out["sequence_length"], out["pretruncation_rank"]
        )
    ]
    out["ko_eligible"] = out["token_class"].eq(KO_OK)
    out["oe_primary_eligible"] = out["token_class"].ne(ILLEGAL)
    out["oe_symmetry_eligible"] = out["token_class"].eq(KO_OK)
    out["illegal"] = out["token_class"].eq(ILLEGAL)
    return out


def summary(frame: pd.DataFrame) -> dict[str, int]:
    return {
        "n_cells": int(len(frame)),
        "n_ko_eligible": int(frame["ko_eligible"].sum()),
        "n_oe_primary_eligible": int(frame["oe_primary_eligible"].sum()),
        "n_truncation": int(frame["token_class"].eq(TRUNCATION).sum()),
        "n_illegal": int(frame["illegal"].sum()),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ledger", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--max-length", type=int, default=4096)
    parser.add_argument("--max-gene-tokens", type=int, default=4094)
    args = parser.parse_args(argv)
    table = audit_ledger(
        pd.read_csv(args.ledger),
        max_length=args.max_length,
        max_gene_tokens=args.max_gene_tokens,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    if args.out.exists():
        raise SystemExit(f"Refusing to overwrite: {args.out}")
    table.to_csv(args.out, index=False)
    stats = summary(table)
    print(stats)
    if stats["n_illegal"]:
        print("ILLEGAL token mismatch — chapter STOPPED")
        return 2
    if stats["n_ko_eligible"] == 0:
        print("TOKEN_UNOBSERVABLE — no KO-eligible cell")
        return 3
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

