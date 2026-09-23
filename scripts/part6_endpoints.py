#!/usr/bin/env python3
"""Frozen endpoint membership and coverage gates for Part 6.

Drops TARGET_GENE from every set before scoring. Primary failure stops
the chapter. Support failure is NOT_ESTIMABLE for that endpoint only.

Usage:
    python scripts/part6_endpoints.py \\
        --sets 00_protocol_manifest/endpoints.yaml \\
        --model-genes 00_input_audit/model_visible_genes.txt \\
        --target TARGET_GENE \\
        --out 00_input_audit/endpoint_coverage.csv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Iterable

import pandas as pd

_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
from stagecraft.io import ensure_repo_on_path as _ensure_repo_on_path  # noqa: E402

_ensure_repo_on_path(__file__)

from stagecraft import EXIT_GATE, EXIT_OK, stop  # noqa: E402
from stagecraft.hashing import sha256_file
from stagecraft.io import load_yaml, publish_new_files


def read_members(path: Path) -> list[str]:
    genes: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        token = line.split("#", 1)[0].strip()
        if token:
            genes.append(token)
    if not genes:
        raise SystemExit(f"empty gene set: {path}")
    return genes


def coverage_row(
    *,
    endpoint_id: str,
    role: str,
    members: Iterable[str],
    model_visible: set[str],
    target: str,
    min_frac: float,
    min_genes: int,
    path: str,
    sha256: str,
) -> dict[str, object]:
    frozen = [g for g in members if g]
    if len(frozen) != len(set(frozen)):
        raise SystemExit("Frozen endpoint contains duplicate genes.")
    dropped = [g for g in frozen if g == target]
    usable = [g for g in frozen if g != target]
    visible = [g for g in usable if g in model_visible]
    n_frozen = len(usable)
    n_visible = len(visible)
    frac = (n_visible / n_frozen) if n_frozen else 0.0
    passed = (frac >= min_frac) and (n_visible >= min_genes)
    if role == "primary":
        status = "PASS" if passed else "STOPPED"
    else:
        status = "PASS" if passed else "NOT_ESTIMABLE"
    return {
        "endpoint": endpoint_id,
        "role": role,
        "path": path,
        "sha256": sha256,
        "n_frozen_including_target": len(frozen),
        "n_target_dropped": len(dropped),
        "n_frozen": n_frozen,
        "n_model_visible": n_visible,
        "coverage_frac": frac,
        "min_frac": min_frac,
        "min_genes": min_genes,
        "status": status,
        "target_dropped": ",".join(dropped),
        "visible_genes": ",".join(visible),
    }


def evaluate_sets(
    specs: list[dict],
    model_visible: set[str],
    target: str,
) -> pd.DataFrame:
    rows = []
    for spec in specs:
        path = Path(str(spec["path"]))
        members = read_members(path)
        digest = sha256_file(path)
        expected = str(spec.get("sha256") or "").strip()
        if expected and expected.lower() != "replace" and digest != expected:
            stop(f"hash mismatch for {path}: {digest} != {expected}", EXIT_GATE)
        rows.append(
            coverage_row(
                endpoint_id=str(spec["id"]),
                role=str(spec.get("role") or "primary"),
                members=members,
                model_visible=model_visible,
                target=target,
                min_frac=float(spec.get("min_frac", 0.70)),
                min_genes=int(spec.get("min_genes", 10)),
                path=str(path),
                sha256=digest,
            )
        )
    return pd.DataFrame(rows)


def primary_blocked(table: pd.DataFrame) -> bool:
    prim = table.loc[table["role"].astype(str).eq("primary")]
    return bool((prim["status"] != "PASS").any()) if len(prim) else True


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sets", type=Path, required=True)
    parser.add_argument("--model-genes", type=Path, required=True)
    parser.add_argument("--target", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)
    raw = load_yaml(args.sets) or {}
    specs: list[dict] = []
    for role in ("primary", "support"):
        for item in raw.get(role) or []:
            item = dict(item)
            item["role"] = role
            specs.append(item)
    if raw.get("direction_control"):
        item = dict(raw["direction_control"])
        item["role"] = "direction_control"
        specs.append(item)
    visible = {line.strip() for line in args.model_genes.read_text(encoding="utf-8").splitlines() if line.strip()}
    table = evaluate_sets(specs, visible, args.target)
    members_out = args.out.with_name(args.out.stem + "_members.csv")
    if args.out.resolve() == members_out.resolve():
        raise SystemExit("Coverage and members output paths must differ")
    long_rows = []
    for row in table.itertuples(index=False):
        for gene in str(row.visible_genes).split(",") if row.visible_genes else []:
            long_rows.append({"endpoint": row.endpoint, "role": row.role, "gene": gene})
    members = pd.DataFrame(long_rows, columns=["endpoint", "role", "gene"])

    def write(partials: list[Path]) -> None:
        table.drop(columns=["visible_genes"]).to_csv(partials[0], index=False)
        members.to_csv(partials[1], index=False)

    publish_new_files([args.out, members_out], write)
    if primary_blocked(table):
        print("PRIMARY COVERAGE FAILED - chapter STOPPED")
        return EXIT_GATE
    print(f"wrote {args.out}")
    return EXIT_OK


if __name__ == "__main__":
    raise SystemExit(main())

