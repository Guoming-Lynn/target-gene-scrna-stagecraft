#!/usr/bin/env python3
"""Validate the required scientific metadata for one formal figure."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

REQUIRED = {
    "estimand",
    "biological_unit",
    "n_definition",
    "n",
    "model_or_test",
    "effect_scale",
    "error_bar",
    "multiple_testing_family",
    "adjustment",
    "claim_ceiling",
}


def validate_manifest(obj: object) -> tuple[list[str], list[str]]:
    """Return missing and empty required fields for one statistics manifest."""
    if not isinstance(obj, dict):
        return sorted(REQUIRED), []
    missing = sorted(REQUIRED - set(obj))
    empty = sorted(key for key in REQUIRED if key in obj and (obj[key] is None or obj[key] == ""))
    return missing, empty


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    args = parser.parse_args(argv)
    obj = json.loads(args.manifest.read_text(encoding="utf-8"))
    missing, empty = validate_manifest(obj)
    if missing or empty:
        raise SystemExit(f"invalid figure manifest: missing={missing} empty={sorted(empty)}")
    if obj["error_bar"] not in {"none", "SD", "SEM", "IQR", "95% CI", "other"}:
        raise SystemExit("invalid error_bar")
    print(f"OK {args.manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
