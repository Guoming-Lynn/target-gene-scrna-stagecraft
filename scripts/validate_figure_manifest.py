#!/usr/bin/env python3
"""Validate the required scientific metadata for one formal figure."""
from __future__ import annotations
import argparse, json
from pathlib import Path

REQUIRED = {"estimand", "biological_unit", "n_definition", "n", "model_or_test",
            "effect_scale", "error_bar", "multiple_testing_family", "adjustment",
            "claim_ceiling"}

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("manifest",type=Path); a=p.parse_args(argv)
    obj=json.loads(a.manifest.read_text(encoding="utf-8"))
    missing=sorted(REQUIRED-set(obj))
    empty=[k for k in REQUIRED if k in obj and (obj[k] is None or obj[k]=="")]
    if missing or empty:
        raise SystemExit(f"invalid figure manifest: missing={missing} empty={sorted(empty)}")
    if obj["error_bar"] not in {"none","SD","SEM","IQR","95% CI","other"}:
        raise SystemExit("invalid error_bar")
    print(f"OK {a.manifest}")
    return 0
if __name__ == "__main__": raise SystemExit(main())

