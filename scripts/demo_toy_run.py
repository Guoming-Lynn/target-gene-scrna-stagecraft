#!/usr/bin/env python3
"""Run helper CLIs on synthetic toy data. Not a biological analysis."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]
_SCRIPTS = Path(__file__).resolve().parent
for path in (_ROOT, _SCRIPTS):
    text = str(path)
    if text not in sys.path:
        sys.path.insert(0, text)

from stagecraft import EXIT_OK  # noqa: E402
from generate_toy_data import main as write_toy  # noqa: E402
from part5_eligibility import main as eligibility_main  # noqa: E402
from part5_pseudobulk import main as pseudobulk_main  # noqa: E402
from part5_source_blocks import main as source_main  # noqa: E402
from simulation_contract import main as simulation_main  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=7)
    args = parser.parse_args(argv)
    out = args.out
    out.mkdir(parents=True, exist_ok=False)
    toy = out / "toy.h5ad"
    write_toy(["--out", str(toy), "--seed", str(args.seed)])
    pb = out / "03_pseudobulk"
    code = pseudobulk_main(
        [str(toy), "--gene", "TARGET_FEATURE", "--group", "cell_type", "--out", str(pb)]
    )
    if code:
        return code
    tables = out / "02_tables"
    code = source_main([str(pb / "metadata.csv"), "--out", str(tables)])
    if code:
        return code
    code = eligibility_main(
        [
            str(tables / "metadata_with_source_block.csv"),
            "--arm-key",
            "cell_type",
            "--categorical",
            "cell_type,dataset",
            "--out",
            str(tables / "eligibility.csv"),
        ]
    )
    if code:
        return code
    code = simulation_main(
        ["--seed", str(args.seed), "--replicates", "100", "--out", str(out / "simulation_manifest.json")]
    )
    if code:
        return code
    report = {
        "status": "PASS",
        "formal_analysis": False,
        "toy": str(toy),
        "artifacts": [
            "toy.h5ad",
            "03_pseudobulk/metadata.csv",
            "02_tables/eligibility.csv",
            "simulation_manifest.json",
        ],
        "note": "Synthetic helper path only. No biological result.",
    }
    (out / "demo_report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"wrote toy demo: {out}")
    return EXIT_OK


if __name__ == "__main__":
    raise SystemExit(main())
