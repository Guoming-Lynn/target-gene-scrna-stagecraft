#!/usr/bin/env python3
"""Part 4 subtype TARGET_GENE catalog + identifiability forecast.

Read-only. Reuses Part 2 geometry with stem prefix F04.

Usage:
    python scripts/part4_figures.py locked_subtype.h5ad \\
        --gene SYMBOL --group subtype --colors colors.yaml --out figures/part4
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from part2_figures import render_catalog  # noqa: E402
from part4_identifiability import identifiability_table, plot_identifiability  # noqa: E402
from plotting_style import apply_publication_style, load_color_map, load_plotting_config  # noqa: E402

try:
    import scanpy as sc
except ImportError:  # pragma: no cover
    sc = None


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("h5ad")
    parser.add_argument("--gene", required=True)
    parser.add_argument("--group", default="subtype")
    parser.add_argument("--colors", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--config", default=str(SCRIPT_DIR / "plotting_config.yaml"))
    parser.add_argument("--min-cells", type=int, default=20)
    parser.add_argument("--donor-key", default="donor_id")
    parser.add_argument("--dataset-key", default="dataset")
    parser.add_argument("--source-key", default="")
    parser.add_argument("--umap-title", default="Locked subtypes")
    args = parser.parse_args(argv)
    if sc is None:
        raise SystemExit("scanpy required")

    adata = sc.read_h5ad(args.h5ad)
    if args.group not in adata.obs:
        raise SystemExit(f"missing obs column {args.group}")
    if adata.obs[args.group].isna().any():
        raise SystemExit(f"{args.group} has NA labels. Part 3 is not locked.")
    if args.gene not in adata.var_names:
        raise SystemExit(f"TARGET_GENE '{args.gene}' is not in var_names. Stop.")
    if "counts" not in adata.layers or "normalized" not in adata.layers:
        raise SystemExit("Need layers['counts'] and layers['normalized'].")
    if "X_umap" not in adata.obsm:
        raise SystemExit("Need this compartment's obsm['X_umap']. Do not borrow the global embedding.")

    colors = load_color_map(args.colors)
    missing = sorted(set(adata.obs[args.group].astype(str)) - set(colors))
    if missing:
        print(
            f"WARNING: {len(missing)} subtypes missing from color YAML; "
            f"plotted gray: {missing[:12]}"
        )
    config = load_plotting_config(args.config)
    out = Path(args.out)
    render_catalog(
        adata,
        gene=args.gene,
        group_key=args.group,
        colors=colors,
        config=config,
        out_dir=out,
        min_cells=args.min_cells,
        donor_key=args.donor_key,
        dataset_key=args.dataset_key,
        stem_prefix="F04",
        umap_title=args.umap_title,
    )

    unit_csv = out / "source_data" / "unit_summary.csv"
    units = pd.read_csv(unit_csv)
    source_key = args.source_key.strip() or args.dataset_key
    if source_key not in units.columns:
        raise SystemExit(f"unit table missing {source_key}")
    table = identifiability_table(units, args.group, source_key=source_key, min_cells=args.min_cells)
    ident_csv = out / "source_data" / "identifiability.csv"
    table.to_csv(ident_csv, index=False, encoding="utf-8-sig")
    apply_publication_style(config)
    plot_identifiability(
        table,
        args.group,
        config,
        out / f"F04_15_{args.gene}_identifiability",
        gene=args.gene,
        input_files=[unit_csv],
    )
    print(f"wrote Part 4 catalog under {out}")
    print("Identifiability flags are Part 5 forecasts, not fitted verdicts.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

