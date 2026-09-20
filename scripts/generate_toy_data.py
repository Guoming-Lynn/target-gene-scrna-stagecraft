"""Generate a small, biologically neutral AnnData fixture for helper demos."""
from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


def build_toy(*, seed: int, n_donors: int, n_datasets: int, n_types: int, cells_per_unit: int):
    rng = np.random.default_rng(seed)
    n_cells = n_donors * n_datasets * n_types * cells_per_unit
    n_genes = 50
    X = rng.poisson(1.0, (n_cells, n_genes)).astype("int32")
    X[:, 0] += rng.poisson(2, n_cells)
    donors, datasets, types = [], [], []
    for dataset in range(n_datasets):
        for donor in range(n_donors):
            for cell_type in range(n_types):
                donors.extend([f"d{donor:02d}"] * cells_per_unit)
                datasets.extend([f"s{dataset}"] * cells_per_unit)
                types.extend([f"t{cell_type}"] * cells_per_unit)
    obs = pd.DataFrame(
        {"donor_id": donors, "dataset": datasets, "cell_type": types},
        index=[f"cell{i:04d}" for i in range(n_cells)],
    )
    var = pd.DataFrame(index=["TARGET_FEATURE"] + [f"G{i:03d}" for i in range(n_genes - 1)])
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.layers["counts"] = X.copy()
    adata.layers["normalized"] = np.log1p(X.astype("float64"))
    adata.obsm["X_umap"] = rng.normal(size=(n_cells, 2))
    return adata


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--n-donors", type=int, default=8)
    parser.add_argument("--n-datasets", type=int, default=2)
    parser.add_argument("--n-types", type=int, default=2)
    parser.add_argument("--cells-per-unit", type=int, default=20)
    args = parser.parse_args(argv)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    if args.out.exists():
        raise SystemExit(f"Refusing to overwrite: {args.out}")
    adata = build_toy(
        seed=args.seed,
        n_donors=args.n_donors,
        n_datasets=args.n_datasets,
        n_types=args.n_types,
        cells_per_unit=args.cells_per_unit,
    )
    adata.write_h5ad(args.out)
    print(
        f"Wrote toy AnnData: {args.out} "
        f"({adata.n_obs} cells, {args.n_donors} donors, {args.n_datasets} datasets)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
