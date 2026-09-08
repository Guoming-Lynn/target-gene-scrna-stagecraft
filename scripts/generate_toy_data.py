"""Generate a small, biologically neutral AnnData fixture for smoke tests."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--out', required=True); ap.add_argument('--seed', type=int, default=7)
    a = ap.parse_args(); rng = np.random.default_rng(a.seed)
    X = rng.poisson(1.0, (200, 50)).astype('int32'); X[:, 0] += rng.poisson(2, 200)
    obs = pd.DataFrame({'donor_id': [f'd{i%10:02d}' for i in range(200)], 'dataset': [f's{i%2}' for i in range(200)], 'cell_type': [f't{i%4}' for i in range(200)]})
    var = pd.DataFrame(index=['TARGET_FEATURE'] + [f'G{i:03d}' for i in range(49)])
    adata = ad.AnnData(X=X, obs=obs, var=var); adata.layers['counts'] = X.copy()
    Path(a.out).parent.mkdir(parents=True, exist_ok=True); adata.write_h5ad(a.out)
    print(f'Wrote toy AnnData: {a.out} (200 cells, 10 donors)')
if __name__ == '__main__': main()

