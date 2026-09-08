"""Stage-specific dependency gate; not a scientific or model parity certificate."""
import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

PROFILES = {
    'smoke': ['numpy', 'pandas', 'scipy', 'anndata', 'yaml'],
    'part1': ['numpy', 'pandas', 'anndata', 'scanpy', 'scrublet', 'harmonypy', 'igraph', 'leidenalg'],
    'part2': ['numpy', 'pandas', 'anndata', 'matplotlib', 'seaborn', 'yaml'],
    'part3': ['numpy', 'pandas', 'anndata', 'scanpy', 'harmonypy', 'igraph', 'leidenalg'],
    'part4': ['numpy', 'pandas', 'anndata', 'matplotlib', 'seaborn', 'yaml'],
    'part5': ['numpy', 'pandas', 'scipy', 'anndata', 'yaml'],
    'part6': ['numpy', 'pandas', 'scipy', 'anndata', 'yaml', 'torch', 'transformers', 'datasets', 'geneformer'],
}

def check(command, timeout, env=None):
    try:
        result = subprocess.run(command, capture_output=True, text=True, errors='replace', timeout=timeout, env=env)
        return dict(status='PASS' if result.returncode == 0 else 'FAILED', exit_code=result.returncode,
                    detail=(result.stdout + result.stderr).strip())
    except (OSError, subprocess.TimeoutExpired) as exc:
        return dict(status='FAILED', detail=str(exc))

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', choices=list(PROFILES), default='smoke')
    parser.add_argument('--rscript', default=os.environ.get('STAGECRAFT_RSCRIPT'))
    parser.add_argument('--timeout', type=float, default=60)
    args = parser.parse_args(argv)
    results = {}
    env = os.environ.copy()
    source = env.get('STAGECRAFT_GENEFORMER_SOURCE')
    if source:
        env['PYTHONPATH'] = source + os.pathsep + env.get('PYTHONPATH', '')
    for module in PROFILES[args.stage]:
        print('Checking ' + module, file=sys.stderr, flush=True)
        code = f"import importlib; m=importlib.import_module({module!r}); print(getattr(m, '__version__', 'import passed'))"
        results[module] = check([sys.executable, '-c', code], args.timeout, env)
    if args.stage == 'part5':
        r = args.rscript or shutil.which('Rscript')
        expr = 'p <- c("Matrix","limma","edgeR","fgsea","jsonlite","yaml","digest","statmod"); for (x in p) { if (!requireNamespace(x, quietly=TRUE)) stop(paste("Missing package",x)); cat(x,as.character(packageVersion(x)),"\\n") }'
        results['R_packages'] = check([r, '--vanilla', '-e', expr], args.timeout) if r else dict(status='FAILED', detail='Set STAGECRAFT_RSCRIPT or add Rscript to PATH')
    if args.stage == 'part6':
        model = Path(env.get('STAGECRAFT_GENEFORMER_MODEL', '__unconfigured__'))
        valid = model.is_dir() and (model / 'config.json').is_file() and any((model / x).is_file() for x in ['model.safetensors', 'pytorch_model.bin', 'model.safetensors.index.json', 'pytorch_model.bin.index.json'])
        results['model_files'] = dict(status='PASS' if valid else 'FAILED', detail='Model load, revision/hash and official parity gates remain required before inference')
    passed = all(row['status'] == 'PASS' for row in results.values())
    print(json.dumps(dict(stage=args.stage, status='PASS' if passed else 'FAILED', results=results), indent=2), flush=True)
    return 0 if passed else 2

if __name__ == '__main__':
    raise SystemExit(main())
