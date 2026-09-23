# Open-source installation

Use Python 3.10 or 3.11.

```bash
python -m pip install -r requirements.txt
python -m pip install -e .
```

`pip install -e .` installs the shared `stagecraft` package. Helper commands stay `python scripts/...`; this install does not put them on `PATH`. Python is capped below 3.12 because that is the validated range for the pinned Scanpy and AnnData stack.

CI runs those Python versions on Linux, macOS, and Windows.

- Runtime: `requirements.txt` (numpy, pandas, scipy, PyYAML, anndata, scanpy, matplotlib, seaborn).
- Tests and lint: `python -m pip install -r requirements-dev.txt`.
- Part 1/3 clustering extras (Scrublet, Harmony, igraph, Leiden): `python -m pip install -r requirements-part1.txt`.

Run `python scripts/check_environment.py --stage part5` before Part 5; use
`part1` through `part6` for the intended stage. Exit 2 blocks that stage.
The default profile is helper smoke only. Dependency imports do not certify
model compatibility, scientific validity or human annotation quality.
`--stage part1` and `--stage part3` require the clustering extras.

If `python` is missing, select an installed interpreter explicitly: `python3`
on Unix, `py -3.11` on Windows, or `conda run -n YOUR_ENV python`.
For PowerShell an executable path requires the call operator: `& '/path/to/python.exe'`.
Conda environments should be activated or launched with `conda run` so their
native DLL directories are configured. No tests ran when the interpreter cannot start.
On Windows, an unactivated Conda environment can import matplotlib successfully
and still crash while rendering (`0xc06d007f`). Activate the environment or use
`conda run -n YOUR_ENV python ...` for the figure smoke check; an import-only
environment PASS does not test native drawing. Do not hard-code another user's
DLL directory in the skill.
All Python subprocesses reuse `sys.executable`; STAGECRAFT_PYTHON is not supported.

## Checks (same as CI)

```bash
python -m pytest -q tests
python quickstart.py --out quickstart_output
Rscript --vanilla tests/audit_boundaries.R
Rscript --vanilla tests/scientific_regression.R
```

Nonzero exit or TIMEOUT is failure, never a skipped pass.
`tests/run_r_integration.py --rscript /path/to/Rscript` wraps the two R files.

## R (required for Part 5)

Install R from CRAN, then:

```r
install.packages(c("BiocManager", "Matrix", "yaml", "jsonlite", "digest", "statmod"))
BiocManager::install(c("limma", "edgeR", "fgsea"))
```

`fgsea` is required by `scripts/check_environment.py --stage part5`.
Set `STAGECRAFT_RSCRIPT` when Rscript is not on PATH.
See [references/r-requirements.md](references/r-requirements.md).

## CLI exit codes

| Code | Meaning |
|---|---|
| 0 | Success |
| 1 | Usage or validation error |
| 2 | Environment/gate failure |
| 3 | Not estimable / no KO-eligible cells |

## Part 6

Part 6 uses a separate Python environment and a user-provided, licensed Geneformer source/checkpoint. Set `STAGECRAFT_GENEFORMER_SOURCE` and `STAGECRAFT_GENEFORMER_MODEL`; do not redistribute model weights with this skill.

Verdict audits require JSON booleans for `blocker`, `smoke_pass`,
`sham_pass`, `run_finished`, and `ko_primary_bh_pass`, plus nonnegative JSON
integers for `ko_eligible_cells`, `n_sign_tests_run`, `family_size`, and
`n_eligible_donors_ko`. Missing or malformed evidence yields `INCONCLUSIVE`.
Successful verdicts additionally require passed smoke/sham gates, a finished
run without blockers, positive KO cell/test counts, at least five eligible
KO donors, and a test count no larger than the frozen family.

## Releases

Release packaging uses the exact file whitelist in `release-files.txt`.
Review and update it when adding source files; files not listed are excluded.
