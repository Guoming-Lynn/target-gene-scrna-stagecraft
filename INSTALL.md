# Open-source installation

Part 6 verdict audits require JSON booleans for `blocker`, `smoke_pass`,
`sham_pass`, `run_finished`, and `ko_primary_bh_pass`, plus nonnegative JSON
integers for `ko_eligible_cells`, `n_sign_tests_run`, `family_size`, and
`n_eligible_donors_ko`. Missing or malformed evidence yields `INCONCLUSIVE`.
Successful verdicts additionally require passed smoke/sham gates, a finished
run without blockers, positive KO cell/test counts, at least five eligible
KO donors, and a test count no larger than the frozen family.

Release packaging uses the exact file whitelist in `release-files.txt`.
Review and update it when adding source files; files not listed are excluded.

Use Python 3.10 or 3.11 with `python -m pip install -r requirements.txt`.
CI runs the Python checks on both versions and on Linux, macOS, and Windows.
For the local test and lint tools, also install `python -m pip install -r requirements-dev.txt`.
Run `python scripts/check_environment.py --stage part5` before Part 5; use
`part1` through `part6` for the intended stage. Exit 2 blocks that stage.
The default profile is helper smoke only. Dependency imports do not certify
model compatibility, scientific validity or human annotation quality.

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

Run `python quickstart.py --out quickstart_output` in a new output directory.
Reports persist there; nonzero exit or TIMEOUT is failure, never a skipped pass.
Run `python -m unittest discover -s tests -v` for Python regressions.
Run `Rscript --vanilla tests/audit_boundaries.R` and
`Rscript --vanilla tests/scientific_regression.R` for R regressions.

R is optional until Part 5. Install R from CRAN, then install Bioconductor packages:
`install.packages(c("BiocManager", "yaml", "jsonlite", "digest", "statmod"))`, then
`BiocManager::install(c("limma", "edgeR", "fgsea"))`.
Set `STAGECRAFT_RSCRIPT` when Rscript is not on PATH.

Part 6 uses a separate Python environment and a user-provided, licensed Geneformer source/checkpoint. Set `STAGECRAFT_GENEFORMER_SOURCE` and `STAGECRAFT_GENEFORMER_MODEL`; do not redistribute model weights with this skill.
