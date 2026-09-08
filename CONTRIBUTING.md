# Contributing

Thank you for improving the workflow.

## Before opening a pull request

- Keep changes scoped to the declared stage or contract.
- Do not add biological identifiers, project data, model weights, or absolute local paths.
- Preserve the distinction between formal protocol gates and empirical scientific calibration.
- Add a focused regression test when changing a data-construction, gate, or verdict behavior.

Run the checks appropriate to the change:

```bash
ruff check .
python -m pytest -q tests
python quickstart.py --out quickstart_output
Rscript --vanilla tests/audit_boundaries.R
Rscript --vanilla tests/scientific_regression.R
```

Part 6 helper tests do not run Geneformer or distribute any checkpoint. Keep model-specific work in a separately pinned, user-provided environment.

## Pull requests

Release archives include only the exact paths in `release-files.txt`.
When adding source files, review their contents and update that manifest.
Never list local data, credentials, analysis outputs, or model weights.

Describe the concrete behavior before and after the change, the applicable scientific boundary, and validation performed. Do not use internal LODO results as external replication.
