# target-gene-scrna-stagecraft

[中文说明](README.zh-CN.md) · English · v2.3.1

[![CI](https://github.com/Guoming-Lynn/target-gene-scrna-stagecraft/actions/workflows/ci.yml/badge.svg)](https://github.com/Guoming-Lynn/target-gene-scrna-stagecraft/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.10 | 3.11](https://img.shields.io/badge/python-3.10%20%7C%203.11-blue.svg)](INSTALL.md)

## What this answers

For one pre-specified gene, can its signal be measured and associated across
multi-study single-cell data at the donor-unit level, with QC, source
sensitivity, and an explicit ceiling on biological claims?

This is a protocol-driven analysis skill, not a generic Scanpy tutorial: it
keeps libraries separate during QC, treats `dataset × donor_id` as the unit,
and forbids cell-level DEG shortcuts.

## What you can run today

| Surface | What exists |
|---|---|
| **AI skill** | Full Part 1–6 specification in `SKILL.md` and `references/` |
| **CLI helpers** | Review tables, figure catalogs, Part 5/6 construction and verdicts |
| **Reporting helpers** | Stage progress, rendered reports, and claim-boundary review |
| **Not bundled** | A turnkey Part 1 Scanpy runner (QC, Harmony, Leiden). Use `scripts/part1_init.py` plus the Part 1 spec |
| **Part 6** | Specified and audited; not turnkey. Needs a separate licensed Geneformer environment |

Start with the [stage route](references/start-here.md) and the [reference index](references/README.md).

**Scientific calibration is incomplete.** Passing smoke tests or a `FROZEN_PASS`
verdict does not establish empirical FDR control, adequate power or independent
replication. See [calibration, precision and local freeze evidence](references/calibration-and-provenance.md).
The 1.9.8 null pilot found elevated BH false discoveries in the repeated-donor
mode. `joint_common_slope` is exploratory only.

## Try it first

Python 3.10 or 3.11:

```bash
python -m pip install -r requirements.txt
python quickstart.py --out quickstart_output
```

Quickstart runs environment checks, a toy helper demo (pseudobulk → source
blocks → eligibility → manifest-only calibration file), and Part 5/6 helper
smokes. It does not claim a biological result or execute Geneformer.

Optional clustering extras for Part 1/3 Scanpy execution:

```bash
python -m pip install -r requirements-part1.txt
```

## Install

Runtime Python packages: `requirements.txt`.
Developer checks: `requirements-dev.txt` (`pytest`, `ruff`).
Part 1/3 clustering libraries: `requirements-part1.txt`.

```bash
python -m pip install -r requirements.txt
python -m pip install -e .
```

Editable install provides the `stagecraft` helpers. Run CLIs as `python scripts/...`. Python stays on 3.10 or 3.11.

Part 5 additionally needs R with `Matrix`, `limma`, `edgeR`, `fgsea`,
`statmod`, `jsonlite`, `yaml`, and `digest`. See [INSTALL.md](INSTALL.md)
and [R packages](references/r-requirements.md).

Interpreter selection, Windows Conda, and Part 6 model paths: [INSTALL.md](INSTALL.md).
The maintainer baseline is [validated-environment.md](references/validated-environment.md);
CI also runs Linux and macOS. Copy `environment.example.yaml` into a stage
manifest and replace placeholders before a formal run.

## Use it as a skill or CLI

**As an AI skill:** mount this directory under `~/.claude/skills/`,
`~/.cursor/skills/`, or the Codex skills directory, then provide the agent a
named target feature, tissue context, and library manifest.

Unzip a release so `SKILL.md` is at `target-gene-scrna-stagecraft/SKILL.md`:

| Agent | Path |
|---|---|
| Codex | `~/.agents/skills/target-gene-scrna-stagecraft/` or `<repo>/.agents/skills/target-gene-scrna-stagecraft/` |
| Cursor | `~/.cursor/skills/target-gene-scrna-stagecraft/` or `<repo>/.cursor/skills/` |
| Claude Code | `~/.claude/skills/target-gene-scrna-stagecraft/` or `<repo>/.claude/skills/` |

Use the repository name `target-gene-scrna-stagecraft`. Do not rename the folder.

**As a CLI workflow:** init a Part 1 stage, execute the Part 1 spec until labels
are locked, then run the helper commands in `SKILL.md`.

```bash
python scripts/part1_init.py --out analysis/01_qc_global_atlas --gene TARGET_GENE
python scripts/demo_toy_run.py --out toy_demo_output
```

## Six parts

| Part | Job | Status | Code in this repo |
|---|---|---|---|
| **1** | Per-library QC, strict-common merge, Leiden, human labels | written spec | `part1_init.py`, `cluster_review_tables.py` |
| **2** | Global `TARGET_GENE` survey (no DEG) | written | `part2_figures.py` |
| **3** | Compartment recluster with new HVG after each DELETE | written | `part3_*` helpers |
| **4** | Subtype survey + identifiability forecast | written | `part4_figures.py`, `part4_identifiability.py` |
| **5** | Donor-unit association, LOO, source-block LODO | written | `part5_*` Python + R |
| **6** | Virtual KO embedding shift | specified / not turnkey | `part6_*` helpers; Geneformer is user-provided |

Ask the agent to open Part 1 with a named `TARGET_GENE` and a GEO/library list.
After labels are locked, open Part 2 (`references/part2-figures.md`). To
recluster a lineage, open Part 3. After subtypes are named, open Part 4, then
Part 5. Part 6 is a separate Geneformer branch. Full map:
[pipeline-map.md](references/pipeline-map.md).

Input layers and objects: [objects-and-layers.md](references/objects-and-layers.md).
Stage directory contract: [stage-layout.md](references/stage-layout.md).
Worked anti-patterns: [examples.md](examples.md).

## Validation (helpers, not scientific certification)

The [scientific validity contract](references/research-validity.md) is
authoritative. Freeze `estimand.mode` before Part 5. `joint_common_slope` cannot
pass formal gates. New executors should read [glossary.md](references/glossary.md).
Structural validators do not judge biological correctness.

```bash
python -m pytest -q tests
python quickstart.py --out quickstart_output
Rscript --vanilla tests/audit_boundaries.R
Rscript --vanilla tests/scientific_regression.R
```

Those commands match CI. `tests/run_r_integration.py` is a wrapper that runs
the two R regression files; it is not a syntax linter.

Figure contracts: [figure-contract.md](references/figure-contract.md),
[figure-statistics-contract.md](references/figure-statistics-contract.md),
[visual-qa-contract.md](references/visual-qa-contract.md).

The cell-level fallback is `scripts/part5_cell_exploratory.py`. It refuses to
run without `I_ACCEPT_CELL_LEVEL_FALSE_POSITIVE_RISK` and cannot enter the
formal verdict. `scripts/simulation_contract.py` writes a **manifest only**;
it does not run 1,000 calibration replicates.

## Release packaging

```bash
python scripts/package_skill.py --out-dir release
```

The skill archive is built from `skill-files.txt` and does not include `.github/`.
`release-files.txt` remains the full repository whitelist. Check the sibling `.sha256.txt` before publishing the zip.

## License

MIT. Cite with [CITATION.cff](CITATION.cff).

## Maintainer

GuomingLin — [GitHub](https://github.com/Guoming-Lynn) ·
[guoming.lin.med@gmail.com](mailto:guoming.lin.med@gmail.com)

[CHANGELOG](CHANGELOG.md) · [CONTRIBUTING](CONTRIBUTING.md) ·
[SECURITY](SECURITY.md) · [CODE OF CONDUCT](CODE_OF_CONDUCT.md)
