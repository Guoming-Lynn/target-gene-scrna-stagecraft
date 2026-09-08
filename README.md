# target-gene-scrna-stagecraft

Start with the [stage route and project arm checklist](references/start-here.md).

**Scientific calibration is incomplete.** Passing smoke tests or a `FROZEN_PASS`
verdict does not establish empirical FDR control, adequate power or independent
replication. See [calibration, precision and local freeze evidence](references/calibration-and-provenance.md)
for the executable null harness and its limits.

The 1.9.8 null pilot found elevated BH false discoveries in the repeated-donor
mode. `joint_common_slope` is consequently exploratory only and cannot pass
formal gates; see the linked pilot table before using that mode.

## Two ways to use this package

**As an AI skill:** mount this directory under `~/.claude/skills/`,
`~/.cursor/skills/`, or the Codex skills directory, then provide the agent a
named target feature, tissue context, and library manifest.

**As a CLI workflow:** install the pinned Python dependencies, generate the
toy fixture with `python scripts/generate_toy_data.py --out toy.h5ad`, then
execute the Part 1–6 commands in `SKILL.md` and the corresponding references.
For a quick environment and smoke-gate demonstration, run `python quickstart.py --out quickstart_output`.

An agent skill for **one pre-specified gene** in a multi-study single-cell atlas.

You already know the gene. This skill stops the usual shortcuts: concatenating
libraries before QC, treating Scrublet failure as "all singlets", clustering on
union-gene zeros, naming Leiden clusters after the gene you care about, and
calling 100k cells a sample size.

## Status

**Part 1 is written:** per-library QC → strict common-gene merge → Leiden
resolution grid (human picks) → top-20 markers + per-cluster QC → complete
naming → locked h5ad → publication figures.

**Part 2 is written:** on that locked atlas, a `TARGET_GENE` figure catalog
(detection vs intensity vs donor-unit endpoints). No DEG. See
`references/part2-figures.md`.

**Part 3 is written:** extract a Part 1 lineage, rebuild the manifold from
counts, DELETE low-quality clusters in rounds with a **new HVG after every
deletion**, and name subtypes only on the last clean round. See
`references/part3-compartment-recluster.md` and
`references/part3-figures.md`.

**Part 4 is written:** the Part 2 `TARGET_GENE` catalog on locked subtypes,
plus an identifiability forecast for a later donor-unit slope. Still no DEG.
See `references/part4-subtype-survey.md` and `references/part4-figures.md`.

**Part 5 is written:** donor-unit association on those locked subtypes.
Target-excluded pseudobulk, limma-voom + edgeR, CAMERA/fgsea, donor LOO
**and source-block LODO in the same protocol**. `NOT_ESTIMABLE` is a result.
See `references/part5-donor-association.md` and `references/part5-figures.md`.

**Part 6 is written:** pinned official Geneformer virtual knockout of
`TARGET_GENE` on a locked subtype. Embedding-axis shift, donor-equal sign
tests, KO/OE unpaired when cell sets differ. Not predicted expression.
See `references/part6-virtual-knockout.md` and
`references/part6-figures.md`.

## Install

Build a versioned release archive with `python scripts/package_skill.py`.
The archive is written next to this folder; verify its adjacent `.sha256.txt`
sidecar before distributing it.
Use `python scripts/package_skill.py --out-dir release` when the repository
parent is not writable.

Unzip it so `SKILL.md` is at `target-gene-scrna-stagecraft/SKILL.md`, then copy that folder to:

| Agent | Path |
|---|---|
| Codex | `~/.agents/skills/target-gene-scrna-stagecraft/` or `<repo>/.agents/skills/target-gene-scrna-stagecraft/` |
| Cursor | `~/.cursor/skills/target-gene-scrna-stagecraft/` or `<repo>/.cursor/skills/` |
| Claude Code | `~/.claude/skills/target-gene-scrna-stagecraft/` or `<repo>/.claude/skills/` |

Folder name may be `strict-scrna-stagecraft`; the frontmatter `name` is
`target-gene-scrna-stagecraft`.

Ask the agent to open Part 1 with a named `TARGET_GENE` and a GEO/library list.
After labels are locked, open Part 2 for the gene figure catalog
(`references/part2-figures.md`). To recluster a lineage, open Part 3
(`references/part3-compartment-recluster.md`). After subtypes are named, open
Part 4 (`references/part4-subtype-survey.md`). After identifiability, open
Part 5 (`references/part5-donor-association.md`). For a virtual knockout of
the same gene, open Part 6 (`references/part6-virtual-knockout.md`). For the
full part map, open `references/pipeline-map.md`.

## License

MIT.

## Maintainer

GuomingLin — [GitHub](https://github.com/Guoming-Lynn) ·
[guoming.lin.med@gmail.com](mailto:guoming.lin.med@gmail.com)

## Dependencies and validation

The package includes [the scientific validity contract](references/research-validity.md)
and the explicitly authorized cell-level fallback.
Freeze `estimand.mode` and choose one subtype for `subtype_specific`; use
`joint_common_slope` explicitly for repeated-subtype common-slope analyses.
The runner emits actual-design diagnostics and marks unblocked repeated-donor
edgeR support unavailable. LODO remains internal sensitivity, and Part 6 sign
p values are marked nominal under shared-axis dependence. QC/label audits and
full scientific simulation calibration are stage responsibilities described in
the contract; they are not claimed complete by helper smoke tests.

The cell-level fallback is `scripts/part5_cell_exploratory.py`. It refuses to
run without `I_ACCEPT_CELL_LEVEL_FALSE_POSITIVE_RISK`, records that acceptance
in JSON, and cannot enter the formal verdict. Use
`scripts/simulation_contract.py` to create a frozen calibration manifest for
null, confounding, collinearity, attrition and failed-holdout scenarios.

Formal figure guidance is centralized in `references/figure-contract.md`,
`references/figure-statistics-contract.md`, and
`references/visual-qa-contract.md`. Use
`scripts/validate_figure_manifest.py` before treating a figure folder as
publication-ready.

New executors should read `references/glossary.md` first. The structural
validators intentionally do not judge biological correctness; they now reject
very short/empty protocols, but a passing validator is not a substitute for
the required evidence worksheet and human review.

Use a dedicated Python 3.10 or 3.11 environment. CI exercises both versions on
Linux, macOS, and Windows. Dependency ranges are installation constraints, not
a guarantee for every combination.
For interpreter selection and the separate Part 6 environment see INSTALL.md:

The validated project baseline is recorded in `environment.example.yaml`.
Copy it to a stage manifest and replace every placeholder with exact values;
formal runs must archive Python package resolution, R `sessionInfo()`, and the
exact model/checkpoint hashes without including biological identifiers.

```bash
python -m pip install -r requirements.txt
python -m pip install -r requirements-dev.txt
python scripts/_part5_smoke.py
python scripts/_part6_smoke.py
python -m unittest discover -s tests -p "test_*.py"
python -m pytest -q tests
Rscript tests/scientific_regression.R
```

Part 5 needs R with `Matrix`, `limma`, `edgeR`, `statmod`, `jsonlite`,
`yaml`, `digest`, and `fgsea` (Bioconductor for limma/edgeR/fgsea).
Install packages in the project's R environment and record `sessionInfo()`.
Part 6 additionally needs the official Geneformer revision, matching PyTorch,
Transformers, model weights and dictionaries frozen in its protocol. The
helper smoke tests do not execute Geneformer or validate a GPU/checkpoint.

Run the legacy R syntax check with (the actual synthetic integration check is
`Rscript tests/scientific_regression.R` above):

```bash
python tests/run_r_integration.py --rscript /path/to/Rscript
```

The current release physically excludes the target from pseudobulk outcomes, aligns
source metadata by unit ID, computes donor LOO and source-block LODO, derives
gene evidence/audit flags from model outputs, and rejects incomplete smoke
metrics and mismatched GMT hashes. See the Part 5 specification for the
remaining pathway robustness and project-specific diagnostic handoff.
