# Changelog

All notable changes to this project are documented here.

## [2.3.1] - 2026-09-23

### Fixed
- Stage reports for `NOT_ESTIMABLE`, `INCONCLUSIVE`, `REPRODUCTION_FAILED`, `COLLINEAR_UNINTERPRETABLE`, and non-passing Part 6 tokens no longer list a directional association as allowed.
- Claim review does not let a negation cross `but` / `但是`, does not treat `非常` as a negation, and flags `driver`, `abolishes`, `mediates`, `required for the phenotype`, and `replicated across datasets`.
- Population labels are flagged only with `claim_lint.py --target GENE`. Lineage labels such as `CD45+ cells` are left alone.
- A malformed `model_audit.json` or evidence table exits 1 and writes no report. Provenance hashes each input once and accepts a table directory outside the stage.
- `stage_status.py` accepts `FROZEN_PROTOCOL.md`, prints the report and claim-lint commands on separate lines, and includes the Part 6 config, control, and smoke steps. Part 6 coverage includes `ko_primary_bh_pass`.
- No statistical output or threshold changed.

## [2.3.0] - 2026-09-23

### Added
- `scripts/claim_lint.py` flags draft sentences that cross the claim boundaries. Findings are review flags unless `--strict` is set.
- `scripts/stage_report.py` renders `06_reports/PARTN_REPORT.md` from `verdict.json` and `model_audit.json`.
- `scripts/stage_status.py` lists completed Part 5 or Part 6 helper steps and prints the next command.
- No statistical output or threshold changed.

## [2.2.1] - 2026-09-23

### Fixed
- Part 5 reads donor, cell, and gene identifiers as text in R, and Part 4 figure tables use the same identity reader.
- Multi-file helpers write partials and rename them only after every file succeeds. `require_new` reserves the path exclusively.
- `part5_source_blocks.py` rejects an incomplete source map before writing. `cluster_review_tables.py` refuses to overwrite review sheets.
- `part6_axes.py` requires a `.npz` path, records skips, and stops when skipped donors exceed `--max-skipped-fraction` (default 0.5).
- Scientific stops for an empty eligible set, a frozen endpoint hash, too few donor axes, and an enlarged sign-test family exit 2.

### Changed
- Script imports go through `ensure_repo_on_path`. Helper commands are listed in `references/cli-catalog.md`.
- Linux Python 3.11 CI installs the built wheel in a clean virtualenv and uploads the quickstart output.

## [2.2.0] - 2026-09-23

### Fixed
- Blocked Part 5 fits re-estimate `duplicateCorrelation` after the second voom and pass that consensus into CAMERA.
- Part 5 pathways require fgsea and read `q_cut`, CAMERA `inter.gene.cor`, and fgsea Monte Carlo settings from config.
- `model_audit.json` records R and package versions. `part5_verdict.py` reports a local freeze receipt when one is present.
- Identity columns stay strings, boolean gates reject non `true`/`false`/`1`/`0` values, and multi-file outputs refuse every sidecar before writing.
- Skipped Part 6 donor axes are written beside the axis archive.
- Python and the Part 5 R engine share one unlikely-arm pattern file.

### Changed
- Dev checks pin ruff 0.16.8.
- Dependabot skips `references/part6-requirements.txt`. Those pins are a historical Geneformer candidate, not the helper runtime.
- The skill archive uses `skill-files.txt` and omits `.github/`.
- CLI exit codes 2 and 3 use the shared `EXIT_GATE` and `EXIT_NOT_ESTIMABLE` constants.
- Linux Python 3.11 CI reports `stagecraft` coverage and builds a wheel. The R job runs a pathway contract.

## [2.1.1] - 2026-09-21

### Fixed
- Windows CI: Part 6 CLI `--help` no longer prints characters that a cp1252
  console cannot encode.
- Status prints from Part 5/6 helpers use ASCII hyphens instead of em dashes.

### Changed
- GitHub Actions use `actions/checkout@v7` and `actions/setup-python@v7`.
- pytest is allowed through 9.x (`>=8,<10`).
- Dependabot ignores Geneformer-pinned packages and major clustering extras.

## [2.1.0] - 2026-09-20

### Added
- Installable `stagecraft` helper package, `pyproject.toml`, and `CITATION.cff`.
- `scripts/part1_init.py` stage scaffold and `scripts/demo_toy_run.py` helper path that consumes toy data.
- `requirements-part1.txt` for clustering extras not imported by helper CLIs.

### Changed
- README, Chinese README, INSTALL, and SKILL now state which parts have CLIs versus specification-only execution.
- Quickstart runs the toy helper demo instead of writing an unused h5ad.
- `simulation_contract.py` labels itself `MANIFEST_ONLY`. `--declared-replicates` (alias `--replicates`) records a count and does not fit models.
- Overwrite refusal, SHA-256, YAML loading, and Excel CSV encoding share one helper module.
- Part 6 status is `specified / not turnkey` in the pipeline map.

### Fixed
- Pathway runner now defaults `set.seed` to 42 when the config omits `random_seed`.
- `check_stage_layout.py` uses argparse.
- Install and test recipes now match CI (`pytest`).
- Chinese protocol template routes Parts 1, 3, and 4.

## [2.0.0]
- Open-source release; scientific calibration remains incomplete.
- Consolidates the 1.9.11 data-construction tests, cross-platform checks,
  and explicit scientific evidence limits.
- Repair the Part 6 smoke fixture to include required cell identities and
  run the documented quickstart in CI, including both helper smoke suites.

## [1.9.11]
- Added focused hash/protocol, review worksheet and rendered plot regressions.
- Hash scanning recognizes bare filenames under inputs/paths and fails when
  no paths are found; protocol section checks now require actual headings.
- Preserve clusters without positive markers in annotation worksheets and
  count dataset-scoped donors in review QC and Part 2 group summaries.
- Fix numeric cluster IDs becoming zero-filled membership heatmaps; require
  both expression layers for Part 3 marker evidence instead of silent fallbacks.
- Plot-only h5ad reads use AnnData; Scanpy marker ranking remains unchanged.
- Statistical engines and scientific evidence ceilings are unchanged.

## [1.9.10]
- Added hand-calculated data construction tests for pseudobulk counts and
  exposure sidecars, source mapping, LODO geometry, endpoint coverage,
  token classification and cell-to-donor eligibility, including actual file I/O.
- Reject duplicate barcodes, missing identities and colliding/split donor keys;
  pooled mixed labels now report POOLED instead of the first subtype.
- Parse source eligibility and token-presence booleans explicitly; reject
  ambiguous LOO mappings, duplicate endpoint members and invalid token counts/ranks.
- Fix the LODO check to permit coincident directions while comparing against
  the normalized training-only mean. Score CSVs now require matching cell_id.
- Part 6 effect tables now require cell_id and reject duplicate effects within
  each donor/target/perturbation/population/endpoint; empty eligible output retains headers.
- These are numerical regression tests, not whole-pipeline scientific calibration.

## [1.9.9]
- Added a short stage route and a project arm inventory retaining missing and
  failed results; within-arm verdicts do not claim project multiplicity control.
- Explicitly separate whole-pipeline calibration and perturbation biological
  validity from protocol execution and embedding consistency.
- Part 3 now rejects mismatched parent barcodes, empty-child deletion and
  missing reviewer evidence; child checkpoints clear stale HVG/PCA/raw state.
- Part 4 now rejects duplicate donor-subtype rows and invalid exposures,
  parses eligibility explicitly and retains zero-eligible subtype slots.
- Added targeted stage handoff and project inventory regressions.

## [1.9.8]
- Restricted joint_common_slope to exploratory evidence after a 100-draw null
  pilot per mode showed a repeated-donor BH false-discovery warning signal.
- Added an executable pseudobulk global-null harness using the real Part 5
  gene engine, donor LOO and source LODO; full scientific calibration remains incomplete.
- Added three-platform R CI with actual model tests, including repeated donors.
- Removed contradictory external-replication examples and reject legacy tokens.
- Added nominal donor-median intervals and forest rendering, with unbounded
  small-sample intervals and shared-axis dependence explicitly labeled.
- Report calibration, precision, power and chronology limitations in outputs.
- Added local protocol freeze receipts and documented their evidence limits.
- Fixed Part 6 figure CLI argument parsing; pinned correctness lint configuration.
- Normalized Markdown line endings and added .gitattributes.

## [1.9.7]
- Added stage dependency gates, donor/block audit boundaries and exploratory
  embedding claim ceilings. Helper smoke tests are not scientific calibration.

## [1.9.6]
- Fixed Part 6 paired/unpaired CLI to use a mutually exclusive group.
- Unmapped source-block datasets now fail fast before LODO.
- Clarified that the supported scope is Parts 1–6.
