
## 2.0.0
- Open-source release; scientific calibration remains incomplete.
- Consolidates the 1.9.11 data-construction tests, cross-platform checks,
  and explicit scientific evidence limits.
- Repair the Part 6 smoke fixture to include required cell identities and
  run the documented quickstart in CI, including both helper smoke suites.

## 1.9.11
- Added focused hash/protocol, review worksheet and rendered plot regressions.
- Hash scanning recognizes bare filenames under inputs/paths and fails when
  no paths are found; protocol section checks now require actual headings.
- Preserve clusters without positive markers in annotation worksheets and
  count dataset-scoped donors in review QC and Part 2 group summaries.
- Fix numeric cluster IDs becoming zero-filled membership heatmaps; require
  both expression layers for Part 3 marker evidence instead of silent fallbacks.
- Plot-only h5ad reads use AnnData; Scanpy marker ranking remains unchanged.
- Statistical engines and scientific evidence ceilings are unchanged.

## 1.9.10
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

## 1.9.9
- Added a short stage route and a project arm inventory retaining missing and
  failed results; within-arm verdicts do not claim project multiplicity control.
- Explicitly separate whole-pipeline calibration and perturbation biological
  validity from protocol execution and embedding consistency.
- Part 3 now rejects mismatched parent barcodes, empty-child deletion and
  missing reviewer evidence; child checkpoints clear stale HVG/PCA/raw state.
- Part 4 now rejects duplicate donor-subtype rows and invalid exposures,
  parses eligibility explicitly and retains zero-eligible subtype slots.
- Added targeted stage handoff and project inventory regressions.

## 1.9.8
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

## 1.9.7
- Added stage dependency gates, donor/block audit boundaries and exploratory
  embedding claim ceilings. Helper smoke tests are not scientific calibration.

## 1.9.6
- Fixed Part 6 paired/unpaired CLI to use a mutually exclusive group.
- Unmapped source-block datasets now fail fast before LODO.
- Clarified that the supported scope is Parts 1–6.
