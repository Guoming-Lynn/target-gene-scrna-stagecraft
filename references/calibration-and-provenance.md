# Calibration, precision and local freeze evidence

## Version 1.9.8 null pilot

Local Windows R 4.5.2, limma 3.66.0, edgeR 4.8.2; 18 donors, three sources,
1,000 genes; 100 independent draws per mode. Seeds 314159-314258 (single
subtype), 414159-414258 (two repeated subtypes per donor). All 200 runs and
their source holdouts completed. Per-draw logs and source hashes were retained
outside the public package; these results can be rerun with the harness above.

| Mode | Raw type-I error (MCSE) | Global-null BH FDR (exact 95% interval) | Marginal 95% CI coverage |
|---|---|---|---|
| subtype_specific | 0.04984 (0.00078) | 0.06 (0.0223-0.1260) | 0.95016 |
| joint_common_slope | 0.05616 (0.00082) | 0.11 (0.0562-0.1883) | 0.94384 |

This is a warning signal for the repeated-donor mode, not proof of its cause.
Version 1.9.8 therefore prevents that mode from passing formal gates and labels
it `EXPLORATORY_ONLY_CALIBRATION_CONCERN`; coefficient estimates remain available.
The simulation used the same fitting code, prior to adding this final verdict
restriction. No fitted p values or thresholds were changed in response to the
pilot. Unblocked edgeR support was already unavailable in that mode.

There were 100 repetitions per mode, **not 1,000**. The interval around 6% is too
wide to certify the single-subtype model. Further null and non-null scenarios,
independent seeds and model investigation remain necessary before validation.

`FROZEN_PASS` describes execution of frozen gates. It is not a claim that
false-discovery rates were empirically calibrated or that the design can detect
the chosen effect floor. The runner and verdict explicitly report
`scientifically_calibrated: false` and `NOT_ESTABLISHED_FOR_THIS_ANALYSIS`.
A small simulation elsewhere cannot change these fields for a new analysis.

## Executable null experiment

The manifest-only helper remains a planning tool. To actually execute the
pseudobulk gene engine, including all donor LOO and source LODO fits:

```bash
Rscript --vanilla scripts/calibrate_part5_null.R null_single 1000 271828 subtype_specific 18 3 1000
Rscript --vanilla scripts/calibrate_part5_null.R null_repeated 1000 371828 joint_common_slope 18 3 1000
```

Arguments are new output directory, repetitions, seed, mode, independent donor
count, source count and gene count. Donors must divide evenly among sources.
The repeated mode has two subtype rows per donor, shares simulated donor
random intercepts, and invokes duplicateCorrelation in the real runner.
Unblocked edgeR support remains unavailable for that mode by design.

The generator freezes scenario parameters, seeds, source hashes and runtime
versions before fitting. Every replicate retains pseudobulk inputs, config,
model outputs and logs. `replicates.csv` records failures without silently
dropping them. `calibration_report.json` reports raw type-I error, global-null
BH FDR, marginal CI coverage, failed-fit and failed-holdout rates with Monte
Carlo uncertainty. Under a global null, FDR is the probability of *any* BH
rejection, not the fraction of rejected genes. Correlated genes are not counted
as independent simulation repetitions.

This experiment starts at pseudobulk counts. It does not simulate raw-cell QC,
label selection, pathway analysis or Geneformer. Its source shifts do not cover
all source-confounding structures. Non-null power, MDE, depth collinearity,
differential attrition and other sample-size regimes require additional frozen
experiments. Never describe this harness or its pilot run as complete pipeline
calibration. Do not tune thresholds against these draws and report them as an
independent validation; freeze new seeds/scenarios after any tuning.

## Precision instead of post-hoc power

Eligibility emits `power_status: NOT_ASSESSED_REQUIRES_DESIGN_SIMULATION` and
an empty MDE. After fitting, the audit reports the median and 90th percentile
of marginal coefficient CI half-widths, and how often those exceed the frozen
effect floor. These quantities describe uncertainty in the measured coefficient
scale. They are neither achieved power nor a BH-adjusted MDE.

Power/MDE need a prospective effect grid, a defined rejection rule, dispersion,
donor/source structure and non-null fraction. Retrospective power calculated
from the observed coefficient is not added. A logFC cutoff after zero-null BH
testing does not test a minimum-effect null (see research-validity.md).

## Local freeze receipts

Before generating tables or pseudobulk files, place the protocol and frozen
configuration in the stage and run:

```bash
python scripts/check_protocol.py stage/00_protocol_manifest/PROTOCOL.md --stage-root stage --freeze
```

Before each subsequent stage command, rerun without `--freeze`:

```bash
python scripts/check_protocol.py stage/00_protocol_manifest/PROTOCOL.md --stage-root stage
```

The receipt hashes the protocol and files in `00_protocol_manifest`, refuses a
new freeze if `02_tables` or `03_pseudobulk` already contains artifacts, and
rejects altered/missing/added manifest files or outputs older than the receipt.
An amendment goes into a new stage with explicit links to the old results;
do not recreate the receipt to conceal changes.

This is an explicitly invoked local consistency check. Calling the R runner
directly does not invoke it, and the runner reports chronology `NOT_VERIFIED`.
File copies, restored timestamps and clock changes can require manual review.
Neither local timestamps, a locally editable receipt nor Git commit dates prove
that an analyst had not already seen results. Trusted preregistration requires
an independently timestamped record and remains outside this helper.
