# Upstream implementation influences

This project uses the following public skill packs as design references. No
upstream source code is copied into the pipeline; the local protocol remains
the authority for estimands, thresholds, exclusions, and claims.

Pinned reference: `anthropics/life-sciences` commit
[`e96556b637b56d6cc3a5ad33987009be9e60aa5c`](https://github.com/anthropics/life-sciences/tree/e96556b637b56d6cc3a5ad33987009be9e60aa5c).
The three referenced directories are distributed under Apache License 2.0.

## `single-cell-rna-qc`

Borrowed patterns:

- separate metric calculation, outlier detection, hard thresholds, filtering,
  and summary reporting into independently testable steps;
- write before/after QC summaries and threshold visualizations;
- make species-specific gene-pattern assumptions explicit rather than hidden
  in code.

Adaptation in this project: QC is still performed one library at a time;
Scrublet runs on the loose-filter population, MAD thresholds are estimated on
singlet candidates, failed calls exclude a library from the primary atlas, and
strict-common genes are frozen before clustering. The upstream skill's generic
single-pass filter must not be substituted for these rules.

## `scvi-tools`

Borrowed patterns:

- require an immutable integer-count layer before modeling;
- validate AnnData compatibility and batch annotations before a run;
- save model/config provenance and evaluate the learned representation after
  training.

Adaptation in this project: scVI/scANVI is an optional exploratory or
integration branch. It cannot change the frozen Part 1 atlas, donor-level
estimand, source-block LODO, or Part 6 claim ceiling. Any latent representation
must carry its model version, seed, input hash, and batch key.

## `nextflow-development`

Borrowed patterns:

- environment preflight before execution;
- a small test profile before real data;
- explicit samplesheet validation, pinned pipeline versions, resumable runs,
  and containerized execution for FASTQ-scale work.

Adaptation in this project: these controls are recommended for a future raw
FASTQ/Part 1 executor. They do not turn the current prose-only Part 1 into an
implemented tool, and `-resume` must never overwrite a frozen artifact.

## Citation and change control

When an implementation adopts one of these patterns, cite this file and the
pinned upstream commit in the stage report. Record any divergence from the
upstream defaults in the stage config and protocol hash.
