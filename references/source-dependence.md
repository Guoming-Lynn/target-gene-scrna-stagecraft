# Source-block dependence

## The failure mode

A multi-GEO atlas can report 20 donors, 3 datasets, LOO 20/20, and q = 1e-18, and still be a **single-study result**.

Typical pattern:

- Two consecutive accessions can share interleaved donor numbering.
- Together they are 17/20 donors.
- The third accession has almost no variance in the exposure, so it cannot identify a slope.
- Donor LOO still looks perfect: each iteration leaves 16 donors from the same study in the model.

Donor LOO measures **leverage of one person**. It does not measure **replaceability of a paper**.

## Procedure (pre-register in the protocol)

### S1. Map sources

Assign `source_block` from repository-internal evidence only (metadata, donor IDs, consecutive GSE). Do not browse the web to "confirm" a paper unless the protocol says so.

Tiers: `SAME_STUDY_LIKELY` | `SAME_STUDY_CONFIRMED` | `UNRESOLVED`. The holdout design does not change with the tier.

Also write exposure mean/sd/IQR per dataset. Flag `NO_EXPOSURE_RANGE` when sd
is below the **modeled exposure** floor frozen in the protocol. Code defaults
(`scripts/part5_analysis_config.example.yaml` and `--sd-floor` / `--min-sd`):

| Scale | Starting floor | Where it applies |
|---|---|---|
| `jeffreys_per_10pct` | 0.5 | Part 5 eligibility and source-block audit |
| log2(CPM+1) | 0.15 | Part 4 forecast and Part 5 alternative exposure |
| detection fraction | 0.02 | Part 4 forecast only |

Do not mix these floors across scales. Amend them in the protocol; do not
silently reuse 0.5 on a 0–1 detection fraction.

### S2. Hold out the whole block

Refit the **parent design** on:

- `FULL`
- drop external block
- drop the suspected sibling block
- drop the largest sibling accession (optional split-half)

If only one `dataset` level remains, drop that term, record a deviation, and still evaluate residual df.

If residual df is below the protocol minimum or the design is rank-deficient: `NOT_ESTIMABLE`. Do not delete covariates until a number appears.

### S3. Within-dataset slopes

Fit inside each accession. Agreement between two halves of one study is split-half consistency, **not** cross-source replication.

### S4. Relabel existing donor LOO

Do not rerun LOO unless the protocol asks. Count what fraction of existing iterations dropped a donor from the dominant block. ≥0.75 → tag `LOO_SOURCE_UNBALANCED`.

## Verdict

Use the canonical ordered table in `scripts/part5_verdict_table.example.yaml`
and the grammar in [evidence-and-verdicts.md](evidence-and-verdicts.md).
Successful source holdouts establish internal sensitivity only. They cannot
emit an external-replication verdict, regardless of sign agreement or q value.
A failed or sign-reversing dominant holdout gives `SINGLE_SOURCE_DEPENDENT`;
incomplete other holdouts prevent `FROZEN_PASS`. A successful holdout alone
does not satisfy the remaining discovery, method and diagnostic gates.

Reproduction of the parent `FULL` coefficient within `abs_tol` (e.g. 1e-4) is a **hard gate**. If it fails, stop; the new subsets are uninterpretable.

## Wording map (must appear in the report)

| Do not write | Write |
|---|---|
| N donors independently replicated | N donors, of whom K come from one source block |
| Cross-dataset / cross-source consistent | Consistent inside one source; external source untestable or uninformative |
| LOO N/N robust | Donor LOO N/N; source-block LODO not estimable (or: performed, result …) |

## What source dependence does *not* overturn

- Absence of duplicate donors
- A slope that is not a between-batch gradient (`dataset` already in the design; low between-dataset variance share)
- Within-block stability

It overturns one phrase: **independent replication**.

A Part 6 donor sign test (17/20 Δaxis medians positive under a pinned
Geneformer) is the same class of sentence. It measures direction
consistency inside this object. It does not replace source-block LODO
and it does not upgrade `SINGLE_SOURCE_DEPENDENT`.


