# Evidence ceilings and verdict grammar

## Ceilings (declare before results)

| Token | Meaning |
|---|---|
| `formal` | Protocol eligibility gates met; discovery is conditional on model assumptions, not empirically certified FDR or adequate power |
| `exploratory` | Underpowered, single-source, or audit-only; no project-level discovery narrative |
| `orthogonal_external` | Different system or cohort. Cannot lift a source-dependence ceiling |
| `geometry_only` | Manifold / PAGA / DPT without a donor-level gene model |

A stage header may also set `can_only_downgrade: true`. Then the verdict table may not contain an upgrade of a parent claim.

## Verdict tokens (closed set; extend only in the protocol)

Use the token in `verdict.json` and in the first sentence of the report. Do not paraphrase it into a stronger English claim.

| Token | Use when |
|---|---|
| `FROZEN_PASS` | Protocol executed; formal gates held; dual-method + LODO rules met as specified |
| `FROZEN_PASS_WITH_SENSITIVITY_CAVEAT` | Primary passed; a named sensitivity (neighbors, PC, threshold) did not |
| `PASS_WITH_LIMITATIONS` | Result is real under the spec, with a limitation that forbids a class of sentences |
| `NOT_ESTIMABLE` | Rank, residual df, single level, or no exposure range. **This is the result.** |
| `INCONCLUSIVE` | Reproduction anchor failed, or no verdict row matches |
| `REPRODUCTION_FAILED` | Could not recover the parent-stage coefficient within tolerance |
| `SINGLE_SOURCE_DEPENDENT` | Dominant source-block holdout is `NOT_ESTIMABLE` or flips |
| `ORTHOGONAL_SUPPORT` | A second readout (e.g. intensity vs detection) agrees; still not a new cohort |
| `PARTIALLY_CONFOUNDED` | Association survives the primary design but a pre-registered confound diagnostic does not clear |
| `MODULE_TOO_SPARSE` | Frozen module lost too many genes in this matrix |
| `COLLINEAR_UNINTERPRETABLE` | Score and exposure not separable |
| `NO_TIER1_AXIS` | Communication / screen produced edges but donor-level models did not |
| `EMBEDDING_SHIFT_CONSISTENT` | Part 6: frozen KO primary sign tests passed BH; still exploratory embedding |
| `TOKEN_UNOBSERVABLE` | `TARGET_GENE` never present in the final V2 sequence |
| `SHAM_DRIFT` | unmodified re-inference exceeded the freeze |
| `KO_OE_UNPAIRED` | tag: KO and OE cell sets differ; not a paired mechanism |
| `CONTROL_MAGNITUDE_NOT_SPECIFIC` | tag: absolute rank among matched controls is not distinctive |

Ordered tables: **first matching row wins**. Do not average verdicts.

## Upgrade and downgrade

| Move | Allowed? |
|---|---|
| Hold the parent wording | Yes |
| Replace "independently replicated" with "within one source block" | Yes (downgrade) |
| Lift `exploratory` to `formal` because a p-value was small | No |
| Lift `SINGLE_SOURCE_DEPENDENT` because an external treatment co-induced genes | No |
| Lift `SINGLE_SOURCE_DEPENDENT` because virtual KO sign-test passed | No |
| Call collinear detection and intensity "two independent datasets" | No |
| Call split-half of one study "cross-dataset replication" | No |

Orthogonal evidence can **change the mechanistic story you are allowed to hint** (e.g. parallel pathway targets rather than one target gene causing a downstream gene) without changing the replication ceiling.

## What a tiny p-value does not mean

Inside one cohort, q = 1e-18 measures signal-to-noise **inside that cohort**. Any confound that is constant across the cohort and varies between donors will pass `dataset` covariates, `duplicateCorrelation`, and donor LOO.

That is why source-block holdout exists, and why this skill treats "N donors" as a wording hazard.

## Report shape

```text
1. One-sentence conclusion with the verdict token
2. Reproduction anchor (if this stage subsets a parent model)
3. Eligibility / coverage
4. Primary numbers (effect, CI, q, n_units / n_donors / n_datasets / rdf)
5. Robustness, including source-block if relevant
6. Wording consequences table
7. Explicitly not shown / not claimed
```

Machine file `05_logs/verdict.json` must include at least:

```json
{
  "status": "SUCCESS",
  "verdict": "SINGLE_SOURCE_DEPENDENT",
  "reason": "...",
  "tags": ["WITHIN_SOURCE_CONSISTENT"],
  "evidence_ceiling": "exploratory",
  "can_only_downgrade": true
}
```



Every Part 5 verdict also reports `calibration_status`, `scientifically_calibrated`,
`precision_status` and `protocol_chronology`. A `FROZEN_PASS` token cannot override
unknown calibration or chronology. Legacy external-replication tokens are
rejected by the matcher, including in custom tables.
