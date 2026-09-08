# Part 5 — Donor-unit association (the statistical spine)

Scientific contract: [research-validity.md](research-validity.md). In particular,
LODO is internal sensitivity and never alone proves independent replication.

Audience: an agent with a **Part 3 locked compartment** and a **Part 4
identifiability table**. Subtype names are frozen. This is the first part
that may use the word **associated**.

It is not a cell Wilcoxon. It is not a new clustering. It does not reopen
Part 3 labels. It does not start because F04_11 looks high.

Stack detail that this chapter operates: [statistical-stack.md](statistical-stack.md),
[source-dependence.md](source-dependence.md),
[evidence-and-verdicts.md](evidence-and-verdicts.md).
Figures: [part5-figures.md](part5-figures.md).

Helpers:

```text
scripts/part5_pseudobulk.py
scripts/part5_source_blocks.py
scripts/part5_eligibility.py
scripts/part5_run_models.R
scripts/part5_run_pathways.R
scripts/part5_verdict.py
scripts/part5_figures.py
```

The R scripts are the gene/pathway engines. Python does not pretend to be
`limma-voomWithQualityWeights`. A Python-only engine is allowed only if the
protocol **names** it and does not claim this stack.

---

## 0. What this part answers

At `dataset × donor_id` (optionally × locked subtype), after declared depth
and confound covariates, is `TARGET_GENE` exposure associated with the rest
of the transcriptome — and what remains after the dominant **source block**
is held out as a whole?

It answers association, enrichment, confound diagnostics, donor leverage,
and source replaceability. It does **not** answer causality, in-vivo time,
or independent replication unless source-block LODO actually ran and passed.

Claim type: associational. Evidence ceiling: `formal` or `exploratory`,
declared **before** looking at q. `RANGE_OK` on F04_15 is a forecast, not
this chapter's verdict.

One Part 5 run = **one estimand** (one population × one exposure × one
outcome family). Another subtype arm, another exposure scale claimed as
discovery, or a clinical contrast is a new protocol, not a renamed CSV.
Matching high vs low cells is a **secondary** track. It never replaces the
continuous primary model.

If a user explicitly requests a cell-level fallback and accepts the exact
string `I_ACCEPT_CELL_LEVEL_FALSE_POSITIVE_RISK`, run
`scripts/part5_cell_exploratory.py`. It must write `EXPLORATORY_ONLY` and a
risk manifest outside the Part 5 verdict/FDR family. It cannot replace a
`NOT_ESTIMABLE` donor result, establish replication, or be described as a
donor-level finding. Without that explicit authorization, refuse execution.

Anti-pattern origin: folder `08` (positive vs undetected cells). That is
not this part. Origin of the stack: `11` v2, `14`, `17`, internals of `29`.
Origin of the source-block gate: `31` — a sequel in the real program, a
**required block of this protocol** in the skill.

---

## 1. Inputs (read-only)

- Part 3 final annotated h5ad. Hash it. Do not write into it.
- Part 4 identifiability CSV (or the protocol names arms a priori and
  says why Part 4 was skipped).
- Frozen subtype color YAML (figures only).
- `layers["counts"]`. Detection, pseudobulk, and models use this layer.
  `layers["normalized"]` is display-only and is not the test.
- Locked subtype key with no NA. `obs`: `dataset`, `donor_id`. Build
  `dataset_donor_id` in memory if absent.
- Frozen local GMT JSON files with SHA-256. Do not download a new library
  to obtain a paragraph.

Assert `TARGET_GENE` is in `var_names`. If not, stop.

Do not recompute HVG, PCA, Harmony, neighbors, or UMAP.
Do not relabel, merge, or split subtypes.
Do not start an arm whose Part 4 flag is `LOW_N` or `NO_EXPOSURE_RANGE`
unless the protocol explicitly overrides the forecast and states the
ceiling as `exploratory`.

---

## 2. Arms (STOP POINT 0 — freeze, then code)

The protocol lists every arm that will be fit. An arm is a row:

| Field | Rule |
|---|---|
| Population | locked subtype name, or a **pre-declared** merge of names |
| Unit | `dataset_donor_id`, or `dataset_donor_id × subtype` if several states share a slope |
| Min cells / unit | 20 (sensitivity 5 / 10 / 50; sensitivity is not a new FDR family) |
| Part 4 flag | quote it; `UNLIKELY_PART5_ARM` stays out of formal discovery |
| Ceiling | `formal` or `exploratory` |

Unresolved / stressed / doublet / debris names do not enter formal arms
just because they have cells.

A merge of subtypes is allowed only if declared here **before** seeing
gene tables. Do not merge after a volcano looks empty.

Human freeze of this table is a stop. The agent does not invent a fifth
arm because a violin was high.

---

## 3. Exposure, unit, eligibility of the *actual* design

### Exposure

Primary (starting freeze):

```text
jeffreys = (n_detected + 0.5) / (n_cells + 1)
jeffreys_per_10pct = 10 × jeffreys
```

`n_detected` is cells with `TARGET_GENE` count > 0 in that unit. The
coefficient is log2 change per **+10 percentage points** of Jeffreys
detection. A raw zero is non-detection in this library, not a biological
off-state.

Alternative scale: unit `log2(CPM+1)` of the same gene, optionally
z-scored on the **full** eligible cohort and then frozen for every
subset (so holdout coefficients stay on one scale).

If Spearman(detection, log2CPM) at eligible units ≥ 0.9, tag
`COLLINEAR_SCALES`. The second model is a scale check. It is not a second
study and cannot lift a source-dependence ceiling.

### Depth and confound

Keep in the design (starting freeze; name them in the protocol):

```text
z_log1p_n_cells
z_log1p_mean_umi          # mean (or median) total UMI per cell, target-excluded
```

Optional, only if the protocol froze the gene lists **before** touching
this matrix: `z_foreign` (environmental / off-lineage panel). Do not build
the panel from genes chosen by looking at the same slope.

z-scores: compute on the full eligible cohort of this arm, stash, reuse
in every LOO / LODO subset. Do not re-standardize inside a holdout.

### Eligibility (formal / exploratory / dead)

Part 4 forecasted n and exposure sd. Part 5 recomputes **residual df of
the design that will actually be fit**.

| Gate | Formal | Exploratory | Dead |
|---|---|---|---|
| Eligible units | ≥ 12 | ≥ 8 | `< 8` → do not fit as discovery |
| Datasets / source_blocks | ≥ 3 | ≥ 2 if protocol allows | 1 unless protocol is explicitly single-source |
| Residual df | ≥ 6 | ≥ 4 | below → `NOT_ESTIMABLE` |
| Design | full rank | full rank | rank-deficient → `NOT_ESTIMABLE` |
| Exposure | IQR / sd large enough to identify a slope | same, or label `NO_EXPOSURE_RANGE` | sd ≈ 0 |

A source with almost no exposure variance cannot rescue a slope, even
with thousands of cells.

Drop a design term when it has a single remaining level. Record the
deviation. Do not keep a constant column. Do not delete covariates until
a number appears.

Run `scripts/part5_eligibility.py` **before** R. If the arm is dead,
write coverage + the gate, skip limma, still emit `verdict.json`.

---

## 4. Target-excluded pseudobulk

From `layers["counts"]` only:

1. One row per unit (`dataset_donor_id`, or × subtype if that is the unit).
2. Sum gene UMIs. Integer counts. No CP10k, no scaled `.X`.
3. Drop `TARGET_GENE` from the **outcome** matrix and from TMM /
   size-factor estimation. Keep its counts in a sidecar as the exposure.
4. Write genes, units, `n_cells`, `n_detected`, `jeffreys_per_10pct`,
   `log2cpm`, depth covariates, `dataset`, `source_block`.
5. `filterByExpr` or CPM ≥ 1 in ≥ max(3, 20% of units), **per arm**,
   after the target is already gone.

Do not:

- pseudobulk on `normalized`
- treat technical libraries as extra units
- put `TARGET_GENE` back into the outcome "so the volcano has a hero gene"
- use cell-level `rank_genes_groups` as the table R is asked to confirm

Command:

```bash
python scripts/part5_pseudobulk.py locked_subtype.h5ad \
    --gene TARGET_GENE --group subtype \
    --unit-key dataset_donor_id --min-cells 20 \
    --out 03_pseudobulk/
```

`genes.csv` is the sole outcome-gene manifest consumed by the R model. It
already excludes `TARGET_GENE`; do not create or substitute a second gene list.

---

## 5. Source blocks (map before any q)

Assign `source_block` from **repository-internal** evidence only
(metadata, donor IDs, consecutive GSE). Do not browse the web to
"confirm" a paper unless the protocol says so.

Default: `source_block = dataset`. Consecutive accessions with
**interleaved** donor numbering are a suggestion, not an automatic merge.
Write `source_block_suggestions.csv`. The protocol must freeze the map.
The agent does not silently collapse two GSE IDs.

Tiers: `SAME_STUDY_LIKELY` | `SAME_STUDY_CONFIRMED` | `UNRESOLVED`. The
holdout design does not change with the tier.

Also write exposure mean/sd/IQR per dataset and per `source_block`. Flag
`NO_EXPOSURE_RANGE` when sd is below the protocol floor (starting freeze:
0.02 on Jeffreys fraction, 0.15 on log2CPM; or 0.5 on
`jeffreys_per_10pct`).

Between-dataset variance share of the exposure is a **confound diagnostic**,
not a replication test. A small share plus `dataset` already in the design
means the slope is not a batch-gradient disguise. It does not mean
independent replication.

```bash
python scripts/part5_source_blocks.py 03_pseudobulk/metadata.csv \
    --map 00_protocol_manifest/source_block_map.yaml \
    --exposure jeffreys_per_10pct \
    --out 02_tables/
```

---

## 6. Gene models

Recommended default (R). Copy the order; do not reorder to chase q.

```text
filter → TMM (target-excluded) → duplicateCorrelation(block=dataset_donor_id)
  → voomWithQualityWeights(block=, correlation=)
  → lmFit(block=, correlation=) → eBayes(robust=TRUE)
```

Typical slope design (amend names, not the idea):

```text
~ subtype + dataset + z_log1p_n_cells + z_log1p_mean_umi + jeffreys_per_10pct
```

When the unit is already `donor × subtype`, the subtype term is still
allowed if several states share one slope. Inside a **single-subtype**
arm, drop `subtype`. Inside a single-dataset subset, drop `dataset`.

If `duplicateCorrelation` cannot run (donors < 3 or no repeated units):
record `NO_BLOCK_FALLBACK`, refit without block, and treat that subset as
exploratory unless the protocol pre-declared the fallback as acceptable
for formal claims.

Support: **edgeR QL**, same units, same filter, same TMM, same design.
`robust_primary` genes need the same **sign**. Support q is not a veto
unless the protocol said it was a gate.

```bash
Rscript scripts/part5_run_models.R path/to/analysis_config.yaml
```

### FDR family (name it)

- Genome-wide arm: BH over all tested genes in **this estimand**. If
  several formal subtypes are concatenated, the protocol must say whether
  `q_state` or `q_project` (all formal gene×state rows) is the discovery
  family. Do not silently pick the friendlier one.
- Panel arm: BH inside the frozen gene list (`q_panel`).
- LOO, LODO, threshold, and alternative-exposure rows **do not** enter
  that family.

### `robust_primary` gene

All of:

1. q < 0.05 in the declared family
2. `|logFC| ≥ 0.10` on the modeled exposure scale
3. support method same sign
4. every **successful** source-block LODO fold same sign (donor LOO is
   additional leverage, not a substitute)
5. not a technical gene (prefixes `MT-`, `RPS`, `RPL`; hemoglobin;
   `MALAT1`; `XIST`)

Technical genes may be measured. They are not strict DEGs.

Non-significant ≠ equivalent. Do not "prove" a bulk signature false
because this slope did not pass q.

---

## 7. Pathways

- Primary: CAMERA on the voom design (competitive).
- Secondary: fgsea multilevel on the primary moderated-t ranking
  (`minSize=10`, `maxSize=500`).
- Libraries: frozen local JSON, SHA-256 checked.
- ORA: only for a pre-specified strict non-technical DEG query of 10–500
  genes and ≤20% of the tested universe; otherwise skip.

`robust_primary` pathway: CAMERA project q < 0.05 **and** fgsea project
q < 0.05, same direction, all successful LODO folds same sign, technical
leading-edge fraction < 0.50.

Identity-like slopes that move thousands of genes often empty CAMERA.
That is expected. Do not switch to a friendlier library to obtain a
paragraph. Cell-level AUCell / ssGSEA q-values are not pathway discovery.

```bash
Rscript scripts/part5_run_pathways.R path/to/analysis_config.yaml
```

---

## 8. Robustness — donor LOO and source-block LODO

Both live in **this** protocol. If LODO was deferred, Part 5 is unfinished.

### Donor LOO

Omit one `dataset_donor_id` at a time on the parent design. Record sign
concordance and the largest absolute effect shift. This measures
**leverage of one person**.

Relabel existing LOO iterations by `source_block`. If ≥75% of iterations
drop a donor from the same block, tag `LOO_SOURCE_UNBALANCED`. Then donor
LOO cannot detect source dependence.

### Source-block LODO (mandatory when more than one study exists)

Refit the **parent design** on:

| Subset | Meaning |
|---|---|
| `FULL` | all eligible units |
| `DROP_EXT` | hold out the smallest / external block |
| `DROP_DOMINANT` | hold out the largest source block as a whole |
| optional split-half | hold out one accession inside a sibling block |

Within each remaining `dataset`, also fit a within-source slope (drop the
`dataset` term). Agreement between two halves of one study is split-half
consistency, **not** cross-source replication.

If residual df dies or the design is rank-deficient: `NOT_ESTIMABLE`.
That is the result. Do not switch to cell Wilcoxon to "recover power".

Reproduction of a parent-stage `FULL` coefficient, when this run subsets
a previous model, is a hard gate (`abs_tol` e.g. 1e-4). Failure →
`REPRODUCTION_FAILED` / `INCONCLUSIVE`. Stop interpreting holdouts.

What a fold must preserve: **sign**, not p-value.

---

## 9. Verdict table (first matching row wins)

### Runner handoff (1.9.3)

The bundled runner executes exactly one frozen arm per stage. It aligns the
source table to the original pseudobulk by `unit_id`, applies the arm labels,
freezes z-scores on that eligible arm, and runs every source-block holdout plus
donor LOO. `robustness.dominant_block` names the dominant block explicitly.
`gene_evidence.csv` has per-gene q/effect, support-direction and LODO checks;
`fold_audit.csv` preserves failed folds. An unrestricted pass requires all
source folds estimable; incomplete folds retain a limited interpretation.

`model_audit.json` derives gene gates from actual outputs. Its gene-family
summary requires all nontechnical q/effect candidates to pass the dual-method
and LODO checks. It does not certify pathways or choose a target gene.
Reproduction is `NOT_APPLICABLE` (null flag) unless `reproduction.required`
is true with a hashed parent coefficient CSV and frozen tolerance/anchors.

Supply project-specific diagnostic booleans in a hashed `audit_flags` JSON
(`confound_cleared`, `sensitivity_pass`, optional `collinear_uninterpretable`).
Missing flags stay unknown, yielding `INCONCLUSIVE` under the example table.
Do not fill them with true simply to obtain a pass. The example sensitivity
caveat row precedes the unrestricted pass and retains every primary gate.

The pathway helper verifies all GMT SHA-256 values itself and adjusts q over
all loaded libraries within each method. It emits `dual_method_candidate`,
with `robust_primary` unset until pathway-level LODO and technical leading-edge
checks are implemented in the stage's integration code. Part 5 remains
unfinished until these declared checks and confound diagnostics are complete.

For input provenance, `hash_inputs.py --verify` compares explicit
`{path, sha256}` entries or sibling `h5ad`/`h5ad_sha256` fields. Relative paths
are stage-relative for configs in `00_protocol_manifest`, otherwise relative
to the manifest. Plain invocation only inventories hashes and does not verify
a freeze. List checkpoint/dictionary files explicitly in an input manifest;
do not treat output paths as frozen inputs.

Declare the ordered table in the protocol before R. Starting freeze for a
genome-wide or panel slope (amend tokens, not the "first hit" rule):

| Order | Token | When |
|---|---|---|
| 1 | `REPRODUCTION_FAILED` | parent `FULL` coefficient not recovered |
| 2 | `NOT_ESTIMABLE` | FULL rank/rdf/exposure gate failed |
| 3 | `COLLINEAR_UNINTERPRETABLE` | protocol asked to separate two scales and could not |
| 4 | `SINGLE_SOURCE_DEPENDENT` | dominant-block holdout `NOT_ESTIMABLE`, sign flip, or CI that admits a large reversal |
| 5 | `PARTIALLY_CONFOUNDED` | primary passes; a pre-registered confound diagnostic does not |
| 6 | `FROZEN_PASS` | formal gates, dual method, LODO sign rules, q and effect floor |
| 7 | `FROZEN_PASS_WITH_SENSITIVITY_CAVEAT` | primary passed; a named sensitivity (threshold, neighbors, PC) did not |
| 8 | `INCONCLUSIVE` | nothing above matched |

Optional tags (may coexist): `WITHIN_SOURCE_CONSISTENT`,
`EXT_NO_EXPOSURE_RANGE`, `LOO_SOURCE_UNBALANCED`, `COLLINEAR_SCALES`.

Do not invent a friendlier English synonym that smuggles extra meaning.
Do not average tokens.

```bash
python scripts/part5_verdict.py 05_logs/model_audit.json \
    --table 00_protocol_manifest/verdict_table.yaml \
    --out 05_logs/verdict.json
```

`verdict.json` required keys: `status`, `verdict`, `reason`, `tags`,
`evidence_ceiling`, `n_units`, `n_donors`, `n_datasets`, `n_source_blocks`,
`rdf`, `drop_dominant_status`.

---

## 10. Directory

```text
analysis/<NN>_<short_name>/
  00_protocol_manifest/   PROTOCOL, analysis_config.yaml, GMT hashes,
                          source_block_map.yaml, verdict_table.yaml
  01_code/                00_run, 01_prepare, 02_run_models.R, 03_integrate
  02_tables/              units, genes, holdout, LOO, CAMERA, fgsea
  03_pseudobulk/          mtx + genes + metadata (target-excluded)
  03_figures/             F05_* PNG/PDF/SVG + source_data/ + sidecars
  05_logs/                input_hashes.json, model_audit.json, verdict.json
  06_reports/             human report; claim boundaries in the same file
  README.md               entry only; not a second protocol
```

Do not overwrite an existing `verdict.json` from a different estimand.
A new contrast is a new folder.

---

## 11. What Part 5 may claim

- In this population, at donor-unit level, exposure is associated with
  gene Y (effect, CI, q, n_units / n_donors / n_datasets / rdf).
- A gene or pathway is `robust_primary` under the frozen definition.
- Dominant-block holdout is `NOT_ESTIMABLE` → `SINGLE_SOURCE_DEPENDENT`.
- Two accessions inside one `source_block` agree in sign → within-source
  / split-half consistency.
- Detection and intensity are collinear at the unit; passing both is not
  two cohorts.

## 12. What Part 5 may not claim

- `TARGET_GENE` drives / causes / regulates Y
- N donors independently replicated (without a successful source-block
  holdout)
- Cross-dataset consistent (when the datasets are one block)
- Donor LOO N/N therefore general
- Cell Wilcoxon / `rank_genes_groups` / MAST as the finding
- Part 4 `RANGE_OK` as a fitted association
- Split-half of one study as external replication
- A second exposure scale as a new study when `COLLINEAR_SCALES`
- Empty CAMERA → try another GMT until a paragraph appears
- Non-significant = equivalent
- Identity DEG tables (Part 7) as "pathways of `TARGET_GENE`"
- Orthogonal bulk / CellChat / STRING as paying this chapter's
  replication debt

---

## 13. Explicitly not done here

- Reclustering or relabeling
- Downloading a new GSE because the report "needs replication"
- Trajectory, CellChat, STRING (Part 8)
- Virtual KO of `TARGET_GENE` (Part 6)
- Frozen module deepdive on a new object (Part 8)
- Clinical grouping / identity contrast (Part 7; own FDR family)
- Post-hoc observed power
- Promoting folder-`08`-style detected vs undetected cells to discovery

---

## 14. Lifecycle

```text
Part 5 progress:
- [ ] Confirm Part 3 lock; read Part 4 identifiability; list arms
- [ ] STOP POINT 0: freeze protocol (estimand, design, GMT hashes,
      source_block map, verdict table, forbidden sentences, ceiling)
- [ ] Hash inputs; abort on mismatch
- [ ] Target-excluded pseudobulk from counts
- [ ] Eligibility of the actual design (rdf, rank, exposure sd)
- [ ] Primary limma-voom QW + support edgeR
- [ ] CAMERA primary + fgsea secondary on hashed GMT
- [ ] Donor LOO
- [ ] Source-block LODO in the same run
- [ ] Collinearity / confound diagnostics
- [ ] verdict.json (first matching row)
- [ ] F05 figures (NOT_ESTIMABLE rows visible; scales not mixed)
- [ ] Report: token + n + rdf + wording consequences
```

Default split: expensive model writes the protocol; cheap model executes
it. Conflict with the protocol → stop and report. Do not silently change
the design to force a number.

