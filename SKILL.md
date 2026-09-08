---
name: target-gene-scrna-stagecraft
description: >-
  Rigorous multi-study scRNA-seq workflow for a pre-specified target gene:
  library-level QC, strict common-gene counts, human-locked clustering,
  descriptive target-gene surveys, compartment cleanup, donor-unit association
  with limma/edgeR and source-block LODO, and pinned Geneformer virtual KO/OE
  as an embedding-shift analysis. Use for multi-library target-gene atlas,
  subtype survey, donor association, or in-silico perturbation requests.
  Triggers include 目标基因, 多GEO, Scrublet, MAD, 共同基因, Harmony, Leiden,
  分室重聚类, 供体关联, 伪批量, limma, LODO, Geneformer, and virtual KO.
license: MIT
metadata:
  version: "2.0.0"
  scope: "Parts 1–5 main path; Part 6 separately specified Geneformer branch."
---

# Target-gene scRNA stagecraft

New executor: read only the [stage route and project arm checklist](references/start-here.md)
and the one specification for the selected Part. Do not preload every reference.
Before freezing any analysis, read [research-validity.md](references/research-validity.md).
This is the authoritative scientific contract for estimands, donor independence,
source sensitivity, selection, QC/annotation audits and Part 6 inference.
Legacy LODO success never by itself establishes independent replication.
Before Part 5 validation or a new freeze, read
[calibration-and-provenance.md](references/calibration-and-provenance.md).
`formal` is a protocol eligibility category; empirical scientific calibration
is incomplete. Report calibration, power and chronology status explicitly.
Before execution, consult [glossary.md](references/glossary.md) for unfamiliar
terms in the selected stage. Resolve uncertainty about biological units,
estimands and mandatory gates before proceeding; do not guess their meaning.

For formal figures, also read [figure-contract.md](references/figure-contract.md),
[figure-statistics-contract.md](references/figure-statistics-contract.md), and
[visual-qa-contract.md](references/visual-qa-contract.md). Figures must expose
the declared estimand, unit, statistics and claim ceiling.

The cell-level fallback is opt-in only and requires the exact user authorization
`I_ACCEPT_CELL_LEVEL_FALSE_POSITIVE_RISK`. It is exploratory and cannot enter
the donor-level verdict or discovery family.

This skill is for a **known gene**, called `TARGET_GENE` in every protocol.
You are not discovering what matters. You are asking whether that gene can be
measured honestly in a multi-study tissue atlas, and which cell states even
express it.

It is not a Scanpy tutorial. Use `scanpy` / `anndata` for APIs.
It is not a genome-wide DEG engine. Do not skip to Wilcoxon-on-cells because
someone is impatient.

If a frozen protocol in the current repo conflicts with this skill, the protocol wins.

## Who this is for

The user can name:

1. One or a few **pre-specified genes** (symbol + species).
2. A tissue / disease / experimental context.
3. Two or more scRNA-seq libraries (often several GEO series).

Then the agent builds a **counts-preserving, donor-aware atlas**, locks human
labels, and — if asked for the gene — draws the Part 2 survey. Expression of
`TARGET_GENE` is descriptive until Part 5.

Do **not** use this skill to:

- fish for "the most interesting gene in the atlas"
- treat 10x barcodes as biological replicates
- name clusters after `TARGET_GENE`
- start CellChat, trajectory, or donor-level DEG before labels are locked
- treat a Geneformer cosine as an expression fold-change

## Pipeline map

Six public parts. Parts 7–8 are future claim types, not part of this release.
The map, folder-to-part table, stop rules,
and anti-patterns: [references/pipeline-map.md](references/pipeline-map.md).

| Part | Job | Status |
|---|---|---|
| **1** | QC + global clustering + human labels | **written** |
| **2** | Global `TARGET_GENE` survey (descriptive; no DEG) | **written** |
| **3** | Compartment recluster, contamination rounds, cluster diagnostic | **written** |
| **4** | Subtype `TARGET_GENE` survey (descriptive; no DEG) | **written** |
| **5** | Donor-unit association, enrichment, **source-block LODO in-protocol** | **written** |
| **6** | Virtual KO of `TARGET_GENE` (pinned Geneformer; embedding shift) | **specified / not turnkey** |

Spine is 1→2→3→4→5. Parts 6–8 are branches: they may hold or downgrade earlier
wording; they may not upgrade it. Part 3 is a **loop** (each lineage, or a
second-tissue merge). Do not execute a later part's rules as if they were Part 1.

Part 1 spec: [references/part1-qc-and-global-atlas.md](references/part1-qc-and-global-atlas.md)  
Part 1 figures: [references/part1-figures.md](references/part1-figures.md)  
Part 2 spec: [references/part2-target-gene-survey.md](references/part2-target-gene-survey.md)  
Part 2 catalog: [references/part2-figures.md](references/part2-figures.md) — this is what
the agent is allowed to draw for `TARGET_GENE`. Do not substitute `sc.pl.umap`.  
Part 3 spec: [references/part3-compartment-recluster.md](references/part3-compartment-recluster.md)  
Part 3 figures: [references/part3-figures.md](references/part3-figures.md) — stricter than
Part 1 Checkpoint B; names only after the last DELETE round.  
Part 4 spec: [references/part4-subtype-survey.md](references/part4-subtype-survey.md)  
Part 4 catalog: [references/part4-figures.md](references/part4-figures.md) — F02 geometry
on locked subtypes, plus an identifiability forecast. Do not fit a slope.  
Part 5 spec: [references/part5-donor-association.md](references/part5-donor-association.md)  
Part 5 figures: [references/part5-figures.md](references/part5-figures.md) — donor-unit
volcano, source-block holdout, `NOT_ESTIMABLE` empty slots. Do not Wilcoxon-on-cells.  
Part 6 spec: [references/part6-virtual-knockout.md](references/part6-virtual-knockout.md)  
Part 6 figures: [references/part6-figures.md](references/part6-figures.md) — donor Δaxis,
no fake KO–OE pairing, axes say embedding shift. Do not predict expression.

## Non-negotiables (Part 1)

1. **Protocol before code.** Freeze `TARGET_GENE`, the library manifest, QC numbers, and Harmony keys. Then run.
2. **One library at a time** through QC. Never estimate MAD or Scrublet on a concatenated soup of studies.
3. **Doublets before MAD.** Scrublet on a *loose* filter; final MAD only on singlet candidates.
4. **Scrublet failure ≠ all singlets.** Failed or too-small libraries leave the primary atlas. Log them.
5. **`n` is already `dataset × donor_id`.** Technical runs (Seq1/Seq2) share a donor. Adjacent/non-index tissue is not a healthy control unless the manifest says so.
6. **Counts stay counts.** This chapter does not overwrite raw UMI. `layers["counts"]` is immutable after QC. Clustering copies from that layer.
7. **Strict common genes for the atlas; union is not for HVG.** Outer-join zeros are missing features, not biological zeros. Keep a gene-presence table.
8. **Assert `TARGET_GENE` is in the common set** before clustering. If it is absent, stop. Do not swap in a homolog or a "similar" gene.
9. **Do not cluster on the story.** Do not Harmony-correct donor, region, condition, disease, or `TARGET_GENE`. Do not name Leiden clusters by `TARGET_GENE`.
10. **Delete-and-recompute.** If a cluster is debris/unresolved, archive those barcodes, reload retained cells from QC counts, and rerun HVG→Leiden. Never reuse the old UMAP.
11. **Leiden is a file-gated human checkpoint.** The agent runs the resolution grid and stops. It may continue only when a non-empty `resolution_choice.yaml` exists with reviewer/date/evidence. Wilcoxon top-20 and per-cluster QC are for the **selected** resolution only.
12. **Every cluster gets a name or a DELETE.** The agent may continue only when the completed KEEP/DELETE CSV has one row per cluster, reviewer/date, evidence and rationale, with no blanks. Mapping keys must equal remaining Leiden IDs. The annotated h5ad is a new file.
13. **No DEG, no communication, no trajectory** until the global label key is locked.

### Implementation aids and provenance

The recommended implementation pattern is modular: calculate QC metrics,
detect outliers, apply hard caps, filter cells, and export an attrition summary
as separate steps. This follows the public `single-cell-rna-qc` skill, while
the ordering and thresholds above remain authoritative. Each library must
produce a machine-readable QC record containing its input hash, metric
parameters, Scrublet status/threshold, singlet count, MAD bounds, hard-cap
results, and final retained count. A failed Scrublet call must be represented
as a failed record and cannot silently become an all-singlet result.

Keep the integer `layers["counts"]` layer immutable and validate its type and
shape before normalization or integration. This is an implementation pattern
adapted from `scvi-tools`; a learned latent space remains an optional branch
and cannot replace the frozen atlas or donor-level estimand.

For future FASTQ-scale execution, use an environment preflight, a tiny test
profile, pinned tool/container versions, validated samplesheets, and resumable
work directories. These controls are adapted from `nextflow-development`;
resuming may reuse completed work but may not overwrite a frozen output.

The exact upstream commit, license, and adaptation boundaries are recorded in
[references/upstream-life-sciences-influences.md](references/upstream-life-sciences-influences.md).

## Part 1 lifecycle

```text
Part 1 progress:
- [ ] Freeze TARGET_GENE, species, manifest (dataset, donor_id, library_id, analysis_role)
- [ ] Per-library read + gene-symbol collapse + integer-count check
- [ ] Loose filter → Scrublet (automatic threshold) → singlet MAD-QC + hard caps
- [ ] Checkpoints hashed (pipeline_version + config_signature)
- [ ] Merge: all-context union / primary union / strict-common counts / optional downsample
- [ ] TARGET_GENE present in strict-common genes (assert)
- [ ] From counts: CP10k+log1p, dataset-aware HVG, PCA, Harmony, neighbors, UMAP, Leiden grid
- [ ] Checkpoint A: human picks one resolution from the overview (agent does not)
- [ ] Export top50/top20 + blank marker worksheet + per-cluster QC table
- [ ] Checkpoint B: human names every cluster or marks DELETE; no blanks
- [ ] If DELETE: archive barcodes, recompute from QC counts, return to A/B
- [ ] Complete mapping YAML → new annotated h5ad (gzip) → re-open verify
- [ ] Lock UMAP (numbers on data, names in legend) PNG/PDF/SVG + sidecar
```

Default split: expensive model writes the protocol; cheap model executes it.
Conflict with the protocol → stop and report.

## Part 2 lifecycle

Run only on a Part 1 locked atlas. Full catalog:
[references/part2-figures.md](references/part2-figures.md).

```text
Part 2 progress:
- [ ] Hash locked h5ad; read-only; colors from Part 1 YAML
- [ ] Assert TARGET_GENE in var_names
- [ ] Group summary + unit × type table (min 20 cells)
- [ ] Render F02_01–14 (UMAP, lollipop, violin+box, donor-unit grid, LODO, depth)
- [ ] One-page survey: which types are on, which endpoints disagree
- [ ] Recommend a Part 3 lineage — or none
```

Hero figure is the **donor-unit four-endpoint grid**, not a cell violin.
No DEG. No new embedding. No `sc.pl.*` as the formal panel.

## Part 3 lifecycle

Run only after Part 1 labels are locked. Part 2 may say which lineage is even
on; it is not required if the human already named the compartment.

House rule: **DELETE low-quality / contaminant clusters in rounds. After every
deletion, recompute HVG from retained counts. Name subtypes only on a round
with zero DELETE.**

```text
Part 3 progress:
- [ ] Freeze extraction types, lineage + contamination panels, HVG/Harmony numbers
- [ ] Round 00 raw-count checkpoint from Part 1 counts (drop global embeddings)
- [ ] Compute: normalize, HVG minus TARGET_GENE, PCA, stash unintegrated PCs,
      Harmony(dataset only), Leiden grid
- [ ] STOP POINT 1: human picks one resolution (agent does not)
- [ ] Markers (strict-positive) + cluster QC + contamination + blank KEEP/DELETE CSV
- [ ] STOP POINT 2: every cluster KEEP or DELETE with a reason; no names
- [ ] If DELETE: prepare-removal → child round from counts → new HVG → STOP POINT 1
- [ ] If two rounds still mixed: STOP the lineage
- [ ] STOP POINT 3: complete mapping YAML → final gzip h5ad → lock UMAP
```

Do not inherit the global UMAP. Do not Harmony-correct donor. Do not put a
cell-type name in an intermediate decision CSV. New Leiden IDs after a recompute
are not the old IDs.

## Part 3 non-negotiables

1. **Rounds first, names last.** Intermediate worksheets are KEEP/DELETE only.
2. **New HVG after every deletion.** Reload retained `layers["counts"]`. Never subset an old UMAP.
3. **Harmony is `dataset` only.** `TARGET_GENE` is out of HVG, PCA, Harmony, neighbors, UMAP, Leiden, and the annotation worksheet.
4. **STOP POINT 1 then 2.** The agent does not pick the resolution, delete clusters, or name them.
5. **Stop the lineage** if two cleanup rounds still leave mixed-lineage / QC debris. Do not invent trajectory or DEG to rescue it.

## Part 4 lifecycle

Run only on a Part 3 **named** lock. A stopped lineage does not get a survey.

```text
Part 4 progress:
- [ ] Hash Part 3 h5ad; read-only; subtype colors from Part 3 YAML
- [ ] Assert TARGET_GENE in var_names; no NA on the subtype key
- [ ] Render F04_01–14 (Part 2 geometry, this compartment's UMAP)
- [ ] Identifiability table + F04_15 (subtype × source)
- [ ] One-page survey: which subtypes are on, which slices cannot identify a slope
- [ ] Recommend Part 5 arms — or none
```

Hero remains the **donor-unit four-endpoint grid**. The extra job is the
identifiability forecast. No DEG. No `associated`. No new embedding.

## Part 5 lifecycle

Run only on a Part 3 **named** lock. Prefer a Part 4 identifiability table.
Do not start because F04_11 looks high.

House rule: **source-block LODO lives in this protocol**. If it was deferred,
the association is unfinished. `NOT_ESTIMABLE` is a result.

```text
Part 5 progress:
- [ ] Confirm Part 3 lock; read Part 4 identifiability; list arms
- [ ] Freeze project-wide arm inventory; disclose all runs, including failures
- [ ] STOP POINT 0: freeze protocol (estimand, design, GMT hashes,
      source_block map, verdict table, forbidden sentences, ceiling)
- [ ] Hash inputs; abort on mismatch
- [ ] Target-excluded pseudobulk from counts
- [ ] Eligibility of the actual design (rdf, rank, exposure sd)
- [ ] Primary limma-voom QW + support edgeR
- [ ] CAMERA primary + fgsea secondary on hashed GMT
- [ ] Donor LOO
- [ ] Source-block LODO in the same run
- [ ] Optional cell fallback: authorization recorded and separate exploratory family
- [ ] Collinearity / confound diagnostics
- [ ] Freeze estimand mode and repeated-donor structure
- [ ] Record QC attrition and annotation-review evidence
- [ ] verdict.json (first matching row)
- [ ] F05 figures (NOT_ESTIMABLE rows visible; scales not mixed)
- [ ] Figure statistics manifest + visual QA at final size
- [ ] Report: token + n + rdf + wording consequences
```

Python helpers build units and figures. R is the gene/pathway engine unless
the protocol names another stack and does not claim this one.

## Part 5 non-negotiables

1. **`n` is `dataset × donor_id`.** Cell Wilcoxon is a screen, never discovery.
2. **Target-excluded counts.** `TARGET_GENE` is the exposure sidecar, not an
   outcome gene and not a TMM factor.
3. **Eat the Part 4 forecast.** `LOW_N` / `NO_EXPOSURE_RANGE` arms are not
   formal. `RANGE_OK` is not a fitted pass.
4. **Dual engine.** limma-voom QW primary; edgeR QL support; CAMERA then fgsea.
   `robust_primary` needs the same sign.
5. **Donor LOO is not source-block LODO.** Both run here. Consecutive GSE IDs
   with interleaved donor numbers are one block until the protocol says otherwise.
6. **First matching verdict row wins.** Do not relax the design to force a number.
7. **LODO is internal sensitivity.** It does not establish independent replication;
   a genuinely held-out cohort is a separate validation claim.
8. **Full rank is necessary, not sufficient.** Record VIF, exposure/depth
   correlations, and condition number before interpreting coefficients.
9. **Cell-level fallback requires explicit authorization.** It is a dependent-
   subsample screen with elevated false-positive risk and cannot upgrade Part 5.

## Defaults (freeze, then change only by amendment)

Documented in [references/part1-qc-and-global-atlas.md](references/part1-qc-and-global-atlas.md).
Do not treat them as universal biology. Treat them as a **starting freeze**.

| Piece | Starting freeze |
|---|---|
| Loose filter | genes ≥ 200, UMI ≥ 500 |
| Scrublet | auto threshold; expected rate 0.05; min 200 cells |
| MAD | log1p genes/UMI: 3 / 4 MAD; MT/HB 3 MAD |
| Hard caps | genes 300–7500, UMI ≥ 500, MT ≤ 15%, HB ≤ 3% |
| Atlas genes | inner join of genes present in every primary library |
| Normalize | target_sum 10,000 + log1p; keep `layers["counts"]` and `layers["normalized"]` |
| HVG | Seurat, `batch_key=dataset`, 3,000 genes |
| PCA | compute 50, choose among {20,30,40,50} |
| Harmony | `dataset` required; `sample_id`/`library_id` allowed if technical splits exist |
| Neighbors / UMAP | 30 neighbors, min_dist 0.3 |
| Leiden | grid 0.2–1.5; **human** picks one after overview |
| Seed | freeze it |

High UMI, ribosome, top50, and complexity are **flags**, not deletion rules, unless the protocol says so.

## What Part 1 may claim

- These libraries passed a stated QC and doublet rule.
- The atlas is built on genes actually shared across primary libraries.
- `TARGET_GENE` is measured in that gene set (or the stage stopped).
- Global clusters were named from lineage markers, not from `TARGET_GENE`.
- A stated fraction of cells have `TARGET_GENE` count > 0 (descriptive).

## What Part 1 may not claim

- `TARGET_GENE` marks cell type X
- `TARGET_GENE` is differentially expressed (no donor-unit model yet)
- Adjacent tissue = healthy control
- N cells = N replicates
- Union-gene zeros = the gene was off
- SoupX / CellBender was done (GEO filtered matrices usually cannot support that claim)
- Independent replication across studies (source-block audit is a Part 5 gate)

Forbidden-sentence lint for later parts still lives in [references/claim-boundaries.md](references/claim-boundaries.md); for Part 1, the list above is enough.

## Part 6 runtime contract

Use a separate Geneformer environment. Read [part6-runtime.md](references/part6-runtime.md).
Before each stage run `scripts/check_environment.py --stage partN` with the intended
interpreter. A nonzero exit blocks that stage. This dependency check does not
replace model hash, official parity, human annotation or scientific gates.

## Part 6 lifecycle

Run only on a Part 3 **named** lock. Endpoints hashed before any CLS.
Quote the parent Part 5 token; do not upgrade it.

House rule: **this is an embedding shift**. KO uses final-token-present
cells. OE, if declared, is a separate population. Official Geneformer is
the engine; smoke vs `InSilicoPerturber` is a gate, not a suggestion.

```text
Part 6 progress:
- [ ] Confirm Part 3 lock; quote Part 5 token or PART5_NOT_RUN
- [ ] STOP POINT 0: subtype, endpoint hashes, KO±OE, model pin, family size
- [ ] Hash h5ad / checkpoint / dictionaries
- [ ] Map counts → Ensembl; TARGET_GENE unique token
- [ ] Endpoint coverage; drop TARGET_GENE from members
- [ ] Token ledger (truncation ≠ biology; illegal mismatch STOPPED)
- [ ] Axis feasibility before inference
- [ ] Freeze matched controls
- [ ] Freeze axis sensitivity variants and negative controls
- [ ] Smoke (official cosine ≤1e-5; OE overflow reference)
- [ ] Original CLS cache; leave-one-donor-out axes
- [ ] Formal KO (±OE); sham
- [ ] Donor eligibility ledger BEFORE summaries
- [ ] Closed-family sign tests; BH with declared family size
- [ ] verdict.json; tag KO_OE_UNPAIRED when populations differ
- [ ] F06 figures (no pairing lines; Δaxis not expression)
- [ ] Figure statistics manifest + grayscale/visual QA
```

## Part 6 non-negotiables

1. **CLS shift, not RNA.** No predicted expression, OCR, ATP, or CRISPR.
2. **KO population = token-present.** Do not force truncated tokens back in.
3. **OE is optional and unpaired** unless the protocol froze the same cells.
4. **Pinned official Geneformer.** No fallback engine. One smoke repair.
5. **Donor medians, equal weight.** No cell-level p-values. BH family does
   not shrink when a cell is NA.
6. **Cannot upgrade Part 5.** Sign test N/N is not source-block replication.

## Later chapters (do not run yet)


## Validators

Environment setup and synthetic tests: [README.md](README.md#dependencies-and-validation).
Run one frozen Part 5 arm per config/stage. The R runner produces measured
gene evidence, all source-block holdouts, and donor LOO. Project-specific
confound diagnostics and pathway robustness must still be completed before
claiming the entire Part 5 stage finished; see its specification.

```bash
python scripts/check_protocol.py path/to/PROTOCOL.md
python scripts/check_stage_layout.py path/to/stage_dir
python scripts/hash_inputs.py path/to/analysis_config.yaml --verify
python scripts/cluster_review_tables.py clustered.h5ad \
    --leiden-key leiden_r0_5 --out tables/Leiden_markers \
    --lineage-genes ACTA2,PECAM1,CD3D,CD79A,LYZ
python scripts/part2_figures.py locked_atlas.h5ad \
    --gene TARGET_GENE --group cell_type \
    --colors config/celltype_colors.yaml --out figures/part2
python scripts/cluster_review_tables.py round_leiden.h5ad \
    --leiden-key leiden_r0_5 --out tables/round_00_initial/markers \
    --strict-positive --exclude-genes TARGET_GENE
python scripts/part3_decision_template.py \
    tables/round_00_initial/markers/leiden_r0_5_cluster_qc.csv \
    --round round_00_initial --leiden-key leiden_r0_5 \
    --out tables/round_00_initial/manual_decision_template.csv
python scripts/part3_prepare_removal.py \
    --parent-raw objects/round_00_initial/raw_counts.h5ad \
    --parent-clustered objects/round_00_initial/leiden.h5ad \
    --leiden-key leiden_r0_5 \
    --decision tables/round_00_initial/manual_decision.csv \
    --removed-out objects/round_00_initial/removed_clusters.h5ad \
    --child-raw objects/round_01_after_removal/raw_counts.h5ad \
    --tables-out tables/round_00_initial
python scripts/part3_figures.py stop1 clustered.h5ad --out figures/round_00_initial
python scripts/part4_figures.py locked_subtype.h5ad \
    --gene TARGET_GENE --group subtype \
    --colors config/subtype_colors.yaml --out figures/part4
python scripts/part5_pseudobulk.py locked_subtype.h5ad \
    --gene TARGET_GENE --group subtype --out 03_pseudobulk/
python scripts/part5_source_blocks.py 03_pseudobulk/metadata.csv \
    --map 00_protocol_manifest/source_block_map.yaml \
    --out 02_tables/
python scripts/part5_eligibility.py 02_tables/metadata_with_source_block.csv \
    --out 02_tables/eligibility.csv
Rscript scripts/part5_run_models.R 00_protocol_manifest/analysis_config.yaml
Rscript scripts/part5_run_pathways.R 00_protocol_manifest/analysis_config.yaml
python scripts/part5_verdict.py 05_logs/model_audit.json \
    --table 00_protocol_manifest/verdict_table.yaml \
    --out 05_logs/verdict.json
python scripts/part5_cell_exploratory.py cells.csv \
    --group subtype --value log1p_expr --group-a A --group-b B \
    --authorization I_ACCEPT_CELL_LEVEL_FALSE_POSITIVE_RISK \
    --out 02_tables/exploratory/cell_level_result.json
python scripts/simulation_contract.py --seed 20260906 \
    --replicates 1000 --out 00_protocol_manifest/simulation_manifest.json
python scripts/part5_figures.py holdout 02_tables/gene_effects.csv \
    --gene TARGET_GENE --out 03_figures
python scripts/part6_endpoints.py --sets endpoints.yaml \
    --model-genes model_visible_genes.txt --target TARGET_GENE \
    --out 00_input_audit/endpoint_coverage.csv
python scripts/part6_token_audit.py 02_tables/token_ledger.csv \
    --out 02_tables/token_audit.csv
python scripts/part6_sign_tests.py 02_tables/donor_effects_eligible.csv \
    --target TARGET_GENE --family-size 2 --out 02_tables/sign_tests.csv
python scripts/part6_verdict.py 05_logs/model_audit.json \
    --table 00_protocol_manifest/verdict_table.yaml \
    --out 05_logs/verdict.json
python scripts/part6_figures.py donors 02_tables/donor_effects.csv \
    --out 03_figures
python scripts/_part6_smoke.py
```

Also: [references/part1-figures.md](references/part1-figures.md),
[references/part2-figures.md](references/part2-figures.md),
[references/part3-figures.md](references/part3-figures.md),
[references/part4-figures.md](references/part4-figures.md),
[references/part5-figures.md](references/part5-figures.md),
[references/part6-figures.md](references/part6-figures.md)

Figure QA helper:

```bash
python scripts/validate_figure_manifest.py figures/F05_08.parameters.statistics.json
```

For a formal figure, pass the complete statistics mapping to
`plotting_style.save_figure(..., statistics=...)`; it writes
`<stem>.parameters.statistics.json`, which is the input to this validator.
Rendering-only `parameters` are provenance metadata and are not a statistics
manifest.
