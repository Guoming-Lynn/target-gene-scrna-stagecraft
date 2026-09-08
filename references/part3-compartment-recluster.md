# Part 3 — Compartment recluster (rounds, then names)

Audience: an agent with a **Part 1 locked atlas**. Part 2 may have said which
lineage is even on. This chapter subsets that lineage and rebuilds the manifold
from `layers["counts"]`.

The house rule, copied from a production multi-round compartment:

> Delete low-quality / contaminant clusters in **rounds**. After every
> deletion, reset to raw counts and **recompute HVG → PCA → Harmony →
> neighbors → UMAP → Leiden**. Numeric IDs after a recompute are new IDs.
> **Name subtypes only on a round that has zero DELETE.**

Do not inherit the global UMAP. Do not name a cluster after `TARGET_GENE`.
Do not Harmony-correct donor, region, condition, disease, or the gene.

Figures: [part3-figures.md](part3-figures.md).

---

## 0. What this part answers

After extracting a frozen set of Part 1 global types, which **numeric**
clusters are debris, doublets, or cross-lineage soup, and which remaining
states can a human name?

It does **not** answer DEG, trajectory, communication, or whether
`TARGET_GENE` marks a subtype.

Claim type: geometry + labels. Evidence ceiling: descriptive.
Part 4 may start after the lock. Part 5 may not.

---

## 1. Extract (round 00)

From the Part 1 annotated h5ad, **read-only**:

1. Freeze `source_label_key` and the list of global types to keep.
2. Subset those barcodes. Prefix nothing; barcodes are already unique.
3. Copy `layers["counts"]` onto `.X`. Drop every PCA / Harmony / neighbor /
   UMAP / Leiden key from the **compartment** object. Do not keep the global
   embedding as if it were this compartment's map. Provenance only:
   `source_global_cell_type`.
4. Assert integer-like counts, `TARGET_GENE` still in `var_names`.
5. Write `objects/round_00_initial/<stem>_00_raw_counts.h5ad` (gzip).
   This file is the round's immutable counts checkpoint.

`n` for coverage tables is already `dataset × donor_id`.

---

## 2. One round's compute (do not reorder)

```text
raw-count checkpoint
  → CP10k + log1p → layers["normalized"]  (counts stay counts)
  → dataset-aware HVG (Seurat), TARGET_GENE forcibly excluded
  → PCA on those HVGs only (assert the gene is absent)
  → Harmony(batch_key=dataset)   # nothing else
  → neighbors, UMAP, Leiden grid
  → STOP POINT 1  (human picks one resolution)
  → markers + QC + contamination on THAT resolution
  → STOP POINT 2  (human: KEEP or DELETE every cluster; no names)
```

HVG freeze (starting numbers; amend in the protocol):

| Piece | Starting freeze |
|---|---|
| Normalize | target_sum 10,000 + log1p |
| HVG | Seurat, `batch_key=dataset`, rank 3,200, **drop `TARGET_GENE`**, keep 3,000 |
| PCA | compute 50, use 30; sensitivity {20,30,40} |
| Harmony | `dataset` only. Prohibit donor, sample, library, region, condition, disease, source cell type, `TARGET_GENE` |
| Neighbors | 30; sensitivity {15,30,50} |
| UMAP | min_dist 0.3, freeze seed |
| Leiden | grid 0.2–1.2; **human** picks one |
| Markers | Wilcoxon on `normalized`; display = strict-positive (below) |
| Annotation-excluded | `TARGET_GENE` plus any protocol extras |

If HVG ranking still listed the gene, drop it **after ranking** and take the
next genes so PCA sees exactly 3,000. Write that audit row. Assert
`TARGET_GENE not in adata.var.loc[adata.var.highly_variable].index` before PCA.

Stash `obsm["X_pca_unintegrated"]` **before** Harmony. Harmony overwrites
`X_pca`. The before/after dataset plot is illegal if the unintegrated PCs
were never saved.

Harmony on `sample_id` is **not** the Part 3 default (stricter than Part 1).
Only an explicit protocol amendment may add a technical key.

Commands (do not invent a fourth):

```text
run-leiden        → STOP POINT 1 (grid only)
run-markers       → STOP POINT 2 (KEEP/DELETE worksheet)
prepare-removal   → archive DELETE clusters; write child raw-count object
                  → then run-leiden on the child
```

There is no `run-names` on an intermediate round.

---

## 3. STOP POINT 1 — resolution only

Use [checkpoint-review-template.md](checkpoint-review-template.md) and stop if
the reviewer cannot state a biological rationale and alternative interpretation.

Agent exports the Leiden grid and **stops**. No Wilcoxon, no DELETE, no names.

Required: resolution-overview table (n_clusters, min size, n units, max
dataset fraction) plus the figure set in [part3-figures.md](part3-figures.md)
§ STOP POINT 1.

A software hint is allowed (cluster count and min size). It is not a choice.
Do not pick the resolution that makes `TARGET_GENE` look cluster-specific.

Wait for `selected_leiden_resolution`.

---

## 4. STOP POINT 2 — KEEP / DELETE, still no names

Every KEEP/DELETE row must carry reviewer/date, evidence paths, rationale,
alternative interpretation, uncertainty and target-gene blinding status.
Missing or target-gene-only justification is not a valid KEEP.

Only now: Wilcoxon on the selected key.

Strict-positive display table (this is what the human reads):

```text
score > 0
logFC > 0
adj P ≤ 0.05
pct_in > pct_out
pct_in ≥ 0.10
```

Drop `TARGET_GENE` from the annotation worksheet even if it ranks. Rank
p-values are cluster diagnostics, not donor-level findings.

Also export:

- top50 (unfiltered) and top20 / top20-by-cluster (strict-positive)
- per-cluster QC (same columns as Part 1 §6.2)
- contamination / off-lineage panel hits
- a **blank decision CSV** (not a naming YAML)

```text
round, selected_leiden_column, cluster_id, n_cells, decision, reason, reviewer, review_date
```

`decision` is `KEEP` or `DELETE` (aliases: keep/retain, delete/remove). Empty
is illegal. **Every** cluster of the selected resolution must have a row.
**Do not put a cell-type name in `decision` on an intermediate round.**

Human options at STOP POINT 2 (agent does none of these unprompted):

- **A** — every cluster `KEEP` and this is the clean round → go to STOP POINT 3 names.
- **B** — listed clusters `DELETE` with reasons → `prepare-removal`, then `run-leiden`.
- **C** — pick a different resolution on **this** round's already-computed grid
  (`run-markers` again). Do not recompute HVG just to fish a prettier cut.

Flag (agent) / delete (human):

- top20 is MALAT1/NEAT1/MT-/RPS/RPL with zero lineage-panel hits
- median %MT or low complexity clearly apart from the rest
- `max_dataset_fraction` > 0.90 and no lineage hits
- contamination-panel genes dominate over the compartment's own panel

A flagged cluster may still be a real state. Write the choice.

Helper:

```bash
python scripts/cluster_review_tables.py round_leiden.h5ad \
    --leiden-key leiden_r0_5 --out tables/round_00_initial/markers \
    --strict-positive --exclude-genes TARGET_GENE \
    --lineage-genes PTPRC,CD3D,PECAM1,ACTA2,CD79A,LYZ
python scripts/part3_decision_template.py \
    tables/round_00_initial/markers/leiden_r0_5_cluster_qc.csv \
    --round round_00_initial --leiden-key leiden_r0_5 \
    --out tables/round_00_initial/manual_decision_template.csv
python scripts/part3_figures.py stop2 round_leiden.h5ad \
    --leiden-key leiden_r0_5 --panels config/marker_panels.yaml \
    --out figures/round_00_initial
```

---

## 5. Ambiguous cluster — diagnostic subroutine

When a cluster is not obvious debris but is biologically mixed (two lineages,
doublet-shaped, one-donor island):

Do **not** name it. Do **not** lower Leiden resolution to hide it.

Run a named diagnostic (own folder `diagnostics/cluster_<id>/`):

- existing Scrublet scores (do not rerun Scrublet on the concat)
- lineage program scores = mean `normalized` of frozen panels; missing genes listed
- raw-count coexpression of the two conflicting markers
- composition by dataset / donor / library
- optional: the same Harmony graph at neighboring resolutions (Jaccard of
  membership) — do not recompute a new graph just to fish

Then the human still returns KEEP or DELETE. This is what a "cluster 3
audit" folder is for. It is not a new part.

---

## 6. DELETE → archive → new round → full HVG

If **any** cluster is `DELETE`:

1. Decision file must name unique cluster IDs, non-empty reasons, matching
   `selected_leiden_column`, and cover every cluster. Hash the file.
2. Partition the **raw-count checkpoint** of this round (not the normalized
   view). Integrity: removed + retained = parent, empty intersection.
3. Archive `objects/<round>/removed_clusters_<ids>.h5ad` plus barcode CSV
   (dataset, donor, library, UMI, n_genes, scrublet, old cluster, reason).
4. Child round name: `round_{nn+1:02d}_after_removal` (or `_after_low_qc_removal`).
5. On retained cells: `.X = counts`, drop `normalized`, drop PCA / Harmony /
   neighbors / UMAP / Leiden. Print `no old embedding is reused.`
6. Re-run **the entire §2 stack**, including **new HVG**. Never subset the
   old UMAP. Never reuse old HVGs.
7. Write a parent→child partition manifest now. Write the old-cluster vs
   new-cluster membership table **after** the child Leiden exists (STOP POINT 1
   of the child). New IDs are not the old IDs.
8. Return to **STOP POINT 1**. Resolution may stay "0.5" as a number; the
   clusters are new.

Overwrite policy: **forbid**. A child round that already exists is an error,
not a clobber.

```bash
python scripts/part3_prepare_removal.py \
    --parent-raw objects/round_00_initial/raw_counts.h5ad \
    --parent-clustered objects/round_00_initial/leiden.h5ad \
    --leiden-key leiden_r0_5 \
    --decision tables/round_00_initial/manual_decision.csv \
    --removed-out objects/round_00_initial/removed_clusters.h5ad \
    --child-raw objects/round_01_after_removal/raw_counts.h5ad \
    --tables-out tables/round_00_initial
python scripts/part3_round_audit.py \
    --parent objects/round_00_initial/leiden.h5ad --parent-key leiden_r0_5 \
    --child objects/round_01_after_removal/leiden.h5ad --child-key leiden_r0_5 \
    --out tables/round_01_after_removal/old_to_new_membership.csv
```

If the human keeps a messy cluster, freeze it only at the **final naming**
step (e.g. `Unresolved / stressed`). Do not delete it later by stealth
inside Part 4/5.

---

## 7. Stop the lineage

If **two** cleanup rounds still leave mixed-lineage or QC-dominated clusters
that the human will not KEEP as a named state: **stop**. Archive the last
object as exploratory. Do not invent trajectory, CellChat, or DEG to rescue
it. Do not keep deleting until a pretty UMAP appears.

That is an allowed terminal result of Part 3.

---

## 8. STOP POINT 3 — names, only when DELETE = 0

Use the checkpoint template for the final mapping. Record independent review and
adjudication; unresolved disagreement is `STOPPED`, not an invented subtype.

A round may be named only if the filled decision file has **no** `DELETE`
rows (or the human explicitly writes `NAME_NOW` after accepting residual
mess as a labeled state).

Then:

- mapping YAML keys = remaining Leiden IDs of **this** round. Extra/missing → abort.
- Do not copy names from a previous round's cluster numbers. Those IDs died
  with the recompute.
- Variant (must be written in the protocol): propagate **cell-wise** labels
  from the last reviewed round onto retained barcodes, and keep the new Leiden
  as a numeric audit track only. Default is **do not propagate**.
- `TARGET_GENE` does not appear in any label string.
- New gzip h5ad. Re-open and assert: no NA labels, counts + normalized present,
  gene still in `var_names`.
- Frozen color YAML. Lock UMAP (numbers on data, names in legend).

Downstream stages do not reopen these labels.

---

## 9. Directory (per round)

```text
analysis/<NN>_<compartment>/
  00_protocol_manifest/     PROTOCOL, extraction list, HVG/Harmony freeze, panels
  objects/round_00_initial/
  tables/round_00_initial/
  figures/round_00_initial/ 01_HVG … 09_QC  (see part3-figures.md)
  reports/round_00_initial/ STOP_POINT_1_*.md  STOP_POINT_2_*.md
  objects/round_01_after_removal/
  …
  objects/final/<stem>_final_annotated.h5ad
  config/subtype_colors.yaml
```

Sibling numbered folders (`12` vs `12_diag`) are different estimands. A
diagnostic archive is not the locked object.

---

## 10. Second-tissue merge (variant)

New libraries still go through **Part 1 QC**. Then this chapter on the concat.

Do not Harmony on tissue. Do not treat the merge as independent replication
of the index atlas. Legacy labels, if kept, are `legacy_*` only — the working
key is the new Leiden until STOP POINT 3.

---

## 11. What Part 3 may / may not claim

May: these barcodes were extracted from named Part 1 types; they passed N
cleanup rounds with archived reasons; remaining states were human-named from
lineage markers; `TARGET_GENE` was excluded from the manifold.

May not: `TARGET_GENE` defines subtype X; adjacent tissue is healthy;
N cells = n; old Leiden 3 = new Leiden 3; SoupX was done; the lineage is
ready for DEG.

---

## 12. Lifecycle

```text
Part 3 progress:
- [ ] Freeze extraction types, lineage + contamination panels, HVG/Harmony numbers
- [ ] Round 00 raw-count checkpoint from Part 1 counts
- [ ] Compute: normalize, HVG minus TARGET_GENE, PCA, Harmony(dataset), Leiden grid
- [ ] STOP POINT 1: human picks one resolution
- [ ] Markers (strict-positive) + cluster QC + contamination + blank KEEP/DELETE CSV
- [ ] STOP POINT 2: every cluster KEEP or DELETE with a reason; no names
- [ ] If DELETE: archive → child round from counts → new HVG → back to STOP POINT 1
- [ ] If two rounds still mixed: STOP the lineage
- [ ] STOP POINT 3: complete mapping YAML → final gzip h5ad → lock UMAP
```

