# Part 1 — QC and global atlas

Audience: an agent building a multi-study scRNA-seq atlas around a **pre-specified**
`TARGET_GENE`. No disease mechanism, no cluster named after the gene, no cell-level
p-value as a finding.

Fill `TARGET_GENE`, species, and a library manifest before any compute.

---

## 0. What this part answers

Can we construct a counts-preserving atlas from several studies such that:

- every primary library was QC'd the same way
- doublets were called before library-specific MAD thresholds
- the gene universe used for clustering is **really shared**
- `TARGET_GENE` is in that universe
- global cell types are named from lineage markers by a human

It does **not** answer whether `TARGET_GENE` is a marker, a driver, or DE.

Claim type: data-foundation + geometry + labels. Evidence ceiling: descriptive.

---

## 1. Manifest before matrices

Every library (usually one GSM) has at least:

```text
dataset, gse, gsm, library_id, sample_id, donor_id,
technical_run_id, possible_technical_replicate,
tissue, region, condition, technology,
analysis_role, include, include_in_primary_atlas,
source_matrix_level, note
```

Rules:

- `donor_id` is the biological person (or animal). Seq1/Seq2 of the same person share `donor_id` and keep distinct `library_id`.
- `analysis_role` examples: `primary`, `paired_context`, `gex_only_cite`, `reference_extension`, `optional_off`.
- Paired adjacent / proximal / contralateral tissue may be kept for context. It does **not** enter `include_in_primary_atlas` unless the protocol says it is the index tissue. Never call it a healthy control without a frozen definition.
- CITE-seq: default GEX only; ADT is a different assay.
- Optional series stay off until someone has read the clinical/experimental context.
- Clinical fields you did not parse stay `unknown`. Do not guess sex, age, or symptom from a title string.

Write the manifest, hash it, and treat it as part of the config signature.

---

## 2. Read, gene names, counts contract

Per library:

1. Read 10x MTX/H5 or a gene×cell text matrix into cells × genes.
2. Prefix barcodes with `library_id` so later concat cannot collide.
3. Map to gene symbols. If a features file marks duplicates with trailing `.1` / `-1`, strip only when that mark is explicit, then **sum counts** for the same base symbol.
4. Keep an audit table of raw symbol → collapsed symbol. Do not rely on `var_names_make_unique()` as the only record.
5. Assert the matrix is non-negative **integer-like** UMI (sample the nonzero values). If it looks log-normalized, stop; this chapter will not cluster log data pretending it is counts.
6. Do **not** drop all-zero genes per library. That silently destroys the common-gene set.

Aggregate 10x runs (one matrix, barcode suffixes = libraries) are split by a frozen suffix map, not by guessing.

Ambient RNA: if the GEO file is already a filtered matrix, **do not claim** SoupX, CellBender, or empty-droplet correction.

---

## 3. Per-library QC order (do not reorder)

```text
raw filtered matrix
  → QC metrics (n_genes, total_counts, pct MT, pct HB, ribo, complexity, top50)
  → loose basic filter
  → Scrublet
  → MAD thresholds estimated on singlet candidates only
  → apply MAD + hard caps
  → write clean raw-count h5ad checkpoint
```

### 3.1 Loose filter (Scrublet input)

Starting freeze: `n_genes ≥ 200`, `total_counts ≥ 500`.

Purpose: remove empty droplets that would break Scrublet, **not** to define biology.

### 3.2 Scrublet

- Minimum cells after loose filter: 200. Below that: skip Scrublet, exclude from primary atlas, status `too_few_cells`.
- Expected doublet rate: 0.05 unless the protocol names a per-library rate.
- **Automatic threshold is the call.** Do not force a quantile so that ~X% of cells are removed.
- If predicted rate is < 0.005 or > 0.15, keep the call, write a warning, inspect the score plot. Do not silently retune.
- Exception or crash: status `failed`. Those cells are **not** singlets. Library out of primary atlas.
- High-confidence flags may be stored; the primary deletion rule is `predicted_doublet == True`.

### 3.3 MAD on singlets

Need ≥ 20 singlet candidates. Else fail the library.

On log1p(n_genes) and log1p(total_counts), scaled MAD (`1.4826 × median |x − median|`):

- low bound: median − 3 MAD, then `max` with hard floor
- high bound: median + 4 MAD, then `min` with hard ceiling (genes only)

On %MT and %HB: median + 3 MAD, capped by hard caps.

Starting hard caps:

| Metric | Floor | Ceiling |
|---|---|---|
| genes | 300 | 7500 |
| UMI | 500 | (flag only; do not delete high UMI by default) |
| %MT | — | 15 |
| %HB | — | 3 |

`pass_final_qc` = loose filter ∧ singlet ∧ not low/high genes ∧ not low UMI ∧ not high MT ∧ not high HB.

**Report-only flags** (do not delete unless protocol amended): high UMI, high ribo, low complexity, high top50, %MT > 10 (sensitivity).

### 3.4 Checkpoints

Each clean h5ad stores `pipeline_version`, `config_signature`, `layers["counts"]`, and QC flags.
If the config changes, old checkpoints are invalid. Do not resume them.

---

## 4. Merge into four objects

After all libraries that the *current* manifest still includes:

| Object | Join | Who is in it | Use |
|---|---|---|---|
| all-context union | outer | every successful clean library | paired-context plots |
| primary union | outer | `include_in_primary_atlas` | inventory, gene-presence |
| **primary common** | **inner** | same cells, genes present in **every** primary library | **clustering input** |
| optional stratified downsample | inner genes | ≤ N cells, stratified by dataset×library | parameter scouting only |

Also write `gene_presence_by_dataset.csv`: present-in-any-library vs present-in-all-libraries of that dataset.

**Union outer-join fills missing genes with 0.** Those zeros are *not measured*. Never compute HVG, PCA, or "this gene is off" on the union object.

**Assert `TARGET_GENE` ∈ common genes.** If not, stop. Do not:

- switch to union so the gene "appears"
- pick a family member
- keep going and mention the gene only in a subset of studies as if it were atlas-wide

Optional downsample (e.g. 100k) is for interactive UMAP scouting. All official numbers come from the full common-gene object.

Assertions on the primary objects: unique obs/var names, integer-like counts, all cells `pass_final_qc`, no predicted doublets, `TARGET_GENE` present.

Clinical audit: if symptom/sex/age are `unknown`, print a warning and do not run phenotype contrasts in this chapter.

---

## 5. Global clustering (compute, then stop)

Input: the **primary common-gene** raw-count object, not the downsample, not the union.

Reload `layers["counts"]` into `.X` before normalize. Then:

1. `normalize_total(target_sum=1e4)` + `log1p`
2. Store the result in `layers["normalized"]`; keep `layers["counts"]` untouched
3. HVG: Seurat flavor, `batch_key=dataset`, top 3000. Consensus across batches.
4. PCA to 50 components (arpack). Audit 20/30/40/50. Record the elbow; pick a candidate and freeze it.
5. Harmony on PCA. Keys: **`dataset` required**. Add `sample_id` or `library_id` only if the manifest has technical splits that should not define biology. **Forbidden keys:** `donor_id`, `region`, `condition`, disease labels, cell type, `TARGET_GENE`.
6. Neighbors on Harmony PCA. Grid 10/15/30/50; freeze one (default 30).
7. UMAP `min_dist=0.3` for display. UMAP is not a clustering input.
8. Leiden **grid**, all of: 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.5. Flavor `igraph`, `n_iterations=2`, frozen seed. Write every key as `leiden_r{res}` with `.` → `_`.

Mixing audit: dataset mixing in uncorrected PCA-UMAP vs Harmony-UMAP on a frozen cell subset. Harmony removes study/library offset, not donor biology.

**Checkpoint A — human picks one resolution. The agent does not.**

Use [checkpoint-review-template.md](checkpoint-review-template.md). The reviewer
must cite marker/QC evidence, record a competing interpretation, and declare
whether the target-gene panel was blinded. A resolution cannot be selected from
cluster count, minimum size, visual separation, or target-gene appearance alone.

Export before asking:

- `Leiden_resolution_overview.csv`: resolution, n_clusters, min / median / max cluster size
- one UMAP per resolution (numeric cluster IDs; `legend_loc="on data"` if ≤25 clusters)
- dataset composition long table: cluster × dataset cell counts and within-cluster fractions
- lineage-marker UMAPs / dotplot (tissue panel frozen in the protocol, **not** `TARGET_GENE`)

A software hint is allowed (example: 12–20 clusters **and** min cluster size ≥ 50). It is not a choice. Do not pick the resolution that makes `TARGET_GENE` look cluster-specific.

Then **stop**. Wait for an explicit `selected_leiden_resolution`.

---

## 6. After the human picks a resolution: markers + per-cluster QC

Only now run Wilcoxon, and only on that resolution (optional: the two neighbors on the grid). Expression layer: `layers["normalized"]`. `use_raw=False`. Method `wilcoxon`. `n_genes=50`, `pts=True`, `tie_correct=True`.

### 6.1 Marker tables (required)

| File | Content |
|---|---|
| `{leiden_key}_top50.csv` | long: resolution, cluster, gene, score, pval, pval_adj, logfoldchange, rank, pts if present |
| `{leiden_key}_top20.csv` | same, `rank ≤ 20` |
| `{leiden_key}_top20_by_cluster.csv` | one row per cluster: n_cells, comma-separated top 20 symbols |
| `{leiden_key}_marker_annotation.md` | worksheet with **blank** proposed-label column |

Worksheet columns, in this order:

```text
cluster | n_cells | top20_markers | proposed_label
```

`proposed_label` stays empty until the human fills it. `rank_genes_groups` is not a cell-type call and not a `TARGET_GENE` finding.

Also a compact heatmap of ~5 genes/cluster and a **lineage** (not target-gene) dotplot.

Helper: `python scripts/cluster_review_tables.py clustered.h5ad --leiden-key leiden_r0_5 --out tables/`

### 6.2 Per-cluster QC table (required, this is how bad clusters get caught)

`{leiden_key}_cluster_qc.csv`, one row per cluster, computed from **cell-level QC columns already on the object** (do not re-filter the atlas here):

| Column | Why it is there |
|---|---|
| `cluster`, `n_cells` | size |
| `n_donors`, `n_datasets`, `n_libraries` | a cluster sitting in one study is a batch island until proven otherwise |
| `median_total_counts`, `median_n_genes` | debris vs real cells |
| `median_pct_mt`, `median_pct_hb`, `median_pct_ribo` | dying / RBC / ribo dump |
| `median_pct_top50`, `median_log10_genes_per_umi` | low complexity |
| `median_scrublet_score` | residual doublets that survived the library call |
| `max_dataset_fraction` | technical monopoly |
| `top20_contains_malat1_neat1` | lncRNA/stress island (from the top20 table) |
| `lineage_panel_hits` | how many frozen lineage genes appear in top20 |
| `target_gene_detection_fraction` | **descriptive only**; do not use to delete or name |

Agent **flags**, human **deletes**. Suggested flag rules (write them in the protocol; do not auto-drop):

- top20 dominated by `MALAT1` / `NEAT1` / `KCNQ1OT1` / `XIST` with **zero** lineage-panel hits
- median %MT clearly above other clusters (e.g. > atlas median + 5 points) and no lineage hits
- `n_cells` tiny relative to the grid (e.g. < 50 when the overview promised otherwise)
- `max_dataset_fraction` > 0.90 **and** no lineage hits

A flagged cluster is labeled `low-quality/unresolved` **or** kept as a real state. Either choice is written down. The agent never silently drops it inside a later DEG stage.

Violin+box of %MT, UMI, and n_genes by cluster (cell display, no stars). See [part1-figures.md](part1-figures.md).

**Checkpoint B — human returns a worksheet:** every cluster is either a proposed name or `DELETE: low-quality/unresolved`. No blanks.

Every row must include reviewer ID, date, evidence paths, biological rationale,
alternative interpretation, uncertainty and target-gene blinding status. Missing
fields stop the stage even when the cluster mapping is syntactically complete.

---

## 7. Unresolved cluster → archive and recompute

If any cluster is `DELETE`:

1. Write `objects/00_removed_<leiden_key>_<clusters>.h5ad` plus a barcode CSV (dataset, donor, library, UMI, n_genes, scrublet score, old cluster, reason).
2. Take **retained** barcodes to the immutable QC common-gene **counts** object.
3. Re-run normalize, HVG, PCA, Harmony, neighbors, UMAP, Leiden from scratch.
4. Do not subset the old embedding.
5. Compare v1 vs v2: cluster-count table, UMAP pairwise-distance Spearman / Procrustes, matched-cluster top20 Jaccard. Stability check, not a license to keep old names.
6. Return to **Checkpoint A** (resolution may stay the same number, but clusters are new IDs). Then Checkpoint B again.

If the human keeps a messy cluster, freeze that name (e.g. `Unresolved / stressed`) and do not delete it later by stealth.

---

## 8. Complete naming, save h5ad, lock

Only after Checkpoint B has a **name for every remaining cluster**.

The mapping is a closed dict. Observed Leiden IDs must **equal** the dict keys. Extra or missing keys → abort.

```yaml
# 00_protocol_manifest/cluster_to_celltype.yaml
leiden_key: leiden_r0_5
celltype_key: cell_type_manual_r0_5
display_key: cell_type_manual_r0_5_display
mapping:
  "0": "Type A"    # replace with the human's names; keys must match remaining IDs
  "1": "Type B"
```

Insertion order = legend order.

Write a **new** object. Do not overwrite the unlabeled clustering h5ad.

Columns:

- keep numeric `leiden_key`
- `celltype_key`: pandas Categorical, categories = mapping values in dict order, no NA
- `display_key`: `"{id}: {name}"` (numbers on the UMAP, full names in the legend)

`uns["manual_annotation"]` stores version, keys, mapping, source path.

Save `compression="gzip"`. Re-open `backed="r"` and assert: both columns present, no NA, `nunique` equals len(mapping), `TARGET_GENE` still in `var_names`, `counts` and `normalized` layers still there.

Write `tables/CellType_annotation/{leiden_key}_manual_annotation.csv` (cluster, cell_type, n_cells, fraction) and a short markdown lock report.

Handoff:

```python
assert "counts" in adata.layers and "normalized" in adata.layers
assert adata.obs[celltype_key].isna().sum() == 0
assert set(adata.obs[leiden_key].astype(str)) == set(mapping)
assert TARGET_GENE in adata.var_names
```

Until this object exists, Part 2 does not start. Figures for the lock: [part1-figures.md](part1-figures.md).

---

## 9. Descriptive TARGET_GENE only

Allowed in Part 1:

- atlas-wide detection fraction (`count > 0`)
- mean log-normalized intensity
- the same two numbers per global cell type
- a UMAP of detection vs intensity (display)

Not allowed:

- cluster markers that are just `TARGET_GENE`
- "specific to cell type X" without a donor-unit model (later part)
- stars, q-values, or GSEA

A raw zero is non-detection in this library, not a certified biological off-state.

---

## 10. Directory sketch

```text
analysis/01_qc_merge/
  00_protocol_manifest/   PROTOCOL, manifest.csv, config signature
  per-library clean h5ad
  tables/qc, doublet, gene_mapping, metadata
  reports/run_report.json

analysis/02_global_atlas/
  00_protocol_manifest/
    PROTOCOL.md
    cluster_to_celltype.yaml          # filled only after Checkpoint B
  objects/
    02_leiden_all_resolutions.h5ad    # unlabeled; never overwrite after lock
    03_after_unresolved_removal.h5ad  # if a recompute happened
    04_*_manual_celltypes.h5ad        # locked names; the handoff
    00_removed_*.h5ad                 # archive
  tables/
    Leiden_resolution_overview.csv
    Leiden_markers/{key}_top20.csv
    Leiden_markers/{key}_top20_by_cluster.csv
    Leiden_markers/{key}_cluster_qc.csv
    CellType_annotation/{key}_manual_annotation.csv
  figures/   PNG+PDF+SVG + .parameters.json + source_data/
  reports/   checkpoint notes + lock report
```

QC outputs are never overwritten by clustering. Clustering never writes back into the QC h5ad. The annotated object never overwrites the unlabeled clustering object.

---

## 11. Forbidden sentences (Part 1)

| Do not write | Write |
|---|---|
| We QC'd the concatenated atlas | Each library was QC'd; then merged |
| Scrublet failed so we kept all cells | Library excluded from primary atlas |
| Adjacent tissue is healthy control | Paired context; not in primary atlas |
| Seq2 is another donor | Same donor_id, distinct library_id |
| Gene is off in dataset B (union object) | Gene not present in B's features |
| Clusters are TARGET_GENE-high vs low | Lineage names; TARGET_GENE is descriptive |
| N cells were independently replicated | N libraries / N donors / N datasets, listed |
| Ambient RNA removed | Filtered GEO matrices; ambient not corrected |
| Ready for DEG / CellChat | Labels locked; those are later parts |

