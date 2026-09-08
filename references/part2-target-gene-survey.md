# Part 2 — Global TARGET_GENE survey

Audience: an agent with a **Part 1 locked atlas**. Labels are frozen. UMAP is
frozen. This chapter only asks where `TARGET_GENE` is detected and how that
description looks at `dataset × donor_id`.

It is not DEG. It is not a new clustering. It is not "the gene marks type X".

Operating figures: [part2-figures.md](part2-figures.md).
Theme / export helper: `scripts/plotting_style.py`.
Catalog renderer: `scripts/part2_figures.py`.

---

## 0. What this part answers

In the locked global types:

1. Which types even contain cells with `TARGET_GENE` count > 0?
2. Is that mostly **prevalence** (many cells weakly on), **intensity** (few
   cells strongly on), or **count contribution** (a large type donates most UMIs)?
3. Do those four donor-unit endpoints still move together when one dataset is
   held out, and are they obvious functions of depth?

It does **not** answer whether the gene is DE, a driver, a marker used for
annotation, or independently replicated across studies.

Claim type: descriptive survey. Evidence ceiling: descriptive.
Part 3 is allowed to start after this; Part 5 is not.

---

## 1. Inputs (read-only)

- The Part 1 annotated h5ad. Hash it. Do not write into it.
- Frozen color YAML from Part 1 (`cell_type: {exact AnnData string: "#hex"}`).
- `layers["counts"]`, `layers["normalized"]`, `obsm["X_umap"]`.
- `obs` keys: the locked label column, `dataset`, `donor_id` (and
  `library_id` if present). If `dataset_donor_id` is absent, build it **in
  memory** from `dataset` + `donor_id`. Do not add columns to the authoritative
  object.

Assert `TARGET_GENE` is in `var_names`. If not, stop. Do not plot a homolog.

Do not recompute HVG, PCA, Harmony, neighbors, or UMAP.

---

## 2. Definitions (freeze these)

| Name | Source | Formula |
|---|---|---|
| Detected cell | `counts` | gene UMI > 0 |
| Detection fraction | counts | n_detected / n_cells |
| Positive-cell intensity | `normalized` | mean(log1p CP10k among detected cells) |
| All-cell mean | `normalized` | mean including zeros |
| Count contribution | `counts` | sum(gene UMI in type) / sum(gene UMI in atlas) |
| Pseudobulk CPM | `counts` | 1e6 × sum(gene UMI) / sum(total UMI) in a unit × type |
| Biological unit | metadata | `dataset × donor_id` |
| Eligible unit × type | counts | n_cells ≥ 20 (sensitivity: 5 / 10 / 50) |

These four endpoints are **not interchangeable**. A type can lead detection
and lose intensity, or lead cell fraction and lose UMI contribution.

A raw zero is non-detection in this library, not a biological off-state.
Do not hide zeros.

---

## 3. What to compute (tables before plots)

Write CSVs first. Plots read tables (except UMAPs, which read the embedding).

1. Cell-type summary (one row per locked type): n_cells, n_donors, n_datasets,
   detection, positive intensity, all-cell mean, count contribution.
2. Unit × type table: the same endpoints plus mean UMI, eligibility flag.
3. Optional LODO table: drop one `dataset` at a time, recompute type-level
   detection on remaining units (cell-weighted). This is a **description**,
   not a replication claim.

Do not put Wilcoxon / Friedman p-values on a figure. If you compute a
cell-level test, park it in `tables/exploratory/` and do not quote it in the
report as a finding.

---

## 4. What Part 2 may claim

- Type A has the highest detection fraction (report the number and n_cells).
- Type B contributes the largest share of `TARGET_GENE` UMIs.
- Among eligible donor-units, type C has the highest median detection.
- Detection and positive-cell intensity do not rank types the same way.
- After holding out dataset D, the remaining ranking still / no longer
  puts type A first (**descriptive**).

## 5. What Part 2 may not claim

- `TARGET_GENE` marks / defines type X
- `TARGET_GENE` is differentially expressed
- Adjacent / proximal tissue is a healthy control
- N cells = N replicates; violin stars; Scanpy `rank_genes_groups` as discovery
- Leave-one-dataset-out = independent replication (source-block is Part 5)
- High-expression tail OR = a new cell type or a mechanism
- Embedding hulls are anatomical borders

---

## 6. Lifecycle

```text
Part 2 progress:
- [ ] Confirm Part 1 lock (label key, color YAML, counts/normalized/UMAP)
- [ ] Hash the input h5ad; open read-only
- [ ] Assert TARGET_GENE in var_names
- [ ] Build group summary + unit × type table (min 20 cells)
- [ ] Render the figure catalog (part2-figures.md)
- [ ] One-page survey: which types are even on, which endpoints disagree
- [ ] Recommend which lineage is worth a Part 3 recluster — or none
```

Default split: expensive model writes the one-page survey; cheap model
runs `part2_figures.py`.

