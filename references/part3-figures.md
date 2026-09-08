# Part 3 figures — rounds first, names last

This catalog is **stricter** than Part 1 Checkpoint B. Formal figures are
matplotlib through `scripts/plotting_style.py`, not `sc.pl.*`. Renderer:
`scripts/part3_figures.py`.

Theme, point size, PNG 300 / PDF fonttype 42 / SVG + sidecar + source CSV:
same contract as [part1-figures.md](part1-figures.md). Do not invent a
second style. Reuse `scripts/plotting_config.yaml`.

During cleanup rounds the hero is **coverage + QC + contamination**, not a
pretty named UMAP. The named UMAP exists only after STOP POINT 3.

```bash
python scripts/part3_figures.py stop1 clustered.h5ad \
    --out figures/round_00_initial --gene TARGET_GENE
python scripts/part3_figures.py stop2 clustered.h5ad \
    --leiden-key leiden_r0_5 --panels config/marker_panels.yaml \
    --out figures/round_00_initial
python scripts/part3_figures.py exclusion parent_clustered.h5ad \
    --partition tables/round_01_after_removal/parent_round_partition_manifest.csv \
    --out figures/round_00_initial
python scripts/part3_figures.py membership \
    tables/round_01_after_removal/old_to_new_membership.csv \
    --out figures/round_01_after_removal
python scripts/part3_figures.py lock annotated.h5ad \
    --group subtype --id-key leiden_r0_5 \
    --colors config/subtype_colors.yaml \
    --out figures/final
```

---

## What is stricter than Part 1

| Rule | Part 1 Checkpoint B | Part 3 |
|---|---|---|
| QC violin x-order | cluster ID, same for every metric | **each panel independently** ordered by **that metric's median, descending** |
| Strip / jitter on QC | not required | **forbidden** (no cell sampling on formal QC) |
| Numeric cluster colors | husl / Scanpy default illegal | frozen `theme.cluster_cycle` in `plotting_config.yaml` |
| Color by cluster ID | — | color follows the **ID**, not the panel's x-order, so cluster 0 stays the same hex when a QC panel reorders |
| Cluster 0 across rounds | — | **not** a biological identity; IDs die with HVG recompute |
| Names on UMAP | after global lock | **only** at STOP POINT 3; rounds show numbers |
| Marker display | top20 as written | **strict-positive** filter; `TARGET_GENE` dropped from the worksheet |
| Harmony covariates | `dataset` required; sample allowed | **`dataset` only** |
| `sc.pl.*` as the panel | banned | banned |
| p-value stars | banned | banned |
| Rotating UMAP | banned | banned |
| Sampling cells to "thin" a UMAP | banned | banned; removed cells live in the archive, not a downsample |

Caption titles: left-aligned, bold. Median stroke: `#B2182B`. Unknown:
`#808080`. Retained / removed: `#0072B2` / `#D55E00`.

---

## Ordering contract

Comparable quantities **descend**.

- Cluster-size bar: cell count descending.
- UMAP legend (rounds): cell count descending; draw smaller groups last.
- Each QC violin: **its own** median descending. Never reuse another panel's
  order. Sidecar must record `panel_cluster_orders`.
- Marker-panel columns: maximum displayed mean log-normalized expression
  descending inside each panel block.
- Marker-panel rows: cluster cell count descending.
- PC index, Leiden resolution, and neighbor-size axes keep their **intrinsic**
  numeric sequence. They are outside the descending-quantity rule.
- UMAP coordinates are never rotated, flipped, or re-fit to make a cluster
  central.

---

## Per-round folder

```text
figures/<round_name>/
  01_HVG/
  02_PCA/
  03_Harmony/
  04_UMAP/
  05_metadata/
  06_Leiden/
  07_markers/          # STOP POINT 2 only
  08_contamination/    # STOP POINT 2 only
  09_QC/               # STOP POINT 2 only
  10_removal/          # after prepare-removal
  source_data/
```

Every panel: PNG 300 dpi, PDF fonttype 42, SVG unless the sidecar says the
scatter exploded, `source_data/<stem>.csv`, `.parameters.json`. Rasterize
UMAP points (`rasterized=True`).

---

## Evidence levels

| Level | When | Job | Hero |
|---|---|---|---|
| **S1** | STOP POINT 1 | Is this resolution usable? | coverage heatmap + resolution landscape |
| **S2** | STOP POINT 2 | KEEP or DELETE? | QC violins + contamination panel + selected UMAP |
| **R** | after DELETE | What left, and did IDs change? | exclusion UMAP + membership heatmap |
| **L** | STOP POINT 3 only | Locked names | numbered UMAP, names in the legend |

Do not lead a cleanup-round report with a named UMAP. There are no names yet.

---

## STOP POINT 1 — required

Stem prefix `F03_S1_`. One round, one embedding. Do not reuse the global atlas
UMAP.

**F03_S1_01 — HVG mean vs dispersion**

- x = gene mean, y = normalized variance / dispersion from the Seurat HVG
  table. Highly variable vs the rest.
- Mark `TARGET_GENE` with an open marker if it is still in `var_names`.
- Title must say it was **excluded from PCA input**. If the gene is inside
  `highly_variable`, that is a failed assert, not a figure.

**F03_S1_02 — PCA variance**

- Bars = `uns["pca"]["variance_ratio"]` after the gene was dropped.
- Intrinsic PC order. Annotate the frozen number of PCs used for Harmony.

**F03_S1_03 — Harmony before / after, colored by `dataset`**

- Left = `obsm["X_pca_unintegrated"]` PC1–PC2. Right = Harmony-corrected
  `X_pca` PC1–PC2. Same point size rule as UMAP.
- If unintegrated PCs were not stashed, **do not fake a before panel**.
  Draw after only and write that failure in the sidecar.

**F03_S1_04 — UMAP by dataset**

- Categorical. Legend if n_datasets ≤ 12; otherwise drop the legend and
  say so.

**F03_S1_05 — metadata UMAPs**

- Required separate panels: `dataset`, `donor_id`, `library_id` (or
  `sample_id` if that is the library key), log10 UMI, n_genes, Scrublet
  score.
- Donor legend: omit if n_donors > 12 (color is still drawn).
- Continuous: viridis; sort points by value so high cells paint last.
- Optional diagnostic UMAP of `TARGET_GENE` normalized expression — title
  must say **diagnostic**, not a cluster definition. The Part 2 catalog is
  not this checkpoint.

**F03_S1_06 — one UMAP per Leiden resolution**

- Numeric IDs **on the data** (median coordinate, bold, white boxed face).
- Color = `cluster_cycle[id % len]`. Legend = `"{id} (n=…)"`, cell-count
  descending.
- No subtype names. No `TARGET_GENE` coloring on this panel.

**F03_S1_07 — resolution landscape**

- From `resolution_summary.csv` (or computed from `leiden_*` columns):
  n_clusters, min size, n donor-units, max dataset fraction vs resolution.
- Resolution axis stays numeric.

**F03_S1_08 — graph sensitivity**

- Heatmap of n_clusters (or a stated metric) at PCs `{20,30,40}` ×
  neighbors `{15,30,50}` for the **candidate** resolution.
- The figure script **plots** a CSV written during `run-leiden`. It does
  not recompute neighbor graphs.

**F03_S1_09 — coverage heatmap**

- Rows = `dataset × donor_id`. Columns = Leiden clusters of **one**
  candidate resolution, or a small grid if the human has not picked yet
  (default: plot the protocol's working candidate, not every resolution).
- Values = cell counts (display may be log1p; sidecar says which).
- This is the S1 hero: a cluster that is one dataset, or one donor, is
  visible here before anyone names it.

---

## STOP POINT 2 — required

Stem prefix `F03_S2_`. Selected resolution only.

**F03_S2_01 — selected-resolution UMAP, numbers only**

- Same construction as F03_S1_06. Title names the Leiden key, not a biology.

**F03_S2_02 — cluster-size bar**

- Descending n_cells. Annotate n_donors on the bar or to the right.
- Color by cluster ID via `cluster_cycle`.

**F03_S2_03 — QC violin + box**

- Required metrics, each its own panel: `%MT`, `log10(UMI+1)`, n_genes.
- Optional extra: Scrublet, `%ribo`, complexity, `%HB` — still own-median
  order.
- Violin α ≈ 0.45, box α ≈ 0.18, **no fliers**, whiskers 5th–95th, median
  `#B2182B`.
- **No strip / swarm / sampled points.** Stage-style "100 random cells"
  is not allowed on this formal panel.
- x-tick labels = cluster IDs in **that panel's** median-descending order.
- Color of the violin body = `cluster_cycle` of that ID, not the x position.

**F03_S2_04 — coverage heatmap** of the selected key (same rules as S1_09).

**F03_S2_05 — lineage-program marker matrix**

- Frozen YAML panels supplied by the protocol. Dot area = raw-count
  detection fraction. Color = mean `normalized`.
- Rows = clusters by cell count descending. Scanpy one-gene or default
  `sc.pl.dotplot` is banned.
- Drop `TARGET_GENE` if it sneaks into the YAML; warn in the sidecar.

**F03_S2_06 — contamination / off-lineage marker matrix**

- Same geometry as 05, different panel block. This is what DELETE
  decisions are allowed to cite.

Do not put Wilcoxon p-values, stars, or "marker of" on these panels.

---

## After DELETE — required

Stem prefix `F03_R_`.

**F03_R_01 — exclusion UMAP**

- Coordinates = **parent** round UMAP. Color = retained vs removed
  (`#0072B2` / `#D55E00`). Draw retained first, removed last.
- Caption: parent-round geometry; the child manifold is recomputed.
  Do not call this the child UMAP.

**F03_R_02 — old → new membership heatmap**

- Rows = parent Leiden IDs. Columns = child Leiden IDs. Values = shared
  barcodes. Annotate integers.
- Title must say new IDs are **not** the old IDs. Empty child columns
  after a recompute are allowed; do not force a 1–1 mapping.

**F03_R_03 — partition counts**

- Bar or table-as-figure: n cells / n units removed vs retained, by
  deleted cluster and by dataset. Source is the partition manifest.

Removed cells are **not** plotted as a faint ghost on the child UMAP.

---

## STOP POINT 3 / lock — required

Stem prefix `F03_L_`. Illegal until the decision file has zero `DELETE`
(or an explicit `NAME_NOW`).

**F03_L_01 — locked subtype UMAP**

- Color = frozen YAML. **Numbers on the data**, full names **only** in
  the legend. Legend order = mapping insertion order (or cell-count
  descending if the sidecar says so). Draw smaller groups last.
- This is the Part 3 lock hero.

**F03_L_02 — optional marker-pair / lineage matrix** on the named key.

Same export contract. Unknown extra category → `#808080`, never a new
random husl.

---

## Marker-panel YAML

```yaml
lineage:
  title: Compartment lineage programs
  genes:
    panel_a: [GENE1, GENE2]
    panel_b: [GENE3]
contamination:
  title: Off-lineage contamination
  genes:
    platelet: [PPBP, PF4]
    erythrocyte: [HBB, HBA1]
```

Do not put `TARGET_GENE` here. Do not put a gene that was used as a Harmony
key (there are none) or as a cluster name.

---

## Caption ceilings

Allowed: "numeric Leiden ID", "retained", "removed", "dataset × donor
coverage", "median %MT", "excluded from HVG/PCA", "diagnostic expression".

Forbidden on the figure: "marker of", "specific to", "upregulated",
"significantly higher", "replicated", "defines a subtype", "tissue
boundary", "old cluster 3 is new cluster 3", a `TARGET_GENE`+ name.

UMAP empty space is unsampled, not a missing intermediate state.

---

## What not to draw as a finding

- `sc.pl.umap` / `sc.pl.dotplot` / `sc.pl.violin` as the formal panel
- Wilcoxon stars on cells
- A pie of cluster abundance
- Jet / rainbow / 3-D bars
- The global atlas UMAP colored by the new compartment Leiden
- A child UMAP that is just the parent embedding with points deleted
- Names on an intermediate round
- Sampling cells so a QC violin looks less overplotted
- Coloring clusters by `TARGET_GENE` and calling that the subtype map

