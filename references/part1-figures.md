# Part 1 figures

Shared contracts: [figure-contract.md](figure-contract.md),
[figure-statistics-contract.md](figure-statistics-contract.md), and
[visual-qa-contract.md](visual-qa-contract.md). Part 1 lock figures must record
their data unit, claim ceiling and final-size visual QA.

Rules distilled from a production atlas plotting stack (`plotting_style.save_figure`,
publication rcParams, UMAP point scaling, annotated UMAP with numeric labels).
Do not invent a second theme inside Part 1.

Formal figures are **not** Scanpy's default `sc.pl.*` PNG dump. Build a matplotlib
figure, then save through one helper.

## Export contract (every formal figure)

```text
figures/<panel_stem>.png          # 300 dpi
figures/<panel_stem>.pdf          # vector; pdf.fonttype = 42
figures/<panel_stem>.svg          # unless the plot is raster-heavy and SVG explodes
figures/<panel_stem>.parameters.json
figures/source_data/<panel_stem>.csv
```

Sidecar JSON includes: UTC time, python/matplotlib versions, output paths, the
plot parameters, and SHA-256 of input h5ad/tables.

`bbox_inches=tight`, `pad_inches=0.04`, white background, not transparent.

If SVG is skipped, say so in the sidecar. Rasterize UMAP points in PDF/SVG
(`rasterized=True` on the scatter) so files stay usable.

## Theme

```text
font: Arial, then Helvetica, then DejaVu Sans
CJK fallback if labels need it: Microsoft YaHei / SimHei (do not download fonts)
base 8 pt, title 9.5, axis 8.5, ticks 7.5, legend 7
axes.linewidth 0.8
spines: top and right off
foreground #202124, grid #D9D9D9 alpha 0.35 (UMAP: no grid)
continuous: viridis (metadata) / magma (expression)
divergent: RdBu_r
never jet / rainbow / 3-D bars / pies / decorative gradients
```

Named sizes (inches): single 3.45×3.15, double 7.10×3.20, wide 8.60×4.00,
multi 9.00×7.20. Annotated global UMAP may be wider (e.g. 14×8) so the legend fits.

## UMAP point size

Deterministic from cell count, not from "what looks pretty":

| n_cells | size | alpha |
|---|---:|---:|
| < 5,000 | 8.0 | 0.85 |
| < 20,000 | 3.5 | 0.72 |
| < 75,000 | 1.25 | 0.58 |
| otherwise | 0.65 | 0.45 |

Do not rotate or flip UMAP to manufacture a narrative.

## Checkpoint A (resolution grid)

Required:

- one UMAP per Leiden resolution, cluster IDs on data
- Harmony before/after for `dataset` (and `sample_id` if it was a Harmony key)
- metadata panel: dataset, donor_id, library/sample, log UMI, n_genes, scrublet score
- lineage-marker UMAP or dotplot (protocol panel)
- optional diagnostic UMAP of `TARGET_GENE` normalized expression — title must say
  diagnostic, not a cluster definition. The **full** gene catalog (detection vs
  intensity vs donor-unit grid) is Part 2, not this checkpoint.

## Checkpoint B (cluster QC)

Required:

- violin **plus** box of `pct_counts_mt`, `total_counts`, `n_genes_by_counts` by cluster
- order x-axis by cluster ID (intrinsic), not by the metric
- no significance stars (cells are not n)
- bar or stacked bar of dataset composition per cluster (source: the composition table)

## Lock figure (after complete naming)

One annotated UMAP:

- color = `display_key` (`"{id}: {name}"`)
- **numbers drawn at the median UMAP coordinate of each cluster**, bold, with a
  white circle / boxed face so they stay readable
- full names only in the legend, right margin
- after lock, colors come from a frozen YAML `cell_type: {name: "#hex"}`
- unknown / extra category → visible gray; do not sample a random husl palette
  for a formal figure
- legend order = mapping insertion order (or cell-count descending if the sidecar
  says so). Draw smaller groups last so they remain visible

Do not claim embedding contours are anatomical tissue boundaries.

## Color YAML

Write `config/celltype_colors.yaml` when names are locked. Keys are the **exact**
AnnData category strings. Colorblind-safe discrete palettes (Okabe–Ito or the
project's existing hex list). Never generate a new random color on re-plot.

## What not to plot as a finding

- Wilcoxon p-values on cells
- `TARGET_GENE` as the only color that "defines" a cluster
- pie charts of cluster abundance
- a second UMAP computed on a subset while still calling it the global atlas

