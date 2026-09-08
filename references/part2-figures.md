# Part 2 figures — what an agent may draw for TARGET_GENE

This is the catalog. Part 2 exists so the agent does **not** dump
`sc.pl.umap(gene)` and call it a survey.

Theme, point size, PNG/PDF/SVG + sidecar: same contract as
[part1-figures.md](part1-figures.md). Do not invent a second style.
Helpers: `scripts/plotting_style.py`, `scripts/plotting_config.yaml`.
Renderer: `scripts/part2_figures.py`.

Part 1 already locked a numbered cell-type UMAP. Part 2 reuses those
coordinates and colors. Numbers-on-data are optional here; the legend
carries names.

Titles on quantitative panels: **left-aligned, bold**. Median stroke in
boxplots: `#B2182B`. Unknown category: `#808080`, never a random husl.

```bash
python scripts/part2_figures.py atlas.h5ad \
    --gene TARGET_GENE \
    --group cell_type \
    --colors config/celltype_colors.yaml \
    --out figures/part2 \
    --min-cells 20
```

---

## Evidence levels (draw in this order)

| Level | Job | Panels |
|---|---|---|
| **A** Localization | Where on the locked UMAP | F02_01–05 |
| **B** Decomposition | Prevalence vs intensity vs UMI share | F02_06–10 |
| **C** Donor-unit | One point = `dataset × donor` × type | F02_11–12 |
| **D** Robustness | Depth, cell number, hold-out dataset | F02_13–14 |
| **E** Optional | High-expression tail; true paired units | only if the protocol names them |

Do not lead the report with a cell-level violin. The **hero** is F02_11
(the four donor-unit endpoints). A and B are context.

---

## Required catalog

Stem = `F02_{nn}_{GENE}_{slug}`. Every panel: PNG 300 dpi, PDF fonttype 42,
SVG unless the sidecar says the scatter exploded, `source_data/*.csv`,
`.parameters.json`.

Rasterize UMAP points (`rasterized=True`). Do not rotate the embedding.

### A. Localization

**F02_01 — locked cell-type UMAP**

- Color = Part 1 YAML. Legend = names, ordered by **descending cell count**.
- Draw smaller groups last so they stay visible.
- No new labels, no new colors. This is orientation, not a finding.

**F02_02 — feature UMAP, full scale**

- Layer: `normalized`. Cmap: `magma`. vmin = 0.
- **Keep zeros.** Sort points by value so high cells paint last.
- Colorbar: `Log-normalized {GENE}`.

**F02_03 — feature UMAP, display-only p99**

- Same as 02, but clip the **display** at the 99th percentile of
  **positive** cells. Underlying values stay untouched.
- Title must say the clip is display-only.

**F02_04 — detection UMAP**

- Layer: `counts`. Detected (count > 0) = `#D73027`. Undetected = `#D9D9D9`.
- Draw undetected first. Legend reports n.
- Caption: non-detection includes dropout.

**F02_05 — pair: cell types | intensity (same coordinates, same axes)**

- Left = F02_01 without the fat legend if it collides; right = F02_03.
- This is the money localization figure. Do not compute a second UMAP
  for the right panel.

### B. Decomposition (cell types)

One-gene **Scanpy dotplots are banned** (a single-column matrix is empty
space). Use a lollipop.

**F02_06 — lollipop**

- y = cell type, ordered by **ascending all-cell mean** (low at bottom).
- x = all-cell mean log-normalized (zeros in).
- Dot area = detection fraction. Stem is a light gray hline from 0 to x.
- Annotate `n=` donor count on the right.
- Color dots from the YAML.

**F02_07 — violin + box, all cells**

- Horizontal. Violin α ≈ 0.42, box α ≈ 0.18, **no fliers**.
- Whiskers = 5th–95th percentile. Median = `#B2182B`.
- **Order this panel by its own median, descending.** Never reuse another
  panel's order.
- Annotate cell counts. No significance stars (cells are not n).

**F02_08 — violin + box, detected cells only**

- Same geometry as 07, subset count > 0. Own median order.
- Together with 07, this stops "the violin looks high because everyone
  is zero except a spike".

**F02_09 — detection vs positive-cell intensity**

- One point per cell type. x = detection fraction, y = mean normalized
  among detected cells. Size ∝ n_cells.
- Title states the point: **these two ranks are allowed to disagree.**
- Direct labels; no legend of 13 overlapping colors without names.

**F02_10 — composition vs count contribution**

- Paired horizontal bars per type: gray = cell fraction of the atlas;
  YAML color = share of all `TARGET_GENE` UMIs.
- Order by **descending UMI contribution**.
- Caption: a large type can dominate UMIs without leading per-cell intensity.

### C. Donor-unit (hero)

**F02_11 — four-endpoint box + strip grid**

2×2, one **point per eligible `dataset × donor` × type** (≥20 cells):

| Panel | Endpoint |
|---|---|
| a | detection fraction |
| b | positive-cell intensity |
| c | all-cell mean (zeros in) |
| d | log2(CPM + 1) on summed counts |

- Horizontal box (no fliers) + jittered points, YAML fill, red median.
- **Each panel ordered by its own median.** n= on the right is donor-units,
  not cells.
- Supertitle states the eligibility threshold.

This is the figure that decides whether a Part 3 lineage is even on.

**F02_12 — unit × type heatmap of detection**

- Rows = `dataset_donor_id`, columns = types, cmap magma, vmin 0.
- Missing type in a unit = NaN / blank, not a biological zero.
- Do not impute.

### D. Robustness (descriptive)

**F02_13 — leave-one-dataset-out heatmap**

- Rows = held-out dataset. Values = remaining-data detection (or rank).
- Caption **must** say this is not independent replication.
- Skip if fewer than two datasets.

**F02_14 — depth diagnostics**

- Scatter: unit detection vs log1p(mean UMI); vs n_cells. Color = type.
- If the cloud is a line, the survey is a depth story and Part 5 will
  need the depth covariates that already exist as defaults.

### E. Optional (protocol must name them)

**High-expression tail forest** — enrichment of cells above a frozen
percentile (e.g. positive-cell top 10%) by type. Odds ratios are
**concentration**, not a new cluster and not a mechanism. If drawn,
the caption repeats that sentence.

**True paired units** — connect the same `dataset × donor` across two
tissues/conditions only when both sides exist. Missing pairs are skipped
and counted, never invented. Adjacent tissue is not a healthy control
unless the Part 1 manifest said so.

**Balanced subsampling** — for a sparse gene, draw a fixed n cells per
eligible unit × type and show that detection still ranks the same.
Useful; not a p-value.

Do not add volcano plots, enrichment bars, CellChat, or a new embedding.

---

## Ordering rules (non-negotiable)

| Plot | Order by |
|---|---|
| UMAP legend | cell count descending |
| Lollipop y | all-cell mean ascending |
| Each violin | median of **that** distribution descending |
| Contribution bars | UMI contribution descending |
| Each donor-unit panel | median of **that** endpoint descending |
| Intrinsic axes (PC, dose, time) | leave in native order |

Ties: stable sort. Record the rule in the sidecar.

---

## Layers

| Need | Layer |
|---|---|
| Detection, UMI sums, CPM, contribution | `counts` |
| Feature UMAP, violin, positive intensity, lollipop x | `normalized` |
| Coordinates | Part 1 `X_umap` only |

If `.X` is what you have, stop and find the layers. Do not guess.

---

## Caption ceilings (copy onto the figure or the report)

Allowed: "detected in", "highest fraction", "largest UMI share",
"donor-unit median", "ranks disagree", "tracks library depth".

Forbidden on the figure: "marker of", "specific to", "upregulated",
"significantly higher", "replicated across studies", "defines a state",
"tissue boundary".

---

## What not to draw as a finding

- `sc.pl.umap` / `sc.pl.dotplot` / `sc.pl.violin` defaults as the formal panel
- Wilcoxon stars on cells
- A pie of cell-type abundance
- Jet / rainbow / 3-D bars
- A second UMAP computed on a subset still labeled "global atlas"
- Hiding zeros so the gene looks uniformly on
- Rotating UMAP to put the high type in the center

---

## Part 4 reuse

Specified: [part4-subtype-survey.md](part4-subtype-survey.md),
[part4-figures.md](part4-figures.md).

Part 4 is this catalog with prefix `F04`, a **locked subtype** key, and
the Part 3 embedding/colors. Same geometry. Do not rewrite the aesthetics.
The extra panel is identifiability (`F04_15`), not a new violin.

