# Part 5 figures — donor-unit association

Shared contracts: [figure-contract.md](figure-contract.md),
[figure-statistics-contract.md](figure-statistics-contract.md), and
[visual-qa-contract.md](visual-qa-contract.md). Every panel carries a
statistics manifest; `NOT_ESTIMABLE` and claim ceilings remain visible.

These panels argue **units, slopes, and holdouts**. They do not argue
cells. UMAP is optional display; it is never the discovery figure.

Theme helper: `scripts/plotting_style.py`. Renderer: `scripts/part5_figures.py`.
Geometry rules that are not restated here: [part1-figures.md](part1-figures.md)
export contract (PNG + PDF + SVG + source CSV + parameter sidecar).

```bash
python scripts/part5_figures.py coverage 02_tables/eligibility.csv --out 03_figures
python scripts/part5_figures.py exposure 03_pseudobulk/metadata.csv --out 03_figures
python scripts/part5_figures.py collinear 03_pseudobulk/metadata.csv --out 03_figures
python scripts/part5_figures.py volcano 02_tables/gene_effects.csv --out 03_figures --gene TARGET_GENE
python scripts/part5_figures.py forest 02_tables/gene_effects.csv --out 03_figures --gene TARGET_GENE
python scripts/part5_figures.py holdout 02_tables/gene_effects.csv --out 03_figures --gene TARGET_GENE
python scripts/part5_figures.py loo 02_tables/loo_with_source.csv --out 03_figures --gene TARGET_GENE
python scripts/part5_figures.py pathways 02_tables/pathway_evidence.csv --out 03_figures
python scripts/part5_figures.py within 02_tables/gene_effects.csv --out 03_figures --gene TARGET_GENE
```

Stem prefix `F05`. One estimand per figure folder. Do not concatenate two
arms onto one volcano and call it a project-level discovery plot.

---

## What is stricter than a Scanpy volcano

| Rule | Why |
|---|---|
| Points are **genes** or **donor-units**, never cells with DE stars | n is not the barcode |
| `NOT_ESTIMABLE` holdout rows stay as **empty slots** with a label | dropping them makes the slope look replicated |
| Two exposure scales are **two columns**, never one mixed axis | collinear scales are not two studies |
| Source colors come from a frozen YAML, not a new husl draw | same source must match Part 4 |
| Caption names `n_units / n_donors / n_datasets / rdf` | q without rdf is a wording hazard |
| Technical genes are gray and unlabelled | they are not strict DEGs |

Do not put a cell UMAP next to a donor-unit forest and imply they share an n.

---

## Required panels

### F05_01 — eligibility / coverage

Stem `F05_01_eligibility`.

Bar or table-as-figure: one row per arm (subtype). Columns the protocol
named: n eligible units, n donors, n datasets, n source_blocks, residual
df, exposure sd, Part 4 flag, Part 5 gate (`formal` / `exploratory` /
`NOT_ESTIMABLE`).

Order arms by **median exposure sd descending**, not by alphabet, unless
the protocol froze another order. Dead arms stay on the figure.

Caption: `n is dataset × donor_id; rdf is of the declared design`.

### F05_02 — exposure range by source

Stem `F05_02_exposure_by_source`.

Unit-level points (one point per eligible unit). x = `source_block` or
`dataset`. y = primary exposure. No cell jitter. Optional box/violin of
**units**. Color by source map.

A source whose sd is below the freeze gets a caption note
`NO_EXPOSURE_RANGE` — not a p-value.

Optional second panel: the alternative scale. Same x order. Do not plot
both scales on one y-axis.

### F05_03 — collinear scales (required if both exposures were fit)

Stem `F05_03_scale_collinearity`.

One point per eligible unit: Jeffreys detection vs `log2(CPM+1)`. Annotate
Spearman. Title: **these two scales may be collinear**. `|ρ| ≥ 0.9` is a
note, not a second cohort.

Skip only if the protocol froze a single exposure and forbade the second.

### F05_04 — gene volcano (donor-unit)

Stem `F05_04_{GENE}_volcano` where `{GENE}` is `TARGET_GENE` (exposure),
not a highlighted outcome.

x = logFC per exposure unit (state the unit in the axis label:
`per +10 pp Jeffreys` or `per 1 SD log2CPM`). y = −log10 p (or q; sidecar
must name which).

Color: `robust_primary` / not / technical. Do not label more than 20
genes. Do not star cells. Do not draw `TARGET_GENE` as an outcome point
(it was excluded).

If the arm is `NOT_ESTIMABLE`, do not draw a volcano of a failed fit.
Draw F05_01 / F05_08 instead.

### F05_05 — TARGET_GENE forest

Stem `F05_05_TARGET_GENE_forest`.

Pre-declared panel only (GR support, negative controls, or the protocol's
list). Horizontal forest: logFC + 95% CI. Order by the protocol list, not
by p. Dual method (limma vs edgeR) may be two symbols; they must share the
x-axis scale.

Genes absent from the filter are blank slots, not zeros.

### F05_06 — pathway concordance

Stem `F05_06_pathway_camera_fgsea`.

CAMERA q vs fgsea q for the frozen libraries, or a two-column bar of
signed NES / effect for pathways that meet the dual-method rule. Technical
leading-edge fraction must be in the source CSV.

Empty CAMERA after a large shift is a valid figure (zero `robust_primary`
pathways). Do not replace it with a prettier ORA.

### F05_07 — donor LOO

Stem `F05_07_donor_loo`.

One point per omitted donor: effect (or signed concordance) for the
pre-declared primary outcome gene, colored by `source_block`. Caption must
say **donor LOO**, not LODO.

If `LOO_SOURCE_UNBALANCED`, the caption states the fraction of iterations
that omitted the dominant block.

### F05_08 — source-block holdout (hero robustness)

Stem `F05_08_source_holdout_forest`.

Rows: `FULL`, `DROP_EXT`, `DROP_DOMINANT`, optional split-half. Two
columns if two exposures were fit — **never mixed on one axis**.
`NOT_ESTIMABLE` rows are drawn as empty positions with that label. Do not
omit them.

Caption: `holdout of a source block, not of one donor`.

This is the panel that forbids "independently replicated".

### F05_09 — within-source slopes

Stem `F05_09_within_source`.

One forest per `dataset` or `source_block`. Same exposure scale as F05_08
left column. A source with `NO_EXPOSURE_RANGE` is an empty slot.

Caption: within-source agreement is not cross-source replication.

### F05_10 — optional display UMAP

Stem `F05_10_display_umap`.

Part 3 embedding, subtype colors, **no DE stars**. Allowed only as
display of which population was modeled. Caption: `display; n is not cells`.

Skip unless the protocol asks. Never the lead figure.

---

## Caption ceilings

Allowed: "associated at donor-unit level", "per +10 pp", "rdf",
"source-block holdout", "NOT_ESTIMABLE", "same sign in dual method",
"collinear scales", "within-source consistent".

Forbidden: "drives", "regulates", "independently replicated",
"cross-dataset validated", "LOO N/N therefore general", "marker of",
"cell-level DEG", "Part 4 already showed association".

Every formal panel's sidecar lists input SHA-256 and the exposure scale.

---

## What not to draw as a finding

- Cell Wilcoxon volcano
- Detected vs undetected cell UMAP as the association
- Dropping `NOT_ESTIMABLE` rows from the holdout forest
- Mixing Jeffreys and log2CPM on one forest axis
- A new husl palette for sources already in Part 4 YAML
- GSEA enrichment plots from a cell-rank file
- High-tail OR as the lead
- Global atlas UMAP colored by Part 5 q-values


