# Part 4 figures — Part 2 catalog on locked subtypes

Do not rewrite the aesthetics. Panels F04_01–14 are F02_01–14 with a new
prefix, the **Part 3** UMAP, and the **Part 3** color YAML.

Full geometry, ordering, layers, caption ceilings:
[part2-figures.md](part2-figures.md). Theme helper unchanged.

The extra panel is identifiability, not a prettier violin.

```bash
python scripts/part4_figures.py locked_subtype.h5ad \
    --gene TARGET_GENE \
    --group subtype \
    --colors config/subtype_colors.yaml \
    --out figures/part4 \
    --min-cells 20
```

Equivalent catalog-only call:

```bash
python scripts/part2_figures.py locked_subtype.h5ad \
    --gene TARGET_GENE --group subtype \
    --colors config/subtype_colors.yaml --out figures/part4 \
    --stem-prefix F04 --umap-title "Locked subtypes"
```

---

## What changes vs Part 2

| Piece | Part 2 | Part 4 |
|---|---|---|
| Object | Part 1 locked atlas | Part 3 locked compartment |
| Group key | global cell type | frozen subtype **name** |
| Colors | Part 1 YAML | Part 3 YAML |
| UMAP | global | this compartment's `X_umap` |
| Stem prefix | `F02` | `F04` |
| Hero | F02_11 four-endpoint grid | F04_11, plus F04_15 as the Part 5 gate |
| Extra | — | identifiability table / heatmap |

Numbers-on-data remain optional. Legend carries names. Unknown extra
category → `#808080`. Do not sample a new husl palette if a subtype is
missing from the YAML — use gray and list it in the sidecar.

Do not draw the global atlas UMAP and color it by the new subtypes.
Do not rotate the embedding.

F04_11 still decides which subtypes are even on. F04_15 decides which
slices can go to Part 5.

---

## Required extra — F04_15 identifiability

Stem `F04_15_{GENE}_identifiability`.

- One cell per **subtype × source** (`source_block` if the protocol has
  it, else `dataset`).
- Display: eligible-unit detection **sd**, or a categorical flag map.
  Sidecar must name which.
- Missing / n < 2 eligible units → blank or gray, not a zero sd.
- Caption: `sd ≈ 0 cannot identify a Part 5 slope`. No p-values.
- Source CSV is the real deliverable: n eligible units, detection
  sd/IQR/min/max, log2(CPM+1) sd/IQR, Spearman of the two scales, flag.

Flag colors (frozen):

| Flag | Hex |
|---|---|
| `RANGE_OK` | `#009E73` |
| `EXPLORATORY_N` | `#E69F00` |
| `WITHIN_SOURCE_RANGE` | `#56B4E9` |
| `NO_EXPOSURE_RANGE` | `#D55E00` |
| `LOW_N` | `#808080` |

Optional **F04_16** — unit-level detection vs log2(CPM+1), one point per
eligible unit, YAML color. Title: **these two scales may be collinear**.
Skip unless the protocol asks. `|ρ| ≥ 0.9` is a note, not a second cohort.

Optional high-tail (Part 2 Level E) stays optional. Caption may not say
the tail is a subtype.

---

## Caption ceilings (additions)

Allowed: "detected in", "highest fraction", "donor-unit median",
"ranks disagree", "no exposure range in source S", "collinear scales".

Forbidden: "marker of", "specific to", "upregulated", "associated with",
"replicated", "defines a state", "Part 5 will pass", "independently
identifiable across studies".

---

## What not to draw as a finding

Everything banned in Part 2, plus:

- Global UMAP with compartment colors
- A named-subtype UMAP computed by deleting Part 3 points from an old
  embedding
- Cell-level Kruskal–Wallis / Wilcoxon stars
- High-tail odds ratios as the lead figure
- Dropping `NO_EXPOSURE_RANGE` rows so the heatmap looks powered

