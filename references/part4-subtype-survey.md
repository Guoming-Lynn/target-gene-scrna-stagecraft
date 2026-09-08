# Part 4 — Subtype TARGET_GENE survey

Audience: an agent with a **Part 3 locked compartment**. Subtype names are
frozen. The compartment UMAP is frozen. This chapter is Part 2 again, on
that label key.

It is not DEG. It is not a new clustering. It is not "the gene marks
subtype X". It does not reopen Part 3 labels.

Operating figures: [part4-figures.md](part4-figures.md).
Methods that are not restated here: [part2-target-gene-survey.md](part2-target-gene-survey.md).
Renderer: `scripts/part4_figures.py` (wraps `part2_figures.py` + identifiability).

---

## 0. What this part answers

In the locked subtypes of **one** Part 3 object:

1. Which subtypes even contain cells with `TARGET_GENE` count > 0?
2. Do detection, positive-cell intensity, all-cell mean, and UMI
   contribution still disagree at the subtype grain?
3. At `dataset × donor_id` × subtype, is there enough **exposure range**
   that a Part 5 slope could even be identified?

It does **not** answer association, enrichment, or replication. Those
words wait for Part 5. A Part 4 flag `NO_EXPOSURE_RANGE` is a forecast,
not a fitted verdict.

Claim type: descriptive survey. Evidence ceiling: descriptive.
Part 5 may start after this only for subtypes/sources the identifiability
table did not kill. A stopped Part 3 lineage (no names) does **not** get
a Part 4.

---

## 1. Inputs (read-only)

- The Part 3 final annotated h5ad. Hash it. Do not write into it.
- Frozen subtype color YAML from Part 3 (`exact AnnData string → #hex`).
- `layers["counts"]`, `layers["normalized"]`, **this object's**
  `obsm["X_umap"]` (the compartment embedding, not the global atlas).
- Locked subtype key (no NA). Numeric Leiden key may exist for
  provenance; grouping is the **name** key.
- `obs`: `dataset`, `donor_id`. Optional `source_block`. If
  `dataset_donor_id` is absent, build it in memory.

Assert `TARGET_GENE` is in `var_names`. If not, stop.

Do not recompute HVG, PCA, Harmony, neighbors, or UMAP.
Do not subset away `Unresolved` / `stressed` states unless the protocol
froze that exclusion in Part 3. Describe them; flag them as unlikely
Part 5 arms.

One Part 3 object per Part 4 run. Another lineage is another Part 4, not
a concatenated survey.

---

## 2. Definitions

Same four endpoints, same unit, same eligibility as Part 2. Copy them;
do not invent a fifth interchangeable mean.

| Name | Source | Formula |
|---|---|---|
| Detected cell | `counts` | gene UMI > 0 |
| Detection fraction | counts | n_detected / n_cells |
| Positive-cell intensity | `normalized` | mean among detected cells |
| All-cell mean | `normalized` | mean including zeros |
| Count contribution | `counts` | subtype gene-UMI / compartment gene-UMI |
| Pseudobulk CPM | `counts` | 1e6 × gene UMI / total UMI in a unit × subtype |
| Biological unit | metadata | `dataset × donor_id` |
| Eligible unit × subtype | counts | n_cells ≥ 20 (sensitivity 5 / 10 / 50) |

A raw zero is non-detection in this library, not a biological off-state.

---

## 3. What to compute

1. **Part 2 tables**, grouped by the locked subtype key.
2. **Identifiability table** (this chapter's extra job), one row per
   subtype × source (`source_block` if present, else `dataset`), plus a
   subtype roll-up across sources.

Starting freeze for the forecast (amend in the protocol):

| Gate | Formal | Exploratory | Dead |
|---|---|---|---|
| Eligible units | ≥ 12 | ≥ 8 | < 8 → `LOW_N` |
| Detection sd | ≥ 0.02 | same | both scales below freeze → `NO_EXPOSURE_RANGE` |
| log2(CPM+1) sd | ≥ 0.15 | same | (with detection) |
| Datasets / source_blocks | ≥ 2 | 1 with range → `WITHIN_SOURCE_RANGE` | — |

First matching flag wins, in this order:

```text
LOW_N
NO_EXPOSURE_RANGE
WITHIN_SOURCE_RANGE   # pooled / ALL row only, when a single source has range
EXPLORATORY_N
RANGE_OK
```

Per-source rows never use `WITHIN_SOURCE_RANGE`; that token is the pooled
forecast ("a slope may exist inside one source, not across sources").

If the subtype name matches unresolved / stressed / doublet / debris
(case-insensitive), append `UNLIKELY_PART5_ARM` in a notes column. That
does not hide the row.

Also write, per subtype on eligible units: Spearman between detection
and log2(CPM+1). `|ρ| ≥ 0.9` → note `COLLINEAR_SCALES`. Part 5 may still
run one primary exposure; the second scale is not a second study.

Do not fit limma. Do not compute residual df of a design that does not
exist yet. Do not Wilcoxon-on-cells as a finding. If a human asks for
Kruskal–Wallis / high-tail OR, park it in `tables/exploratory/` and do
not quote it in the one-page survey.

---

## 4. What Part 4 may claim

- Subtype A has the highest detection fraction (number + n_cells).
- Subtype B contributes the largest share of `TARGET_GENE` UMIs in this
  compartment.
- Among eligible donor-units, subtype C has the highest median detection.
- Detection and positive-cell intensity do not rank subtypes the same way.
- After holding out dataset D, the remaining ranking still / no longer
  puts subtype A first (**descriptive**).
- Source S has detection sd below the freeze in subtype C, so a Part 5
  slope in that slice is forecast `NO_EXPOSURE_RANGE`.

## 5. What Part 4 may not claim

- `TARGET_GENE` marks / defines subtype X
- `TARGET_GENE` is differentially expressed or **associated** with a
  transcriptome (that word is Part 5)
- Leave-one-dataset-out = independent replication
- High-expression tail OR = a new subtype or a mechanism
- `RANGE_OK` = the Part 5 model will pass
- `NO_EXPOSURE_RANGE` on one source, therefore the gene is uninteresting
- Adjacent tissue is a healthy control
- N cells = n
- The global atlas UMAP is this compartment

---

## 6. Lifecycle

```text
Part 4 progress:
- [ ] Confirm Part 3 lock (named subtypes, color YAML, counts/normalized/UMAP)
- [ ] Hash the input h5ad; open read-only; assert TARGET_GENE in var_names
- [ ] Render F04_01–14 (same geometry as F02, new prefix, subtype colors)
- [ ] Identifiability table + F04_15 (subtype × source)
- [ ] One-page survey: which subtypes are on, which endpoints disagree,
      which slices cannot identify a Part 5 slope
- [ ] Recommend Part 5 arms — or none
```

Default split: cheap model runs `part4_figures.py`; expensive model writes
the one-page survey from the tables, not from a cell violin.

