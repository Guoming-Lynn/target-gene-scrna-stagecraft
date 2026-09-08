# Statistical stack

Read [research-validity.md](research-validity.md) before selecting a model.
Its estimand and inference qualifications govern the defaults below.

## The unit

Independent biological unit: `dataset × donor_id` (string `dataset_donor_id`).

Cells, droplets, technical libraries, and repeated wells are subsamples. Plot them. Do not star them.

When a donor contributes several cell states, the unit may be `dataset_donor_id × cell_type` with `duplicateCorrelation(block = dataset_donor_id)`. That is still donor-aware. It is not cell-iid.

## Default formal gates

| Gate | Formal | Exploratory | Insufficient |
|---|---|---|---|
| Eligible units | ≥12 | ≥8 | below |
| Datasets | ≥3 | ≥2 if protocol allows | 1 unless protocol is explicitly single-source |
| Residual df | ≥6 | ≥4 | below |
| Design | full rank | full rank | rank-deficient → `NOT_ESTIMABLE` |
| Slope exposure | IQR / sd large enough to identify a slope | same, or label `NO_EXPOSURE_RANGE` | sd ≈ 0 |

A source with almost no exposure variance cannot rescue a slope, even with thousands of cells.

## Primary gene model

Recommended default (R):

1. `filterByExpr` or CPM ≥ 1 in ≥ max(3, 20% of units)
2. TMM `calcNormFactors` on the **target-excluded** count matrix
3. `duplicateCorrelation(..., block = dataset_donor_id)`
4. `voomWithQualityWeights(block=, correlation=)`
5. `lmFit(block=, correlation=)` → `eBayes(robust=TRUE)`

Python-only is acceptable if the protocol names the engine and does not pretend it is this stack.

Support model: edgeR QL, same design. `robust_primary` genes need the same sign.

Typical identity-contrast design:

```text
~ donor + arm + z_log1p_n_cells + z_log1p_mean_umi
```

Typical continuous-exposure design:

```text
~ cell_type + dataset + z_depth + z_foreign + exposure
```

Drop `dataset` from the formula when only one dataset level remains. Record the deviation. Do not keep a constant column.

## Pathways

- Primary: CAMERA on the voom design (competitive).
- Secondary: fgsea multilevel on the primary moderated-t ranking (`minSize=10`, `maxSize=500`).
- Libraries: frozen local JSON, SHA-256 checked (Hallmark / KEGG / GO BP or the protocol's list).
- ORA: only for a pre-specified strict non-technical DEG query of 10–500 genes and ≤20% of the tested universe.

`robust_primary` pathway: CAMERA project q < 0.05 **and** fgsea project q < 0.05, same direction, all successful LODO folds same sign, technical leading-edge < 0.50.

Identity DE that moves thousands of genes often empties CAMERA. That is expected. Do not switch to a friendlier library to obtain a paragraph.

Cell-level AUCell / ssGSEA q-values are not pathway discovery.

## Technical genes

Prefix `MT-`, `RPS`, `RPL`; hemoglobin genes; `MALAT1`; `XIST`.

Measure them. Do not call them strict DEGs. If they dominate a leading edge (>50%), downgrade the pathway.

## Exposures

Dropout-prone gene detection (Jeffreys):

```text
p = (n_detected + 0.5) / (n_cells + 1)
report per +10 percentage points
```

Alternative: donor-unit log2CPM of the same gene.

If Spearman(detection, intensity) > 0.9 at the unit level, passing the second model is not orthogonal support for "not a dropout artifact" beyond a weak statement. Protocol must say so.

## Matching and other secondaries

1:1 matching of high vs low exposure within strata is a **secondary** track. It never replaces the primary continuous model. If matching coverage is poor, retire it in an amendment; do not keep a dead track in the README.

## Multiple testing

Name the family: all genes in this estimand, or a frozen panel (`q_panel`).

Do not:

- merge hit lists from two stages and speak of one discovery
- put LODO / threshold / leave-one-donor rows into the project FDR family
- use a sensitivity model's q as a veto unless the protocol said it was a gate

Non-significant ≠ equivalent. Do not "prove" an external signature false because an atlas contrast did not pass q.

## Collinearity and scores

If adding a module score as a covariate, pre-register:

- Spearman(exposure, score) and VIF gates (e.g. Spearman > 0.95 or VIF > 10 → `COLLINEAR_UNINTERPRETABLE`)
- a shrink fraction for the coefficient of interest
- negative-control genes that should not collapse

Do not build the module from genes chosen by looking at the same atlas matrix.

## Virtual knockout (Part 6)

Unit remains `dataset × donor_id`. The test is a two-sided exact sign test
on **donor medians** of Δaxis, not a cell mixed model.

- Family size is declared (primary endpoints × declared perturbations).
- BH uses that size even when a cell of the family is NA.
- Near-zero `|median| ≤ 1e-6` is excluded from the sign count.
- Controls are matched comparators, not a permutation null.

Δaxis is a unit-sphere projection onto a leave-one-donor-out axis. It is
not a log2FC and not a pathway NES.


