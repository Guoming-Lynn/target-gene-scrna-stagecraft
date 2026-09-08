# Scientific figure contract

Read this before creating or revising a formal Part 1–6 figure. A figure is a
claim-bearing output, not a decorative export.

## Before plotting

Freeze the scientific question, estimand, biological unit, outcome/exposure
scale, evidence ceiling, panel list, source table, color semantics, journal
size, and caption wording. State one primary conclusion per figure. If a panel
does not support that conclusion or an audit, move it to diagnostics.

Plots must read locked tables or a locked embedding. Do not silently recompute
labels, embeddings, filters, or statistics inside a plotting script. Preserve
missing, `NOT_ESTIMABLE`, and failed-fold rows visibly.

## Export

Every formal panel uses `scripts/plotting_style.py` and produces PNG at 300 dpi,
PDF with embedded Type 42 fonts, SVG unless the sidecar records why it was
skipped, source CSV, and a JSON sidecar. The sidecar records input paths and
SHA-256, configuration hash, software versions, figure size, formats, panel
order, sorting rules, colors, and claim ceiling. UMAP points are rasterized in
vector outputs when dense.

Use one frozen backend for a figure folder. Use colorblind-safe semantic colors,
no random palettes, no rainbow/jet, no 3-D bars, pies, decorative gradients,
or dual y axes. Continuous color requires a labeled colorbar. Categorical
colors must remain stable when panel sorting changes. Include a grayscale check
when color encodes a scientific distinction.

## Statistical display

The figure statistics manifest must state `n`, independent unit, model/test,
error-bar definition, multiple-testing family, effect scale, interval type,
and exact claim ceiling. Error bars are never unlabeled. Do not put cell-level
stars on donor-level figures. Do not display a p value without its test and
family in the caption or sidecar.

## Final QA

Render at final physical size and check: text clipping, missing glyphs, legend
overflow, tick overlap, panel alignment, colorblind/grayscale distinguishability,
empty-slot labels, correct input hashes, and consistency with `verdict.json`.
Any failed check remains in the sidecar and blocks a publication-ready claim.

