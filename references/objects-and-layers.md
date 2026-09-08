# Objects, layers, compartments

## Authoritative object

Every compartment has exactly one current input. Downstream stages name it, hash it, and do not replace it with a "similar" h5ad from an earlier round.

Minimum load assertions (adapt numbers to the object):

```python
assert "counts" in adata.layers and "normalized" in adata.layers
assert adata.obs[label_key].isna().sum() == 0
assert adata.obs["dataset_donor_id"].nunique() >= 1
# never: adata.X as counts
```

Keep both a human label (`*_final_cell_type`) and the numeric Leiden key for traceability. Downstream stats use the human label.

## Layer contract

| Layer | Role |
|---|---|
| `layers["counts"]` | Raw UMI after QC. Detection, composition, pseudobulk, limma/edgeR. |
| `layers["normalized"]` | `normalize_total(1e4)` + `log1p`. Color, violins, module display scores. |
| `.X` | Whatever was last written. Unknown until inspected. |
| `obsm["X_pca"]`, Harmony, UMAP | Geometry of **this** object. Not counts. Not transferable after a new subset. |

If `.X` was overwritten by scaled data, still prefer named layers. If a layer is missing, stop; do not silently use `.X`.

## Compartment isolation

Global atlas labels are a **candidate pool**, not a manifold.

1. Subset by the global (or previous) label.
2. Drop contaminants in audit rounds (platelet genes, RBC, off-lineage markers, low-QC islands).
3. From `layers["counts"]`, recompute HVG (batch-aware), PCA, Harmony on `dataset` only, neighbors, Leiden.
4. Human-lock the final labels.
5. Write a handoff that states shape, keys, what was deleted and why, and the hash.

The operating loop (DELETE rounds → new HVG → names last) is
[part3-compartment-recluster.md](part3-compartment-recluster.md).

Do not:

- interpret a compartment on the global UMAP
- Harmony-correct donor or the target gene
- carry DC / cycling cells into a monocyte–LINEAGE_C continuum "to see if they connect" and then treat disconnection as proof of the exclusion
- reopen labels in a DEG stage

## Target exclusion

The exposure gene is not a feature of the manifold and not an outcome of its own slope.

Exclude it from: HVG, PCA, Harmony, neighbors, UMAP, PAGA, root choice, TMM / size-factor estimation, the outcome count matrix of a slope model, and enrichment leading edges.

Part 5 builds that target-excluded matrix from `layers["counts"]` at
`dataset × donor_id` (see [part5-donor-association.md](part5-donor-association.md)).
Do not pseudobulk `normalized`. Do not put the gene back so the volcano has
a hero point.

Part 6 tokenizes `layers["counts"]` into Geneformer and scores frozen
endpoints on **unperturbed** `layers["normalized"]`. There is no
perturbed-expression layer. `TARGET_GENE` is the deleted token, not an
axis member (drop it from endpoint gene sets).

If the stage tests the gene as an **identity DEG** (e.g. state A vs state B), add its raw counts back after building the target-excluded matrix, and say so in the protocol.

## Detection vs biology

A raw count of zero is non-detection under this library's depth. It is not a certified biological negative. High-expression tails, detection fractions, and log2CPM are different readouts; collinear readouts are not independent replications.

## Gene universe

Multi-GEO projects should freeze a **strict common-gene** universe for the atlas that all later stages inherit. Do not silently switch to union genes in a downstream test unless the protocol is a mapping stage.

Gene symbols must be unique. Collapse duplicates at ingest, not at plotting time.

