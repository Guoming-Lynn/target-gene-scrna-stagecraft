# Beginner glossary

Use these definitions before running a stage. If a term is still unclear, stop
and record `AUTHOR_INPUT_NEEDED`; do not infer a substitute.

| Term | Meaning in this skill | Common mistake |
|---|---|---|
| biological unit / `n` | The independent unit used for inference, normally `dataset × donor_id` | Counting cells as independent replicates |
| counts layer | Raw non-negative integer UMI counts preserved in `layers["counts"]` | Running models on normalized values |
| normalized layer | Log-normalized display values, not the pseudobulk outcome matrix | Treating normalized values as counts |
| MAD | Median absolute deviation rule for within-library QC | Estimating one threshold after mixing libraries |
| Scrublet | Doublet score method run before final MAD filtering | Treating a failed run as all singlets |
| strict common genes | Genes measured in every primary library | Treating union-filled zeros as biological zeros |
| HVG | Highly variable genes used to build a manifold | Selecting HVGs using the target gene story |
| Leiden resolution | A clustering granularity parameter, chosen at a human checkpoint | Choosing the resolution with the prettiest target-gene pattern |
| Harmony | Batch integration of a declared key, normally dataset | Correcting away donor or condition biology |
| KEEP/DELETE | Part 3 review decisions before subtype naming | Naming a cluster and deleting it later while keeping the old UMAP |
| source block | A frozen provenance group used for source holdout | Treating every accession as an independent study |
| LOO | Leave one donor out; leverage/sensitivity analysis | Calling it replication |
| LODO | Leave one source block out; internal sensitivity analysis | Calling it external validation |
| residual df / `rdf` | Information left after fitting the declared design | Reporting q without checking rdf |
| VIF | Variance inflation diagnostic for collinearity | Assuming full rank means stable coefficients |
| CAMERA | Competitive gene-set test using the fitted expression model | Running it on cell-level ranks as discovery |
| fgsea | Gene-set enrichment on a frozen ranked statistic | Shopping gene sets after seeing the result |
| `robust_primary` | A predeclared evidence annotation combining q/effect/support/holdout rules | Treating it as a new FDR-controlled test |
| `NOT_ESTIMABLE` | The declared design cannot answer the question | Replacing it with a convenient cell-level test |
| CLS | Geneformer embedding representation used for geometry | Calling an embedding shift expression fold-change |
| Δaxis | Projection of an embedding shift onto a frozen biological axis | Calling it a pathway score or causal effect |
| nominal p value | A p value whose independence/calibration assumptions are not established | Presenting it as definitive inference |
| sidecar | JSON metadata next to a figure or result | Delivering a plot without inputs, hashes, or parameter record |

The glossary does not replace the stage protocol. It prevents silent vocabulary
errors; the protocol still defines the actual estimand and stop rules.

