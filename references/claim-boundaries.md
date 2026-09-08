# Claim boundaries (do not write)

Use this as a lint list while drafting reports and README bullets. If a sentence matches a left-hand pattern, replace it.

## Replication and n

| Forbidden | Allowed |
|---|---|
| Independently replicated in N donors | N donor-units; source structure stated |
| Cross-dataset consistent | Within-dataset or within-source; holdout result stated |
| LOO N/N therefore general | Donor LOO vs source-block LODO distinguished |
| Cells n = 100,000 therefore powered for DE | Power is in donor-units and residual df |
| Three GEO accessions = three studies | Source-block map required |

## Causality and time

| Forbidden | Allowed |
|---|---|
| Target gene drives / causes / regulates Y in tissue | Association between exposure and Y at donor-unit level |
| In situ differentiation proven | Transcriptional order on a snapshot graph |
| State A is older / later in hours | Higher median DPT; fate probability |
| RNA velocity as evidence without unspliced | Do not run it |
| Virtual KO predicts expression change | Embedding shift under a frozen model, cell sets named |
| KO and OE as a paired mechanistic mirror | Separate populations; `KO_OE_UNPAIRED` when n differs |
| Sign test N/N donors = independently replicated | Direction consistency inside this object |
| Geneformer failed, so we ran scTenifoldKnk | `STOPPED`; no engine fallback |
| Homegrown perturber without official smoke | Wrapper cosine ≤ 1e-5 vs `InSilicoPerturber` |

## Methods dressed as biology

| Forbidden | Allowed |
|---|---|
| Cell Wilcoxon / `rank_genes_groups` as the finding | Screen only; donor pseudobulk is the test |
| Authorized cell fallback as discovery | `EXPLORATORY_ONLY`; authorization and false-positive warning required |
| UMAP distance as relatedness | PAGA connectivity at a frozen threshold |
| STRING hub = mechanism | Topology on a frozen set; hub audit after clustering |
| CellChat edge = experimental axis | Recurrence bar + donor-level estimability |
| Module score from the same matrix that is tested | Module frozen in an external cohort |

## Evidence laundry

| Forbidden | Allowed |
|---|---|
| External treatment lifts source dependence | Orthogonal pharmacological readout; ceiling unchanged |
| Detection and intensity as two datasets | Two scales; collinearity stated |
| High-tail OR as a new subtype | Exploratory concentration; labels stay frozen |
| Part 4 `RANGE_OK` as a fitted association | Forecast only; Part 5 still has to run |
| Source-block omitted, still "independently replicated" | Association unfinished until LODO |
| Empty CAMERA → swap GMT for a paragraph | Empty CAMERA is allowed; do not shop libraries |
| Split-half of one study as external replication | Split-half consistency |
| Non-significant = equivalent / falsifies bulk | Effect, CI, q; bulk already used once |
| Merging FDR hit lists across stages | Separate estimands, optional qualitative overlap |
| Incomplete KO as a clean negative control | Residual expression / exon capture noted |

## Identity and contamination

| Forbidden | Allowed |
|---|---|
| Reopening frozen labels inside a DEG stage | New stage, new object, human lock |
| Global UMAP as compartment proof | Recomputed manifold |
| Algorithm disconnected DC, therefore exclusion was justified | Exclusion was a priori |
| Old Leiden 3 is new Leiden 3 | New IDs after HVG recompute |
| Subtype names on a DELETE-round worksheet | KEEP/DELETE only; names at STOP POINT 3 |
| Child UMAP = parent embedding minus points | Full recompute from retained counts |
| `TARGET_GENE`+ as a subtype name | Lineage markers; gene stays out of labels |

If the user asks to "just quickly" do one of the forbidden moves, refuse the wording even if you run the plot they wanted.


