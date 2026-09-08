# Worked patterns

These are abstract worked patterns from multi-study scRNA-seq analyses. Numbers are illustrative of *shape*, not a license to copy results into another paper.

## 1. Zero at the donor, hits at the cell

A symptom-status contrast on myeloid donor-pseudobulk returned **zero** strict DEGs. A cell-level Wilcoxon screen on the same object produced hits.

Keep the screen in an `exploratory` folder. Do not promote it in the README as the finding. The biological n never became the cell.

## 2. Identity DEG is not the gene-of-interest story

Two donor-paired identity contrasts (foam-like vs inflammatory monocyte; SPP1+ vs inflammatory monocyte) produced hundreds of `robust_primary` genes. The gene of interest was **not** a strict DEG in either.

That is a result. Do not recycle the identity table as "pathways of gene X". Do not merge the two FDR families. They are parallel estimands on the same donors.

## 3. Geometry is a gate, not a gene list

A monocyte–LINEAGE_C continuum dropped cDC and cycling cells **before** HVG. PAGA had a frozen threshold. Failure at that threshold would have stopped pseudotime rather than lowered the cutoff.

Donor-level late-vs-early pseudobulk is a **different** FDR family from the identity contrasts. Draft folders that only recomputed UMAP/PAGA have no gene claims.

## 4. Tiny q, one study

A tissue slope had q ~ 1e-18, donor LOO 20/20, and three accessions. Most donors were in two consecutive accessions with interleaved donor IDs. Holding out that block left residual df 0. The external accession had no exposure range.

Verdict: `SINGLE_SOURCE_DEPENDENT`. The slope inside the block did not shrink. What died was the sentence "independently replicated".

## 5. Orthogonal bulk does not pay the replication debt

An external treatment co-induced two genes across independent bulk series. That changes the **mechanistic prior** (parallel pathway targets rather than A causes B). It does not change `SINGLE_SOURCE_DEPENDENT` on the tissue slope.

A later module-score test must freeze the module in the external data *before* touching the atlas, and must still refuse to call itself a new cohort.

## 6. Collinear exposures

Detection fraction and log2CPM of the same gene at donor-unit level can have Spearman > 0.9. Then a second model on intensity is a scale check, not a new study. If negative-control genes move even more than the highlighted gene on that scale, the narrative belongs on the arm where the controls are weaker.

## 7. Communication without an estimable donor model

A CellChat screen produced many recurrent edges. Every pre-registered donor-level model was `NOT_ESTIMABLE`. Verdict: no Tier 1 axis. Do not pick an experiment from the pretty network.

## 8. Virtual perturbation is an embedding shift

In silico KO/OE moved a frozen embedding. Significant ends were OE-only, and
KO vs OE used different cell sets. That is not a paired mechanistic mirror
and not an expression prediction.

Virtual KO is Part 6, not a Part 5 sequel. It does not rerun limma.
Truncation cells are technical provenance. BH family size is frozen; a
missing test does not shrink it. A 17/20 donor sign test does not pay the
source-block debt from Part 5.

## 9. Incomplete genetic negative controls

A perturbation RNA-seq arm still showed residual target counts and residual induction. Do not treat that arm as a clean falsifier of the unperturbed co-induction. Label it incomplete.

## 10. Protocol for a cheap model

The designer writes the protocol (question, gates, verdict table, forbidden sentences). The executor hashes inputs, reproduces the parent coefficient, writes `NOT_ESTIMABLE` when df dies, and does not download a new GSE because the report "needs replication".

## 11. Rounds first, names last

A lineage cleanup that names clusters, then deletes two of them, then keeps
the old UMAP with holes, is not a recluster. The surviving numeric IDs after
an HVG recompute are new IDs. Name them only when the KEEP/DELETE sheet has
zero DELETE rows. If two rounds still leave mixed-lineage debris, stopping
the lineage is a result.

## 12. Identifiability before the slope

A subtype survey that reports a high cell violin and then fits limma is two
parts glued together. Write the unit-level exposure sd first. If a source
has detection sd ≈ 0, say `NO_EXPOSURE_RANGE` before anyone opens R. That
forecast is not a fitted `NOT_ESTIMABLE` verdict.

## 13. Donor LOO is not source-block LODO

A slope can have q ~ 1e-18, donor LOO 20/20, and three GEO accessions, and
still be one study. Consecutive accessions with interleaved donor IDs are
one `source_block`. Holding that block out can leave residual df 0. Write
`SINGLE_SOURCE_DEPENDENT`. The within-block slope may not have moved. What
died is the sentence "independently replicated".

## 14. Folder 08 is not Part 5

Positive vs undetected cells at the barcode is the anti-pattern that had
to be rebuilt as a donor-continuous model. Keep any cell Wilcoxon in an
exploratory folder. The finding is the donor-unit table, including
`NOT_ESTIMABLE`.


