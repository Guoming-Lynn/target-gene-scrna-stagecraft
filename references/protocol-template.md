# Protocol template (freeze this, then code)

Copy into `00_protocol_manifest/PROTOCOL.md` (or `PROTOCOL_CN.md`). Fill every heading. Empty headings are not optional.

Protocol version: `x.y.z`
Frozen on: YYYY-MM-DD
Part: `1_qc_global_atlas` | `2_target_gene_survey` | `3_compartment_recluster` | `4_subtype_survey` | `5_donor_association` | `6_virtual_knockout` | (later parts)
`TARGET_GENE`: SYMBOL (species: …). Pre-specified. Not a discovery target.
Evidence ceiling: `descriptive` | `formal` | `exploratory` | `orthogonal_external` | `geometry_only`
`can_only_downgrade`: true | false
Executor: cheap model must treat this file as the only spec.

For Part 1, follow [part1-qc-and-global-atlas.md](part1-qc-and-global-atlas.md) and stop at human-locked global labels.
For Part 3, follow [part3-compartment-recluster.md](part3-compartment-recluster.md). STOP POINTS are human gates; do not name clusters until DELETE = 0.
For Part 4, follow [part4-subtype-survey.md](part4-subtype-survey.md). Reuse the Part 2 catalog; do not fit a slope.
For Part 5, follow [part5-donor-association.md](part5-donor-association.md). Source-block LODO is in this protocol. `NOT_ESTIMABLE` is a result.

Before modeling, freeze the estimand mode (`subtype_specific` or
`joint_common_slope`), repeated-donor handling, source-block interpretation,
and pre-fit VIF, correlation, and condition-number review thresholds.
For Part 6, follow [part6-virtual-knockout.md](part6-virtual-knockout.md). This is an embedding shift, not predicted expression. KO and OE are unpaired when cell sets differ.

---

## 0. Question

One sentence. Name the population, the contrast or slope, and the unit.

Claim type: associational / geometric ordering / orthogonal pharmacological / audit-only.

Allowed sentences (max 4).
Forbidden sentences (max 6). These are part of the spec, not a style note.

## 1. Hard boundaries

- Upstream directories: read-only. List them.
- Will this stage recluster, relabel, or recompute embeddings? Default no.
- Will this stage download new data? Default no.
- What this stage is *not* allowed to upgrade.

## 2. Inputs (read-only)

Table: purpose, path, expected shape / required columns, SHA-256 after freeze.

Abort if a required column is missing. Do not invent a substitute gene or a substitute annotation key.

## 3. Universe and labels

Which frozen annotation key. Which labels enter, which are excluded, and **why the exclusion is a priori** (not "the algorithm disconnected them").

If a lineage does not belong on the manifold, drop it **before** HVG/PCA, then recompute.

## 4. Unit and eligibility

Unit key. Minimum cells per arm/unit. Formal / exploratory / insufficient gates (units, datasets, residual df, exposure range).

Depth covariates. Block for `duplicateCorrelation`.

## 5. Frozen lists

Genes, negative controls, GMT files with SHA-256. Discovery universe versus pre-specified panels.

The target/exposure gene is excluded from outcome matrices and leading edges unless it is the *tested* gene in an identity contrast (protocol must say which).

## 6. Models

Primary: design formula, software, robust eBayes yes/no.
Support: same design, different engine.
Display-only: cell plots, UMAP, GAM curves. No discovery stars on cells.

One FDR family. Name it. Sensitivity rows do not enter that family.

## 7. Gates that stop interpretation

Identity recovery, PAGA threshold, root rule, reproduction anchor of a previous coefficient, collinearity VIF/Spearman, module sparsity.

Failure → stop biological language; still write coverage and the gate that failed.

## 8. Robustness

LODO by **source_block** if more than one study exists; donor LOO; cell-count thresholds; neighborhood/PC only for geometry.

State what a fold must preserve (sign, not p-value).

## 9. Verdict table

Ordered. First match wins. Include `NOT_ESTIMABLE` and `INCONCLUSIVE`.

Wording consequences: a table mapping current phrases → required replacement phrases.

## 10. Outputs

Tables, figures (PNG/PDF/SVG + source_data + sidecar), `verdict.json` required keys, report path.

## 11. Explicitly not done

Bullet list. "Would be nice" items stay here. They are not executed.

## 12. Run order

Numbered. Checkpoints that require a human look (e.g. UMAP of inherited labels on a new manifold) **stop** until the human passes them.

