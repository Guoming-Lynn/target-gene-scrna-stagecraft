# Research validity contract (1.9.3)

Read before freezing Parts 1, 3, 5, or 6. This contract supersedes older
statements implying that a passing LODO establishes independent replication.
Thresholds below are diagnostic conventions, not universal validity theorems.

## Estimand and biological unit

Freeze population, exposure scale, outcome universe, coefficient/contrast,
covariates and their causal rationale, biological unit, missing-data rule,
multiple-testing family, effect threshold, and permitted conclusion.
Use a donor crosswalk to identify the same person across accessions. A
`dataset_donor_id` string only namespaces IDs; it does not establish independence.
Technical libraries and repeated tissues must not create new independent donors.

Separate subtype slopes have one observation per biological donor within each
subtype. A joint common-slope model assumes the slope is shared across subtypes;
it does not estimate subtype-specific slopes. A joint interaction model must
name the reference subtype, exposure main effect and interaction contrasts.
Within-donor and between-donor exposure effects are different estimands; when
appropriate, separate donor-mean exposure from within-donor deviations.

The bundled Part 5 runner supports `subtype_specific` and explicitly declared
`joint_common_slope`. Interactions require stage-specific code and tests.
Repeated subtype rows require donor correlation modeling. Unblocked edgeR QL
is not valid independent support for repeated donor rows; it stays unavailable
until a justified support design is implemented. Never add donor fixed effects
to a donor-constant exposure model without checking identifiability.

## Identifiability and confounding

Inspect dataset/source by condition, tissue, sex and other prespecified factors.
Perfect confounding means the effects cannot be separated from these data.
Batch correction cannot solve this. Choose adjustment sets from subject-matter
knowledge; mediators and colliders must not be added solely because correlated.

Check missing/non-finite data, rank, donor count, residual df, exposure SD/IQR,
within-source exposure support, exposure VIF and design condition number.
Report the condition number of centered, unit-length non-intercept columns;
raw-unit condition numbers are scale dependent. VIF 10, condition number 30
and absolute Spearman rho 0.90 are review flags, not automatic proof of bias.
Stop on rank failure or no exposure range. Freeze whether a diagnostic flag
requires stopping or documented sensitivity analysis. Do not zero-impute a
missing design value or remove covariates merely to pass a gate.

Detection fraction has binomial measurement uncertainty and depth dependence;
Jeffreys smoothing does not eliminate either. Specify minimum cells before
analysis, report detected/total cells, and assess a frozen alternative exposure
scale or measurement-error sensitivity. Correlated alternative scales are not
independent evidence. Do not use post-hoc observed power; use effect intervals
and design-stage simulation for precision planning.

## Sources and holdouts

LODO is internal sensitivity analysis on overlapping training sets. Its folds
are dependent, and their sign count is not a binomial replication test. No
universal number of sources converts LODO into external validation. One source
offers no source holdout; few sources yield weakly constrained generalization.
Record accession provenance, source independence, donor overlap and unresolved
relationships. Consecutive accessions are clues, not proof of a common study.

For every fold report remaining donors/sources, rank, rdf, exposure SD/IQR,
coefficient, CI, absolute coefficient change and sign. Relative change is
undefined near a zero full coefficient: report NA using a frozen denominator
floor. An inestimable fold is unresolved evidence, not evidence of no effect.
Independent validation requires a genuinely held-out cohort, prespecified
analysis and no reuse for selection, labeling, endpoint or threshold tuning.

## Discovery and selection

Freeze the project-wide number and identity of subtype arms, including planned
alternative specifications, before Part 5. Each new protocol remains part of
that project's multiplicity and selection history. Report all attempted arms,
failures, inconclusive results and amendments; never publish only passing arms.
Use [start-here.md](start-here.md) for the inventory helper and its limitations.
Within-arm BH does not establish project-wide FDR or the probability of at least
one false passing arm. Pass tokens are not p values; shared donors make arm
results dependent. A project discovery claim requires a prespecified joint or
hierarchical testing procedure, not a union of per-arm discovery lists.

Keep discovery (whole-family p/q), practical effect magnitude, and robustness
as separate columns. Filtering a BH discovery set by direction or effect size
does not automatically preserve FDR in the filtered subset. `robust_primary`
is an evidence annotation, not a newly FDR-controlled discovery family.
If the claim concerns effects larger than a minimum, use a prespecified
threshold-null procedure such as TREAT with its own validated implementation;
a zero-null p-value plus an absolute-effect cutoff does not test that claim.

Retain all tested rows, failure slots and original q values. State whether
families span subtypes, genes, pathways and libraries. Do not pool sensitivity
tests into discovery or select the best-performing specification. limma and
edgeR use the same observations, so agreement is method sensitivity, not
replication. Technical-gene exclusions must be defined before inspection;
they are context-dependent, not a claim that these genes lack biology.

## QC, feature coverage, integration and labels

Maintain a sequential attrition ledger by library and biological donor:
input, loose-filter, doublet-tested, singlet, final QC and final annotation
cells; include failure reasons, QC distributions, and measured/detected target
counts at each stage. For an unmeasured gene, detection is NA, not zero.
Compare retention across source/condition and track donors lost entirely.
Differential attrition motivates a selection-bias discussion and prespecified
sensitivity; retention is not automatically an adjustment covariate.

Record library complexity, doublet-method assumptions, loading/library units
and reasons for failure. Inspect sensitivity to plausible QC thresholds without
optimizing the target result. Preserve feature-presence tables and distinguish
unmeasured from measured-zero genes. Report the common-gene testing universe
and gene loss; it is not necessarily the whole transcriptome. Alternative
coverage analyses must not zero-fill missing genes.

Compare integrated and unintegrated representations for marker conservation
and cell-state stability, alongside batch mixing. Strong mixing alone is not
success. Completely confounded biology/batch remains unidentified.
Target exclusion from HVG cannot remove all correlated biological information.

Keep reviewers blinded to target-gene results during label selection. Prefer
two independent reviews before adjudication when feasible; record whether the
same evidence was seen, raw agreement, category prevalence and disagreements.
Kappa is optional for comparable labels and enough clusters, with uncertainty
and its prevalence limitation; it is not annotation accuracy. Single-reviewer
work must be disclosed. Keep uncertain labels and test reasonable alternative
annotations rather than deleting uncertainty to make a clean manifold.

## Part 6 geometry and inference

Report four distinct levels: implementation parity, model geometry, donor
consistency under the analysis assumptions, and external biological validation.
Parity with the official engine establishes only implementation agreement.

Freeze a small justified sensitivity set before inference: for example
quartile versus tercile axes on the same original/perturbed CLS cache. Report
all variants and primary geometry first; transforms define different geometries
and cannot be selected by the smallest p value. Preserve cell/donor alignment.
Reject non-finite embeddings, missing scores and zero/cancelling axes. Constant
scores cannot define high versus low biology via barcode tie-breaking.

Include sham, matched descriptive gene controls, frozen unrelated endpoints
and seeded random-axis comparisons. Random axes and matched genes are not a
calibrated null distribution without justified exchangeability assumptions.
Report effect sizes and uncertainty, not sign counts alone.

Leave-one-donor-out axes share training donors. Resulting donor effects need
not have independent signs, so the nominal exact binomial sign p value is not
automatically an exact inferential p value. Label it nominal until calibrated
or evaluate on donors using an independent fixed training axis. Bootstrap must
resample independent donors/source blocks as justified and rebuild axes inside
each replicate. Bootstrapping cached donor effects alone is conditional on
fitted axes and omits axis uncertainty. Few source blocks limit both methods.

## Calibration and structured evidence

Before calling a new model mode validated, simulate null and non-null count
data through the actual pipeline: donor correlation, source confounding,
depth/exposure collinearity, one-donor influence, dropout, differential QC
retention and failed holdouts. Vary donor/source counts. Freeze seeds, scenario
parameters and rejection rules, and report repetitions and Monte Carlo SE or
binomial intervals for false-positive rate, FDR, power and interval coverage.
Null scenarios must respect the tested estimand; a marginal effect can exist
when a conditional effect is zero. Generic OLS simulation does not validate
limma, Geneformer or the complete QC-selection procedure.

The included synthetic checks verify implementation behavior only. They do not
estimate type-I error or establish scientific calibration.
The separate limited null pilot estimates error only for the simulated
pseudobulk gene-model scenarios. QC choices, annotation decisions, source
grouping, pathways and Part 6 have no pipeline-wide error guarantee from it.
Reviewer agreement, blinding and attrition logs provide audit evidence, not a
numerical false-positive calibration. `formal` never certifies the whole chain.
Every report must carry an evidence matrix with estimand, biological unit,
identifiability, discovery family, effects/CI, holdouts, selection/confounding,
external validation and claim ceiling. Unknown checks remain unknown.
A summary token cannot override missing evidence in this matrix.


Executable pseudobulk-null experiments, precision reporting and local freeze
receipts: [calibration-and-provenance.md](calibration-and-provenance.md).
The null harness does not replace the full scenario suite specified above.
