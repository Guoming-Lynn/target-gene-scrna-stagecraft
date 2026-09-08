# Part 6 — Virtual knockout of TARGET_GENE (embedding shift)

Read [research-validity.md](research-validity.md), especially shared-axis dependence.
Binomial sign p values below are nominal unless independence or calibration is
established; official implementation parity does not establish biological validity.

Even consistent embedding shifts across donors do not validate the biological
effectiveness of the Geneformer perturbation module itself. Neither a sign-test
pass nor matched-control separation demonstrates prediction of experimental
perturbation responses. That requires an independent experimental benchmark
matched to the tissue, task and perturbation, with suitable baselines and no
evaluation leakage. The verdict records
`perturbation_biological_validity: NOT_ESTABLISHED_BY_EMBEDDING_SHIFT`.

Audience: an agent with a **Part 3 locked subtype** (exact string) and a
**hashed, pre-specified endpoint gene set**. This is the first branch after
the 1–5 spine that is still about **that one gene**.

It is not a DEG. It is not an expression predictor. It is not a real
knockout. It does not reopen labels. It does not pay a Part 5 replication
debt.

Figures: [part6-figures.md](part6-figures.md).
Verdict grammar: [evidence-and-verdicts.md](evidence-and-verdicts.md).
Claim lint: [claim-boundaries.md](claim-boundaries.md).

Helpers:

```text
scripts/part6_endpoints.py
scripts/part6_token_audit.py
scripts/part6_controls.py
scripts/part6_axes.py
scripts/part6_eligibility.py
scripts/part6_sign_tests.py
scripts/part6_smoke_gate.py
scripts/part6_verdict.py
scripts/part6_figures.py
```

The transformer is **official Geneformer**, pinned. Python helpers do not
reimplement KO/OE token surgery. A homegrown perturber is not this chapter
until a smoke table against official `InSilicoPerturber(emb_mode="cls")`
passes the frozen cosine gate.

Origin: folders `19` (myeloid) and `21` (EC). Those runs taught the rules
below. Copy the **rules**, not the PATHWAY_ENDPOINT numbers.

---

## 0. What this part answers

In **one** locked subtype, after deleting the `TARGET_GENE` token from the
official rank-value sequence (virtual KO), does the last-layer CLS
embedding move along a **pre-frozen, donor-cross-fitted** state axis in a
direction that is consistent across `dataset × donor_id`?

Optional second estimand, only if the protocol lists it **before** smoke:
virtual overexpression (OE) — move/insert the token to the front of the
gene-token ranks. OE is **not** a required mechanistic mirror of KO.

Claim type: orthogonal / geometric (embedding). Evidence ceiling:
`exploratory` or `geometry_only`. `can_only_downgrade: true`.

This chapter may **not** use the word **associated** as if it were Part 5.
It may not use **causes**, **regulates**, **predicts expression**,
**knockout phenotype**, **OCR**, **ATP**, **respiration**, or **in-vivo**.

Anti-patterns this chapter exists to block:

- treating cosine displacement as up/down without a frozen axis
- drawing KO–OE pairing lines across different cell sets
- cell-level p-values on embeddings
- swapping the perturber when Geneformer is inconvenient
- using a significant embedding shift to upgrade
  `SINGLE_SOURCE_DEPENDENT`

---

## 1. When it may start

Required:

- Part 3 named lock. Exact subtype string. No NA on that key.
- `TARGET_GENE` in `var_names`. Unique symbol.
- Endpoint membership files hashed in the protocol **before** any CLS is
  written.
- A compute budget (wall-clock and disk) written in the protocol.

Recommended: a finished Part 5 on the same subtype, so the endpoints are
not shopped from a volcano on the same object. If Part 5 did not run,
endpoints must be an **external** hashed GMT / literature list, and the
ceiling stays `exploratory`.

Do not start Part 6 because F02/F04 violins look high.
Do not start it to "confirm" a Part 5 gene hit. The axis is a **gene set**,
not a second test of `TARGET_GENE` itself.

A stopped Part 3 lineage does not receive Part 6.

---

## 2. Inputs (read-only)

- Part 3 final annotated h5ad. Hash it. Do not write into it.
- Exact subtype value (character-for-character).
- `layers["counts"]` — tokenization, detection, control matching.
  Non-negative integer UMI. `.raw` is not a substitute.
- `layers["normalized"]` — **unperturbed** baseline endpoint scores only.
  Never a predicted-expression layer.
- `obs`: `dataset`, `donor_id`. Build `dataset_donor_id` in memory if
  absent. Assert it equals `dataset × donor_id`.
- Frozen endpoint member files + SHA-256.
- Optional: Part 5 `verdict.json` (quoted as a ceiling, not as a parent
  coefficient to reproduce).

Do not recompute HVG, PCA, Harmony, neighbors, UMAP, or Leiden.
Do not relabel, merge, or split subtypes.
Do not read a removal archive to "add cells".
Do not download a friendlier checkpoint.

---

## 3. STOP POINT 0 — freeze, then install nothing else

The protocol lists one row:

| Field | Rule |
|---|---|
| Population | exact locked subtype string |
| Unit | `dataset_donor_id` |
| Min cells / unit (cohort) | 10 (starting freeze; amend to 20 if Part 5 used 20) |
| Perturbations | `KO` required. `OE` only if listed here |
| Primary endpoints | 1–3 hashed gene sets. Not `TARGET_GENE` |
| Support endpoints | optional; descriptive; not in the FDR family |
| Direction control | optional named support set (e.g. a complex expected to move with the primary, or a set expected not to) |
| Model | official Geneformer, pinned revision + checkpoint + dictionary hashes |
| Ceiling | `exploratory` or `geometry_only` |
| Parent Part 5 token | copy it, or `PART5_NOT_RUN` |

Human freeze of this table is a stop. The agent does not add a third
endpoint after looking at Δaxis.

If parent Part 5 is `NOT_ESTIMABLE` or `SINGLE_SOURCE_DEPENDENT`, write
that token into the Part 6 header. Part 6 cannot replace it.

---

## 4. Model contract (write these hashes, then do not wander)

Starting freeze (amend only with new hashes; do not "upgrade" mid-run):

| Piece | Starting freeze |
|---|---|
| Source | official `ctheodoris/Geneformer` |
| Revision | protocol pin (example origin: `04c2b2e84da7c0f385c3f9ad8f3ec24bab6650e5`) |
| Checkpoint | `Geneformer-V2-104M` unless the protocol names another **one** |
| Mode | pretrained inference, `eval()`, no train, no fine-tune |
| Sequence | V2 max 4,096 including `<cls>` / `<eos>` → 4,094 gene tokens |
| Representation | last-layer CLS; resolved hidden-state index must be recorded (12-layer 104M → 11) |
| Seed | 42 |
| `nproc` | 1 |
| Batch | GPU 8; **one** allowed resource drop to 2. CPU 2 if no GPU |
| Input IDs | Ensembl, from the **pinned** name→ID dictionary |
| Counts | unnormalized raw UMI from `layers["counts"]` |

Prohibited without a new protocol (a new folder):

- V2-316M or any second checkpoint "to see if it agrees"
- layer scan (`emb_layer` shopping)
- hyperparameter search
- combinatorial multi-gene perturbation
- whole-transcriptome perturbation
- fine-tuning on this atlas
- GEARS / scGen / CellOracle / scTenifoldKnk as a **fallback**

`scTenifoldKnk` and friends are a different claim type. If the protocol
says `NOT EXECUTED`, write that README line and stop. A Geneformer smoke
failure is `STOPPED`, not an invitation to switch engines.

Input mapping gates (starting freeze):

- `TARGET_GENE` maps to exactly one Ensembl ID and one vocabulary token
- mapped genes in model vocabulary ≥ 70% of input genes
- median per-cell retained raw UMI fraction ≥ 80%
- `.X` of the Geneformer input object is non-negative integer counts
- duplicate symbols→same Ensembl: sum raw counts per cell, log the
  collapse; do not pick the alias that looks better
- no network alias search

The Geneformer `gene_mapping_file` is Ensembl→Ensembl. It is not a
symbol→Ensembl table. Use the pinned `gene_name_id_dict.pkl` for symbols.

---

## 5. Cohort

1. Subset the exact subtype string.
2. Recompute `n_cells` per `dataset_donor_id` **from this object**.
3. Keep units with `n_cells >=` the freeze. Log excluded units. Do not
   lower the threshold to rescue n.
4. Tokenize only the formal cohort.

If a cross-check metadata table from Part 5 exists, compare unit counts
cell-for-cell. Conflict → stop. Do not pick the larger table.

---

## 6. Token observability (KO and OE are different populations)

After official V2 tokenization, every formal cell gets a ledger row:

| Column | Meaning |
|---|---|
| `raw_count` | `TARGET_GENE` UMI in `layers["counts"]` after mapping/collapse |
| `raw_detected` | `raw_count > 0` |
| `final_token_present` | token in the **final** 4,096-length sequence |
| `sequence_length` | official length including specials |
| `pretruncation_rank` | rank of the gene token before the 4,094 cap |
| `cap_status` | at cap or not |

Rules. These are estimand definitions, not QC trivia.

| Observation | Action |
|---|---|
| raw+ and token present | KO eligible; OE eligible; OE-symmetry eligible |
| raw+ and token absent, sequence at 4,096, pretruncation rank > 4,094 | **truncation**. KO ineligible. OE primary eligible. Not a biological subgroup |
| raw+ and token absent, any other reason | `STOPPED` |
| raw− and token present | `STOPPED` |
| raw− and token absent | KO ineligible; OE primary eligible |

Do not force the token into an original sequence to "save" KO cells.
Do not change the 4,096 cap.
Do not treat truncation cells as a `TARGET_GENE`-low subtype.

**KO population:** final-token-present cells only.

**OE primary** (if declared): every successfully tokenized formal cell.

**OE symmetry** (if OE declared): the KO population, descriptive only. It
does not add a test to the FDR family.

If KO n_cells ≠ OE n_cells, the protocol must say `KO_OE_UNPAIRED`. Figures
must not draw pairing segments between those two columns.

Helper: `scripts/part6_token_audit.py`.

---

## 7. Frozen endpoints (before any perturbation)

Primary endpoints are gene **sets**, hashed, frozen, and few.

Coverage, starting freeze:

| Role | Gate | Failure |
|---|---|---|
| Each primary | ≥ 70% of frozen members model-visible **and** ≥ 10 genes | whole chapter `STOPPED` |
| Support | ≥ 70% and the protocol's min (default ≥ 3) | that endpoint `NOT_ESTIMABLE` |
| Named direction control (e.g. 4/4 subunits) | exact freeze | loses direction-control status; `NOT_ESTIMABLE` |

`TARGET_GENE` is **removed** from every endpoint before scoring, even if
the GMT listed it. Record the drop. The gene is the perturbation, not a
coordinate of the axis.

Do not:

- add a leading-edge from Part 5 after seeing Δaxis
- migrate a myeloid-only core into an EC run "for continuity"
- drop an endpoint because it did not move
- score endpoints on perturbed expression (there is none)

Baseline score, unperturbed `layers["normalized"]`, model-visible members
only: unweighted mean. No AUCell. No ssGSEA p-values.

Helper: `scripts/part6_endpoints.py`.

---

## 8. Matched controls (frozen before KO/OE results)

Controls are **descriptive comparators**. They are not a permutation null.
They do not enter the FDR family. They do not prove specificity.

Pool: model-vocabulary genes in the formal cohort, raw-count detection
and mean. Exclude:

- `TARGET_GENE`
- the union of all frozen endpoint members (after target drop)
- `MT-`, `RPS`, `RPL`, `MRPS`, `MRPL`
- cell-cycle core: `MKI67, TOP2A, PCNA, MCM2–7, CDK1, CCNA2, CCNB1,
  CCNB2, UBE2C, BIRC5, CENPF, TYMS`
- strong TF/stress families matching
  `^(JUN|FOS|STAT|IRF|NFKB|CEBP|KLF|ZNF|GATA|SPI|RUNX|MYC|EGR|ATF|CREB|FOXO|HIF|REL)`
- technical extremes: `ACTB, GAPDH, B2M, MALAT1, XIST, TMSB10, TMSB4X,
  EEF1A1`
- additional project-specific exclusions: freeze a symbol list and pass
  `--exclude-genes`; do not embed a project's target family in the public defaults
- mean-count top/bottom 1% in the candidate pool

Match window (starting freeze):

1. `|Δ detection| ≤ 0.05` and `|log2 mean-count ratio| ≤ 0.5`
2. If fewer than 10 genes, **one** expansion: `0.10` / `1.0`
3. No second expansion

Distance, z-scored on the eligible pool:

```text
x1 = logit(detection clipped to [1e-6, 1-1e-6])
x2 = log1p(mean raw counts / cell)
distance = Euclidean(z(x1), z(x2)) to TARGET_GENE
```

Take up to 10 by distance, then symbol. Seed 42 is a documented no-op
unless a true tie remains.

| n controls | Consequence |
|---|---|
| 10 | full descriptive ranks |
| 5–9 | continue; disclose reduced resolution |
| < 5 | TARGET_GENE KO continues; control ranks `NOT_ESTIMABLE` |

Do not replace a control after seeing Δaxis. Each control uses **its own**
token-present set for KO.

Helper: `scripts/part6_controls.py`.

---

## 9. Smoke gate (unlocks the scientific run)

Before the formal cohort is perturbed:

1. Take the 3 formal donors with the most cells. Seed 42. Sample ≤ 100
   cells per donor. Total ≤ 300. Only `TARGET_GENE` KO (and OE if
   declared).
2. Token semantics: KO deletes the token; OE places it first after
   `<cls>` without duplicating it; specials intact.
3. Embeddings finite; cell order preserved; shape = n × hidden.
4. Deterministic rerun, max absolute embedding difference ≤ `1e-6`.
5. Independent wrapper vs official `InSilicoPerturber(emb_mode="cls")`:
   max absolute cosine difference ≤ `1e-5` for **both** KO and OE (if
   OE is in the protocol).
6. OE comparator **must** use the official overflow reference:
   `overexpress_tokens` + `calc_n_overflow` +
   `truncate_by_n_overflow_special` on the unperturbed sequence before
   CLS. A wrapper that skips this will fail OE parity while KO looks
   fine. That is a harness bug, not biology. Origin: Amendment 03.

Allowed repair: **one** input/API/environment fix, or batch 8→2, then
exactly one rerun. Second failure → `STOPPED`. Do not relax `1e-5`.
Do not switch checkpoint. Do not drop OE from the protocol to make
smoke pass unless the human amends the protocol **and** restarts.

If smoke wall-clock extrapolates to > 24 h or > 25 GB new disk for the
formal run → compute `STOPPED`. Do not subsample biology to sneak under
the cap.

Helper: `scripts/part6_smoke_gate.py` (judges a parity table; it does
not download the model).

---

## 10. Cross-fitted axes (no leakage, no post-hoc direction)

Official cosine tells you **that** the embedding moved, not **which
way** relative to a biological set. The axis is frozen on **unperturbed**
cells before KO/OE.

For each endpoint E that passed coverage, for each donor d with ≥ 10
tokenized cells:

1. Baseline score from unperturbed `normalized`, model-visible members,
   unweighted mean.
2. Inside d, sort by score; ties by `cell_id`. Top and bottom size
   `ceil(0.25 × n)` , each ≥ 3. If that fails, donor d is
   `NOT_ESTIMABLE` for E.
3. L2-normalize each cell CLS.
4. `a(d,E) = L2[ mean(CLS_top) − mean(CLS_bottom) ]`
5. For cell i in donor d:
   `a(−d,E) = L2[ mean_{j≠d} a(j,E) ]`
   Training axes ≥ 10. Primary endpoints need ≥ 11 valid donor axes
   before any perturbation. Otherwise `STOPPED` (do not merge donors).

Sparse `normalized` must be indexed as `X[rows, :][:, columns]`
(two-stage). Pairing `X[rows, cols]` on CSR is a silent wrong axis.
Origin: prior block-indexing audit.

Cell-level effect:

```text
Δaxis(i,E) = dot( L2(CLS_perturbed,i) − L2(CLS_original,i) , a(−d,E) )
```

- `Δaxis > 0` — embedding moved toward the **unperturbed baseline-high**
  state of E
- `Δaxis < 0` — toward baseline-low
- The number is a unit-sphere projection. It is not log2FC, not OCR,
  not ATP, not a probability of differentiation

Also store `1 - cosine(CLS_original, CLS_perturbed)` as unsigned QC.
It is not the primary endpoint.

Axis construction must not use perturbed CLS. If it did, stop; the
direction was fitted on the answer.

Helper: `scripts/part6_axes.py`.

---

## 11. Formal perturbations

Targets: `TARGET_GENE` + frozen controls. One gene at a time.

Cache original CLS **once**. Hash the cache. Perturb against that cache.

Keep failed rows. Do not drop them and pretend n is the successes.
Keep ineligible KO rows (token absent) as `INELIGIBLE`, not as Δaxis = 0.

All `Δaxis` on `RUN` rows must be finite.

Sham: unmodified sequences, one re-inference of the whole formal
cohort. Max embedding absolute difference and max cosine distance
≤ `1e-6`. Failure → do not interpret Δaxis as biology; write
`SHAM_DRIFT`. Do not zero-fill.

---

## 12. Donor eligibility, then summaries (never the reverse)

Build the eligibility ledger **before** across-donor medians, ranks, or
p-values. Origin: prior correction of an earlier implementation
gap.

| Population | Donor minimum |
|---|---|
| KO | ≥ 5 final-token-present **successful** perturbations |
| OE primary | ≥ 10 successful tokenized perturbations |
| OE symmetry | ≥ 5 token-present successful perturbations |

Ineligible donors stay in the ledger with a reason. They do not enter
medians, sign counts, control ranks, or p-values. Do not borrow cells
from another donor.

Inside an eligible donor: median / Q1 / Q3 / IQR of cell Δaxis.
Across donors: **equal weight**. Do not weight by n_cells.

Near-zero donor median: `|median| ≤ 1e-6`. Excluded from the sign test,
still counted in n_evaluable.

Helper: `scripts/part6_eligibility.py`.

---

## 13. Primary tests (closed family)

The only inferential tests are:

```text
each primary endpoint × each declared perturbation of TARGET_GENE
```

KO only, 2 primaries → family size 2.
KO + OE, 2 primaries → family size 4.

That integer is frozen. If one cell of the family is `NOT_ESTIMABLE`,
**do not shrink** the family. BH uses the declared size. Origin: Stage
an earlier script that BH-adjusted only the finite p-values.

Test: two-sided exact binomial sign test on **donor medians** (not
cells). H0: P(positive) = 1/2 among non-near-zero donors. Require ≥ 5
non-near-zero eligible donors; else NA + reason.

Support endpoints, controls, sham, OE symmetry: **description only**.

No cell-level p-values. No Wilcoxon on barcodes. No mixed model that
treats cells as iid.

Directional rank vs controls (descriptive): KO more positive = better
rank if the protocol's pre-registered KO direction is positive; OE more
negative = better rank if the pre-registered OE direction is negative.
Absolute rank is `|across-donor median|`. Report **both**. A gene can
rank 2/11 on direction and 11/11 on magnitude. That is the result, not
a license to hide the absolute rank. Origin: prior implementation audit.

Pre-registered directions (if any) are for **caption checking**, not for
one-sided tests, not for dropping endpoints.

Helper: `scripts/part6_sign_tests.py`.

---

## 14. Source structure (say it; do not launder it)

A 17/20 donor sign test is **direction consistency inside this object**.
It is not source-block replaceability. Consecutive GEO accessions with
interleaved donor IDs remain one study.

Part 6 does **not** re-run Geneformer under source-block LODO by default
(that is a new, expensive protocol). It **does** require:

- a source_block column on the donor ledger (from Part 5's map, or a
  freeze in this protocol)
- a sentence in the report: sign consistency ≠ independent replication
- if parent Part 5 is `SINGLE_SOURCE_DEPENDENT`, Part 6 verdict cannot
  be paraphrased as "now independently replicated in the model"

---

## 15. Pre-registered verdict table (first matching row wins)

Starting table (copy, then amend). Tokens are closed; do not invent
English synonyms in the report.

| Order | Token | When |
|---|---|---|
| 1 | `STOPPED` | hash / label / mapping / primary coverage / smoke / axis / identity blocker |
| 2 | `TOKEN_UNOBSERVABLE` | `TARGET_GENE` has no KO-eligible cell after audit |
| 3 | `NOT_ESTIMABLE` | primary family has no test with ≥ 5 non-near-zero donors |
| 4 | `SHAM_DRIFT` | sham gate failed |
| 5 | `KO_OE_UNPAIRED` | OE was declared **and** KO vs OE cell sets differ **and** the report drew them as a paired mechanism (lint). Prefer: keep this as a **tag** on an otherwise valid KO result, not a veto of KO |
| 6 | `CONTROL_MAGNITUDE_NOT_SPECIFIC` | optional tag when absolute rank of TARGET_GENE is last among controls; does not veto a sign test |
| 7 | `EMBEDDING_SHIFT_CONSISTENT` | declared family executed; at least the pre-registered KO primary tests pass BH < 0.05 in the pre-registered direction, or the protocol's weaker descriptive rule |
| 8 | `PASS_WITH_LIMITATIONS` | scientific run finished; unpaired populations, control ranks thin, renderer fallback, or parent Part 5 already limited |
| 9 | `INCONCLUSIVE` | no row matched |

`KO_OE_UNPAIRED` as a **tag** is the honest default when OE primary uses
all tokenized cells and KO uses token-present cells. Do not delete OE to
avoid the tag. Do not claim a paired KO↔OE mechanism.

Negative, opposite, or near-zero Δaxis is not `STOPPED`.

Helper: `scripts/part6_verdict.py` (ordered table, first hit).

---

## 16. Directory

```text
analysis/<NN>_virtual_TARGET_GENE_<subtype>/
  00_protocol_manifest/   PROTOCOL, analysis_config.yaml, endpoint hashes,
                          control freeze (after selection, before KO table)
  00_input_audit/         h5ad hash, cohort manifest, mapping, coverage
  02_environment/         revision, checkpoint hashes, CUDA/cuDNN, DLL notes
  03_geneformer/          input, tokenized, smoke, original CLS, axes,
                          perturbations, sham
  04_optional_engines/    README saying NOT EXECUTED (no silent fallback)
  05_controls/            candidate table + frozen controls + SHA-256
  02_tables/              token ledger, donor eligibility, donor effects,
                          sign tests, control ranks
  03_figures/             F06_* PNG/PDF/SVG + source_data/ + sidecars
  05_logs/                smoke_parity.json, sham.json, model_audit.json,
                          verdict.json
  06_reports/             claim boundaries in the same file as the numbers
```

Do not overwrite a frozen smoke-fail directory with a passing rerun.
Amendments are append-only. Prior stopped runs retain their evidence for
Amendments 01–03.

---

## 17. What Part 6 may claim

- In this locked subtype, under this pinned model, virtual deletion of
  the `TARGET_GENE` token shifted last-layer CLS along frozen axis E,
  with N eligible donor medians in the stated direction (exact sign p,
  BH p, family size stated).
- OE, if run, moved (or did not) on its **own** population.
- Truncation left K cells KO-ineligible; they are technical provenance.
- Sham drift was below the freeze (or was not).
- Control-relative directional rank and absolute rank, described as
  ranks among ≤ 10 matched genes, not as a genome-wide specificity test.

## 18. What Part 6 may not claim

- Virtual KO predicts post-perturbation expression
- This is a real knockout / overexpression / CRISPR phenotype
- The shift is OCR, ATP, membrane potential, ROS, complex activity
- EndMT / differentiation / time / fate
- `TARGET_GENE` causes or regulates the endpoint
- KO and OE are a paired mechanistic mirror when cell sets differ
- Cell-level p-values
- Sign test N/N donors = independently replicated / cross-dataset
- Embedding shift upgrades Part 5 (`SINGLE_SOURCE_DEPENDENT` stays)
- Controls form an empirical null or a permutation p-value
- A second foundation model "confirmed" this run
- Incomplete KO cells (token truncated) are a biological low-expressor
  class

## 19. Explicitly not done here

- Part 5 donor-unit limma (already written; do not rerun it here)
- Frozen **module dose-response** on a new object (Part 8 claim type;
  was historically folders `18` / `29-11`)
- Identity / clinical contrasts (Part 7)
- CellChat, STRING, PAGA, virtual KO of a **different** gene
- Shopping a new GMT because Δaxis was small

---

## 20. Lifecycle

```text
Part 6 progress:
- [ ] Confirm Part 3 lock; quote parent Part 5 token or PART5_NOT_RUN
- [ ] STOP POINT 0: freeze subtype, endpoints+hashes, KO±OE, model pin,
      compute cap, verdict table, forbidden sentences
- [ ] Hash h5ad / dictionaries / checkpoint; abort on mismatch
- [ ] Cohort n_cells gate; do not read removal archives
- [ ] Map symbols → Ensembl on counts; TARGET_GENE unique token
- [ ] Endpoint coverage; drop TARGET_GENE from members
- [ ] Token ledger; STOP on illegal mismatches
- [ ] Axis feasibility (≥11 donor axes) BEFORE inference
- [ ] Freeze matched controls (or mark ranks NOT_ESTIMABLE)
- [ ] Smoke (official parity ≤1e-5, determinism ≤1e-6, OE overflow)
- [ ] Cache original CLS; build leave-one-donor-out axes
- [ ] Formal single-gene KO (±OE) on TARGET_GENE + controls
- [ ] Sham
- [ ] Donor eligibility ledger BEFORE summaries
- [ ] Closed-family sign tests; BH with declared family size
- [ ] verdict.json (first matching row); tags include KO_OE_UNPAIRED
      when populations differ
- [ ] F06 figures: no fake pairing lines; axes say embedding shift
- [ ] Report: token + n_donors + family size + wording consequences
```

Default split: expensive model writes the protocol; cheap model executes
it. Conflict → stop and report. Do not delete OE, shrink BH, or switch
to scTenifoldKnk to force a paragraph.
