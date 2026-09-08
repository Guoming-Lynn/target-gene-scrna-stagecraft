# Pipeline map — parts vs numbered folders

This skill is a **pattern language**, not a replica of one project's 32 folders.

A real program accumulates false starts, stopped lines, draft siblings, and
audits that arrived too late. The parts below are the order you would run **if
you already knew those lessons**. Folder numbers in `Your folders` are the
origin of each pattern. They are not a required directory scheme.

`TARGET_GENE` is the pre-specified gene. Never name a cluster after it.

---

## Design rules

1. **One part = one job.** QC is not DEG. Descriptive detection is not a
   donor-unit slope. Identity DEG is not the gene-of-interest story.
2. **Spine first, branches later.** Parts 1–5 are sequential. Parts 6–8 are
   optional named claim types. A branch may **hold or downgrade** earlier
   wording. It may not upgrade it.
3. **Reuse the loop, do not clone the tissue.** Smooth-muscle, T/NK,
   endothelium, myeloid, and a blood–tissue merge are the **same Part 3**
   with a different subset. Write the loop once.
4. **Source-block audit lives in Part 5**, not as a late cleanup chapter.
   A slope with q ~ 1e-18 and donor-LOO 20/20 can still be
   `SINGLE_SOURCE_DEPENDENT`.
5. **Numbered sibling directories are different estimands.** `26` vs `26A`
   vs `26B`, and a diagnostic archive next to a locked object, do not share
   a verdict.

---

## The six public parts

Parts 7–8 below are historical/future claim types only; they are not
implemented or included in the public workflow.

| Part | Job | Your folders (origin) | Status |
|---|---|---|---|
| **1** | Per-library QC, strict common-gene merge, global Leiden, human labels, locked atlas | `01` QC + `02` global clustering | **written** |
| **2** | Global `TARGET_GENE` survey (detection / intensity / donor-unit description). No DEG. | `03` | **written** |
| **3** | Compartment subset → independent manifold from counts → contamination rounds → human-locked subtypes | `04`+`06`+`07` LINEAGE_A; `09` T/NK; `12`+`12_diag` EC; `15` myeloid; `28` second-tissue merge | **written** |
| **4** | Same survey as Part 2, now on **locked subtypes**. Still no DEG. | `10` T/NK; `13` EC; `16` myeloid | **written** |
| **5** | Donor-unit association, enrichment, confound diagnostics, **mandatory source-block LODO** | `11` v2; `14`; `17`; internals of `29`; **`31` as a gate, not a sequel** | **written** |
| **6** | Virtual KO / optional OE of `TARGET_GENE` on a locked subtype (embedding shift) | `19`; `21` | **written** |
| **7** | Estimands that are **not** the `TARGET_GENE` story (identity, clinical grouping) | `24`; `25`; `27` | later |
| **8** | Named orthogonal claims: trajectory, CellChat, STRING, frozen module, second exposure scale, orthogonal bulk | `18`; `20` `22` `23` `26` `30` `31a` `32`; `14` mito / `29-11` module | later |

Cross-cutting rules (protocol, `n`, layers, verdicts, forbidden sentences)
stay in `SKILL.md`. They are not a ninth part.

---

## Lifecycle (what may start)

```text
Part 1  locked global atlas
   │
   ▼
Part 2  global TARGET_GENE survey  →  choose which compartment is even worth it
   │
   ▼
Part 3  compartment recluster      →  may loop: another lineage, or a second tissue
   │                                →  STOP the line if contamination cannot be cleaned
   ▼
Part 4  subtype TARGET_GENE survey →  is a donor-unit slope even identifiable?
   │
   ▼
Part 5  donor-unit association     →  source-block LODO in the same protocol
   │
   ├── Part 6  virtual KO of TARGET_GENE (embedding shift)
   ├── Part 7  identity or clinical contrast (own FDR family)
   └── Part 8  geometry / communication / topology / frozen module / bulk
```

Do not start Part 5 because someone is impatient after Part 1.
Do not start Part 6 Geneformer before Part 3 labels are locked.
Do not start Part 8 CellChat before Part 3 labels are locked.
Do not treat Part 7 identity tables as "`TARGET_GENE` pathways".
Do not treat a Part 6 sign test as paying the Part 5 source-block debt.

---

## Part-by-part contract (sketch)

Write the full operating spec only when that chapter is opened. This section
is the freeze for **scope**, so later writing cannot quietly merge jobs.

### Part 1 — Atlas foundation

Already specified: [part1-qc-and-global-atlas.md](part1-qc-and-global-atlas.md),
[part1-figures.md](part1-figures.md).

Handoff: gzip annotated h5ad, `layers["counts"]` immutable, human label key,
`TARGET_GENE` present in the common gene set, descriptive detection only.

### Part 2 — Global gene survey

Specified: [part2-target-gene-survey.md](part2-target-gene-survey.md),
[part2-figures.md](part2-figures.md).

Handoff: group + unit tables, F02_01–14, a one-page survey naming which
locked types even express `TARGET_GENE`. No p-values as findings. The next
legal move is a Part 3 recluster of a lineage that is actually on — or stop.

### Part 3 — Compartment recluster

Specified: [part3-compartment-recluster.md](part3-compartment-recluster.md),
[part3-figures.md](part3-figures.md).

Question: after subsetting one (or a merged) lineage from the Part 1 atlas,
what subtypes survive independent HVG → PCA → Harmony → Leiden, and which
clusters are debris / doublets / cross-lineage soup?

This is a **loop**. Smooth muscle, T/NK, endothelium, myeloid, and a
second-tissue merge are repeats, not new methods.

House rule: DELETE in rounds; **recompute HVG from retained counts after
every deletion**; **name only when DELETE = 0**. Intermediate worksheets
are KEEP/DELETE, not names. New Leiden IDs after a recompute are not the
old IDs. Harmony is `dataset` only. `TARGET_GENE` is excluded from HVG/PCA
and from the annotation worksheet.

Handoff: gzip annotated compartment h5ad, `layers["counts"]` +
`layers["normalized"]`, no NA on the final label key, `TARGET_GENE` still
in `var_names`, archived removal reasons. Stop-the-line is an allowed
handoff (exploratory archive, no names).

Variant: second-tissue merge still needs Part 1 QC on the new libraries.
Name-propagate from a pre-removal round is a named protocol variant, not
the default.

### Part 4 — Subtype gene survey

Specified: [part4-subtype-survey.md](part4-subtype-survey.md),
[part4-figures.md](part4-figures.md).

Same methods as Part 2, new label key, Part 3 UMAP and colors. Still
descriptive. Still no DEG.

Handoff: F04_01–14, group + unit tables, identifiability table (F04_15)
forecasting which subtype × source slices cannot identify a Part 5 slope
(`LOW_N`, `NO_EXPOSURE_RANGE`, `WITHIN_SOURCE_RANGE`). Collinear detection
vs log2CPM is a note, not a second study. A stopped Part 3 lineage does
not receive Part 4.

### Part 5 — Donor-unit association

Specified: [part5-donor-association.md](part5-donor-association.md),
[part5-figures.md](part5-figures.md).

Question: at `dataset × donor_id` (optionally × subtype), is `TARGET_GENE`
exposure associated with the rest of the transcriptome, after declared
depth covariates — and what remains after the dominant source block is
held out as a whole?

This is the first part that may use the word **associated**.

Must include in the **same** protocol, not as a later stage:

- biological unit, eligibility gates, residual df of the *actual* design
- target-excluded pseudobulk
- primary limma-voom + support edgeR
- CAMERA primary / fgsea secondary on a hashed GMT
- donor LOO **and** source-block LODO
- collinearity check if a second exposure scale (detection vs log2CPM) is run
- pre-registered verdict table (first matching row wins)
- forbidden sentences (no "independently replicated" without LODO)

`08` (positive vs undetected cells) is **not** this part. It is the
anti-pattern: a cell-level contrast that later had to be rebuilt as `11` v2.

`31` is not "Part 9". If Part 5 omitted source-block, the association is
unfinished.

Handoff: `verdict.json`, gene/pathway tables, F05 holdout forest with
`NOT_ESTIMABLE` slots visible. Part 6 may start only on a named Part 3
object; it may not upgrade this ceiling.

### Part 6 — Virtual knockout of TARGET_GENE

Specified: [part6-virtual-knockout.md](part6-virtual-knockout.md),
[part6-figures.md](part6-figures.md).

Question: in one locked subtype, after deleting the `TARGET_GENE` token
in a **pinned** official Geneformer, does last-layer CLS move along a
pre-frozen, donor-cross-fitted gene-set axis — consistently across
`dataset × donor_id`?

Claim type: embedding shift. Ceiling: `exploratory` / `geometry_only`.
`can_only_downgrade: true`.

Must include in the **same** protocol:

- hashed endpoints before any CLS; `TARGET_GENE` dropped from members
- tokenization from `layers["counts"]`; baseline scores from unperturbed
  `layers["normalized"]`
- KO = final-token-present cells only; cap truncation is technical
- official wrapper smoke (cosine ≤ 1e-5), including OE overflow reference
  if OE is declared
- leave-one-donor-out axes (no self-leakage)
- donor eligibility **before** summaries
- closed-family exact sign test; BH **does not shrink**
- `KO_OE_UNPAIRED` when the two cell sets differ
- no fallback to scTenifoldKnk / GEARS / a second checkpoint

Handoff: `verdict.json`, F06 donor Δaxis panels with embedding-shift
axis labels. A significant KO shift does not upgrade Part 5.

### Part 7 — Non-target estimands

Identity (state A vs state B in the same donors) and clinical grouping
(symptom vs not) are **different FDR families** from Part 5.

Rules:

- `TARGET_GENE` not being a strict DEG **is** a result; do not recycle the
  identity table as "pathways of the gene"
- do not merge hit lists across `24`/`25`/`27`
- a donor-pseudobulk null plus a cell-Wilcoxon hit list: keep the screen in
  an exploratory folder; the finding is the null
- zero results stop mining

### Part 8 — Named orthogonal claims

Each method has a claim type. None of them pay a Part 5 replication debt.

| Method | Allowed | Forbidden |
|---|---|---|
| PAGA / DPT | snapshot order; fate *probabilities* | hours, "older", in-vivo time |
| CellChat | recurrent edges under a frozen bar | experimental axis if donor models are `NOT_ESTIMABLE` |
| STRING | topology on a frozen set | mechanism; post-hoc "we knew this hub" |
| Frozen module score | hashed set; `TIER3_FORCED` allowed to stop | shopping a GMT; upgrading Part 5 |
| Second exposure scale | scale check; `ORTHOGONAL_SUPPORT` only if protocol says so | a new study |
| Orthogonal external assay | pharmacological / mechanistic *prior* | independent in-tissue replication |

Virtual KO/OE lives in **Part 6**, not here. A second foundation model is
still not a Part 8 "replication" of Part 6.

Draft geometry folders (`26A`/`26B`) have **no gene claims**.

---

## What is not a part

| Your folder | Why it is not its own part |
|---|---|
| `05` | absorbed into `04` v2; numbering holes are not missing methods |
| `08` | superseded cell-level DE; cited as anti-pattern inside Part 5 |
| `12_EC_r04_deep_annotation` | diagnostic archive inside Part 3, not a second EC atlas |
| `26A`, `26B` | unfinished drafts; warning inside Part 8 |
| `31` | source-block **gate of Part 5** |
| `29/11`–`29/14` | sub-analyses of one Part 5/6 stage, not new parts |

---

## Suggested writing order

Keep writing one part at a time. Do not draft 6–8 in one pass.

1. **Part 7–8** — remaining optional claim types, each with a hard ceiling

Default split remains: expensive model writes the part spec; cheap model
executes a frozen protocol.

