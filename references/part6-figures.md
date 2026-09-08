# Part 6 figures — virtual knockout (embedding shift)

Shared contracts: [figure-contract.md](figure-contract.md),
[figure-statistics-contract.md](figure-statistics-contract.md), and
[visual-qa-contract.md](visual-qa-contract.md). Shared-axis sign tests are
marked nominal unless their calibration and dependence assumptions are stated.

These panels argue **donor medians of Δaxis**. They do not argue cells.
They do not argue expression. UMAP is not in this catalog.

Theme helper: `scripts/plotting_style.py`. Renderer: `scripts/part6_figures.py`.
Export contract from [part1-figures.md](part1-figures.md): PNG + PDF + SVG +
source CSV + parameter sidecar.

```bash
python scripts/part6_figures.py observability 02_tables/token_ledger.csv --out 03_figures
python scripts/part6_figures.py eligibility 02_tables/donor_eligibility.csv --out 03_figures
python scripts/part6_figures.py donors 02_tables/donor_effects.csv --out 03_figures
python scripts/part6_figures.py forest 02_tables/sign_tests.csv --out 03_figures
python scripts/part6_figures.py ranks 02_tables/control_ranks.csv --out 03_figures
python scripts/part6_figures.py support 02_tables/support_summary.csv --out 03_figures
python scripts/part6_figures.py symmetry 02_tables/donor_effects.csv --out 03_figures
python scripts/part6_figures.py sham 05_logs/sham.json --out 03_figures
```

Stem prefix `F06`. One subtype × one model pin per figure folder.

---

## What is stricter than a Geneformer screenshot

| Rule | Why |
|---|---|
| y-axis says **embedding-axis shift (Δaxis)**, never expression / OCR / ATP / function | the estimand is CLS projection |
| KO and OE are **two columns**. Pairing ticks only if `analysis_population` is the same and donors match | token-present KO vs all-tokenized OE is `KO_OE_UNPAIRED` |
| `NOT_ESTIMABLE` donors stay as empty slots or a labelled gap | dropping them makes 5 donors look like 20 |
| n on the panel is **eligible donor-units**, not cells | cells are subsamples |
| Truncation cells appear on F06_01 as technical, not as a subtype | Amendment 01 |
| Directional rank **and** absolute rank both drawn | Example: rank 2/11 direction, 11/11 magnitude |
| Sham is a QC panel, not a biological result | drift vs biology |

Do not put a cell UMAP next to a donor forest and share an n.
Do not star barcodes.
Do not use a new husl palette for KO/OE; freeze two colours in the YAML.

---

## Required panels

### F06_01 — token observability

Stem `F06_01_token_observability`.

Counts of: raw+/token+, raw+/token− (truncation), raw−/token−, illegal
mismatches (should be zero or the run is STOPPED). Optional: pretruncation
rank vs cap for mismatch cells.

Caption: `KO uses final-token-present cells only`.

### F06_02 — donor eligibility

Stem `F06_02_donor_eligibility`.

One row per `dataset_donor_id`. Bars or a table-as-figure: n cells
tokenized, n KO-eligible, n OE-eligible, pass/fail vs the freeze
(KO ≥ 5, OE ≥ 10). Failed rows stay visible. Color by `source_block`
if the column exists.

Caption: `eligibility before summaries`.

### F06_03 — donor Δaxis distributions (primary)

Stem `F06_03_donor_delta_axis`.

For each primary endpoint, two columns: KO | OE (omit OE if not
declared). Points are **donor medians**. Optional donor IQR as thin
lines. **No segments connecting KO to OE** unless the protocol froze
`KO_OE_PAIRED` (same `analysis_population`, same donors).

A horizontal line at 0. y-label: `embedding-axis shift (Δaxis)`.

Order donors by KO median, and reuse that order in the OE column (donors
missing on one side = gap, not a zero).

### F06_04 — donor-summary forest

Stem `F06_04_primary_forest`.

One row per primary endpoint × perturbation. Point = across-donor
median of donor medians. Interval = central order-statistic nominal 95%
interval for the population median, using all eligible donors including zeros.
`part6_sign_tests.py` reports bounds and achieved coverage under independent,
identically distributed donors. It is not a Walsh-average pseudomedian interval.
With fewer than six donors the 95% interval is unbounded; label that fact and
do not replace it by the observed range. Ties can make coverage conservative.
Shared fitted axes can invalidate the independence assumption: the panel and
sidecar must state `NOMINAL_SHARED_AXIS_DEPENDENCE_NOT_CALIBRATED`.
These marginal intervals do not provide simultaneous family coverage and need
not invert the near-zero-filtered sign test, whose donor subset is different.

Empty row if `NOT_ESTIMABLE`. Annotate n_positive / n_negative /
BH p / family size.

### F06_05 — control ranks

Stem `F06_05_control_ranks`.

Two panels, same genes: **directional rank** and **absolute rank**.
`TARGET_GENE` marked. Controls unlabeled or small. If ranks are
`NOT_ESTIMABLE` (< 5 controls), draw the empty panel with that token.

Caption: `≤10 matched genes; not a permutation p-value`.

### F06_06 — support heatmap (if support endpoints exist)

Stem `F06_06_support_heatmap`.

Rows: support endpoints (direction control tagged). Columns: KO / OE
across-donor median Δaxis. Diverging scale centered at 0. Missing
coverage = blank + `NOT_ESTIMABLE`.

Do not add p-value stars. Support is descriptive.

### F06_07 — KO vs OE symmetry (only if OE was run)

Stem `F06_07_ko_oe_symmetry`.

Points: donors present in **both** KO and OE-symmetry (token-present)
populations. x = KO median, y = OE-symmetry median. Quadrants labelled.
Donors only in one population are listed in the sidecar, not forced onto
the diagonal.

Caption: `opposite sign is descriptive; not an extra FDR test`.

Skip this file if OE was not in the protocol. Do not substitute OE
primary (all tokenized) here.

### F06_08 — sham QC

Stem `F06_08_sham`.

Report max |Δ embedding| and max cosine vs the freeze `1e-6`. A pass/fail
bar is enough. If sham failed, this panel stays in the catalog; biology
panels keep their numbers but the report token is `SHAM_DRIFT`.

---

## Optional display (not discovery)

A single CLS-UMAP of **unperturbed** cells coloured by baseline endpoint
score may live in a supplement folder named `display_only/`. It does not
get an F06 stem. It does not get p-values.

---

## Sidecar

Every stem writes `*.parameters.json` including: subtype string, model
revision, checkpoint, hidden index, family size, KO/OE population names,
n eligible donors, `KO_OE_UNPAIRED` true/false, parent Part 5 token.

