# Stage layout and audit

## Tree

```text
analysis/<NN>_<short_name>/
  00_protocol_manifest/
    PROTOCOL_CN.md | PROTOCOL.md
    analysis_config.yaml
    input hashes after freeze
    frozen gene lists / GMT pointers
  01_code/
    _common.py
    00_run.py          # orchestrator, no one-off diagnostics
    01_prepare.py
    02_run_models.R    # if the stack is limma
    03_integrate.py
  02_tables/           # or 03_tables; pick one and keep it
  03_figures/
    F<NN>_*.{png,pdf,svg}
    source_data/
    *.parameters.json  # sidecar
  05_logs/
    input_hashes.json
    model_audit.json
    run_audit.json
    verdict.json
    deviations
  06_reports/
    *_REPORT.md        # claim boundaries in the same file as the numbers
  README.md
```

Amendments live in a dated subfolder (`09_reports/amendment_v1_1_0/`) or a new stage. Do not overwrite a frozen report in place and keep the old verdict in the README.

Draft siblings (`26A`, `26B`) that stopped at geometry **must** say so in the folder name or README. They do not inherit `FROZEN_PASS` from the completed sibling.

## Hashes

Before compute, hash every path listed in the protocol. After compute, hash them again. If an upstream file moved, the stage is invalid until the protocol is amended.

`input_hashes.json` maps path → sha256 → role.

## Deviations

If execution cannot follow the protocol (missing column, R package absent, single factor level):

1. Write `DEVIATION_00k_*.md` with the reason
2. State whether the run stops or follows a pre-written fallback
3. Never invent a fallback that changes the estimand

Fallbacks that *were* pre-written (e.g. drop `dataset` when one level remains) are not surprises; still log them.

## Figures

Formal figures: PNG + PDF + SVG, a CSV in `source_data/`, and a sidecar with the exact plotting parameters.

`NOT_ESTIMABLE` rows stay visible (empty forest slot, explicit label). Do not drop them so the plot looks complete.

Do not mix two exposure scales on one axis.

UMAP is display. Empty space is unsampled, not a missing intermediate state.

## Handoff

A completed stage writes what the next stage may read, what it must not reopen, and which sentences are now illegal.

The project README is updated **after** verification, not used as a lab notebook during the run.

Part 3 compartment stages may use round-scoped trees (`objects/round_00_initial/`,
`tables/round_00_initial/`, `figures/round_00_initial/`) beside
`00_protocol_manifest/`. That still counts as a stage layout. Do not clobber
an existing child round.

