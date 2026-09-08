# Choose the next stage

Start from the artifacts you already have. Open only the linked stage contract
and the relevant glossary entries before execution. Missing prerequisites mean
return to their producing stage, not substitute a cell-level test.

| Current input and question | Next step | Required handoff |
|---|---|---|
| Library counts; no locked atlas | [Part 1](part1-qc-and-global-atlas.md) | Human-selected clustering and complete labels |
| Locked atlas; where is the target measured? | [Part 2](part2-target-gene-survey.md) | Descriptive donor-unit summary; no testing |
| Named compartment needs cleanup/subtypes | [Part 3](part3-compartment-recluster.md) | Human KEEP/DELETE, recompute, then lock names |
| Locked subtypes; which have usable exposure range? | [Part 4](part4-subtype-survey.md) | Forecast, not an association verdict |
| Locked subtypes and a prespecified association question | [Part 5](part5-donor-association.md) | Project arm inventory, frozen design, donor/source audits |
| Locked subtypes; model embedding response to virtual KO/OE | [Part 6](part6-virtual-knockout.md) | Separate pinned environment, model parity, exploratory ceiling |
| Other estimands or mechanistic extensions | Separate protocol | Parts 7/8 are not implemented here |

For a new executor, first run `python quickstart.py --out new_smoke_output`.
Before the chosen stage run `python scripts/check_environment.py --stage partN`.
Smoke success proves neither a working Geneformer model nor scientific validity.
For a Part 4 forecast, provide exactly one biological donor row per subtype;
use `--donor-key` on `part4_identifiability.py` for a frozen crosswalk identity.
The default is `dataset_donor_id`; a namespaced string does not prove that
the same donor was not sampled across accessions. Select `--source-key
source_block` when the input contains a frozen source map. Dataset counts alone
do not establish source independence. Missing/invalid exposures stop the check.
See [INSTALL.md](../INSTALL.md) for interpreter selection. Human label decisions
remain required; the short route does not authorize an agent to invent them.

## Before the first Part 5 arm

Freeze a project-wide list, including planned subtype arms and alternative
specifications, in `project_arms.yaml` at the project root:

```yaml
arms:
  - id: ARM_A
    verdict: analysis/arm_a/05_logs/verdict.json
  - id: ARM_B
    verdict: analysis/arm_b/05_logs/verdict.json
```

Hash this file in each arm's frozen inputs. After running arms:

```bash
python scripts/project_arm_inventory.py project_arms.yaml --out project_arm_report.json
```

The report retains every declared slot; missing outputs give exit 2. Preserve
failed and inconclusive verdicts. Disclose added arms and protocol amendments,
and report the inventory alongside any selected result. The helper cannot
discover undeclared runs elsewhere or prove the manifest was frozen in advance.

Within-arm BH does not control project-wide FDR or the probability that any arm
passes. Do not combine q values, treat pass tokens as p values, or apply an
independence formula to overlapping arms. A project-wide discovery claim needs
a prospectively defined joint/hierarchical testing strategy and implementation;
this helper does not implement that correction. Without one, conclusions stay
within arm and project multiplicity remains `NOT_ESTABLISHED`.
