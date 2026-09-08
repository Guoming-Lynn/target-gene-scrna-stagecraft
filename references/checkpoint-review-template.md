# Human checkpoint review template

The reviewer must fill this table before the agent continues. A blank, copied,
or target-gene-only rationale is a stop condition.

| Field | Required entry |
|---|---|
| stage/round | exact stage and round identifier |
| reviewer_id / second_reviewer_id | two independent reviewers when feasible |
| review_date | ISO date |
| selected_resolution | one grid value, or `STOPPED` |
| evidence_files | marker/QC/contamination/coverage table paths and hashes |
| biological rationale | marker programs and QC evidence in plain language |
| alternative interpretation | plausible competing label or reason unresolved |
| target-gene blinding | `BLIND`, `UNBLINDED_WITH_REASON`, or `NOT_APPLICABLE` |
| decision | `KEEP`, `DELETE`, `NAME`, or `STOPPED` |
| uncertainty | confidence and what would change the decision |

For Part 1 Checkpoint B every cluster needs one row. For Part 3 STOP POINT 2
every cluster needs KEEP or DELETE with a reason, and no subtype name is legal.
For STOP POINT 3 the mapping must be complete and DELETE count must be zero.
If reviewers disagree, preserve both rows, adjudicate explicitly, and run the
prewritten alternative-label sensitivity if the disagreement changes a Part 4
or Part 5 population. An agent may not invent a biological label to fill a blank.

