#!/usr/bin/env python3
"""Apply the Part 6 ordered verdict table. First matching row wins.

Reuses the Part 5 matcher. Extra audit keys (smoke, sham, unpaired) pass through.

Usage:
    python scripts/part6_verdict.py 05_logs/model_audit.json \\
        --table 00_protocol_manifest/verdict_table.yaml \\
        --out 05_logs/verdict.json
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from part5_verdict import apply_table, load_audit, load_table  # noqa: E402

PASS_THROUGH = [
    "blocker",
    "ko_eligible_cells",
    "n_sign_tests_run",
    "sham_pass",
    "ko_primary_bh_pass",
    "run_finished",
    "ko_oe_unpaired",
    "parent_part5_verdict",
    "family_size",
    "smoke_pass",
    "n_eligible_donors_ko",
    "evidence_ceiling",
    "can_only_downgrade",
]


def guarded_verdict(rows, audit):
    verdict = apply_table(rows, audit)
    bools = ('blocker', 'smoke_pass', 'sham_pass', 'run_finished', 'ko_primary_bh_pass')
    counts = ('ko_eligible_cells', 'n_sign_tests_run', 'family_size', 'n_eligible_donors_ko')
    invalid = [k for k in bools if type(audit.get(k)) is not bool]
    invalid += [k for k in counts if type(audit.get(k)) is not int or audit[k] < 0]
    token, reason = None, None
    if invalid:
        token, reason = 'INCONCLUSIVE', 'Missing or invalid audit fields: ' + ', '.join(invalid)
    elif audit.get('blocker') is True or audit.get('smoke_pass') is False:
        token, reason = 'STOPPED', 'Blocker or failed smoke gate'
    elif audit.get('sham_pass') is False:
        token, reason = 'SHAM_DRIFT', 'Sham gate failed'
    elif not audit['run_finished']:
        token, reason = 'INCONCLUSIVE', 'Run has not finished'
    elif audit['ko_eligible_cells'] == 0:
        token, reason = 'TOKEN_UNOBSERVABLE', 'No KO-eligible cells'
    elif audit['n_sign_tests_run'] == 0 or audit['n_eligible_donors_ko'] < 5:
        token, reason = 'NOT_ESTIMABLE', 'Insufficient tests or eligible KO donors'
    elif not 0 < audit['n_sign_tests_run'] <= audit['family_size']:
        token, reason = 'INCONCLUSIVE', 'Invalid frozen test family'
    elif verdict['verdict'] == 'EMBEDDING_SHIFT_CONSISTENT' and not audit['ko_primary_bh_pass']:
        token, reason = 'INCONCLUSIVE', 'KO primary BH gate did not pass'
    if token:
        verdict.update(verdict=token, reason=reason, matched_row=None)
    return verdict


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("audit", type=Path)
    parser.add_argument("--table", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)
    audit = load_audit(args.audit)
    verdict = guarded_verdict(load_table(args.table), audit)
    verdict["evidence_class"] = "exploratory_embedding_only"
    verdict["causal_inference"] = False
    verdict["expression_prediction"] = False
    verdict["perturbation_biological_validity"] = "NOT_ESTABLISHED_BY_EMBEDDING_SHIFT"
    for key in PASS_THROUGH:
        if key in audit:
            verdict[key] = audit[key]
    tags = list(verdict.get("tags") or [])
    if audit.get("ko_oe_unpaired"):
        if "KO_OE_UNPAIRED" not in tags:
            tags.append("KO_OE_UNPAIRED")
    if audit.get("control_ranks_not_estimable"):
        if "CONTROL_RANKS_NOT_ESTIMABLE" not in tags:
            tags.append("CONTROL_RANKS_NOT_ESTIMABLE")
    verdict["tags"] = tags
    verdict["evidence_ceiling"] = "exploratory_embedding_only"
    verdict["can_only_downgrade"] = True
    verdict["ko_oe_unpaired"] = audit.get("ko_oe_unpaired") is not False
    if verdict["ko_oe_unpaired"] and "KO_OE_UNPAIRED" not in tags:
        tags.append("KO_OE_UNPAIRED")
    if args.out.exists():
        raise SystemExit(f"Refusing to overwrite: {args.out}")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(verdict, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"verdict={verdict['verdict']} row={verdict['matched_row']} tags={tags}")
    print("Embedding shift does not upgrade a Part 5 source-block ceiling.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

