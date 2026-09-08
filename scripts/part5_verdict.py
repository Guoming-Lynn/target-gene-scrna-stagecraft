#!/usr/bin/env python3
"""Apply a pre-registered ordered verdict table. First matching row wins.

Does not refit models. Reads model_audit.json (and optional extra flags).

Usage:
    python scripts/part5_verdict.py 05_logs/model_audit.json \\
        --table 00_protocol_manifest/verdict_table.yaml \\
        --out 05_logs/verdict.json
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Mapping

try:
    import yaml
except ImportError:  # pragma: no cover
    yaml = None


def _require_new(path: Path) -> None:
    if path.exists():
        raise SystemExit(f"Refusing to overwrite: {path}")


def _as_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    if isinstance(value, (int, float)):
        return bool(value)
    return str(value).strip().lower() in {"true", "yes", "1", "success"}


def load_audit(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise SystemExit("model_audit.json must be an object")
    return payload


def load_table(path: Path) -> list[dict[str, Any]]:
    if yaml is None:
        raise SystemExit("PyYAML required: pip install pyyaml")
    raw = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    rows = raw.get("rows", raw) if isinstance(raw, dict) else raw
    if not isinstance(rows, list) or not rows:
        raise SystemExit("verdict table needs a non-empty 'rows' list")
    return rows


def matches(when: Mapping[str, Any] | None, audit: Mapping[str, Any]) -> bool:
    if not when:
        return True
    for key, expected in when.items():
        observed = audit.get(key)
        if observed is None:
            return False
        if isinstance(expected, bool) or str(expected).lower() in {"true", "false"}:
            if _as_bool(observed) != _as_bool(expected):
                return False
            continue
        if observed is None:
            return False
        if str(observed) != str(expected):
            return False
    return True


def apply_table(rows: list[dict[str, Any]], audit: Mapping[str, Any]) -> dict[str, Any]:
    # Validate every row, including unreachable rows, before evaluating matches.
    # These legacy tokens falsely promoted internal holdouts to external evidence.
    forbidden = {"EXTERNALLY_REPLICATED", "EXTERNALLY_CONSISTENT_UNDERPOWERED"}
    if any(str(row.get("token", "")).strip() in forbidden for row in rows):
        raise SystemExit("Internal sensitivity cannot emit external-replication tokens")
    for index, row in enumerate(rows, 1):
        token = str(row.get("token") or "").strip()
        if not token:
            raise SystemExit(f"verdict row {index} has no token")
        if matches(row.get("when") or {}, audit):
            return {
                "status": "SUCCESS",
                "verdict": token,
                "reason": str(row.get("reason") or audit.get("reason") or f"first matching row: {token}"),
                "matched_row": index,
                "tags": list(audit.get("tags") or []),
                "evidence_ceiling": audit.get("evidence_ceiling"),
                "n_units": audit.get("n_units"),
                "n_donors": audit.get("n_donors"),
                "n_datasets": audit.get("n_datasets"),
                "n_source_blocks": audit.get("n_source_blocks"),
                "rdf": audit.get("rdf"),
                "drop_dominant_status": audit.get("drop_dominant_status"),
                "full_status": audit.get("full_status"),
                "full_reproduced": audit.get("full_reproduced"),
            }
    return {
        "status": "SUCCESS",
        "verdict": "INCONCLUSIVE",
        "reason": "no verdict row matched",
        "matched_row": None,
        "tags": list(audit.get("tags") or []),
        "evidence_ceiling": audit.get("evidence_ceiling"),
        "n_units": audit.get("n_units"),
        "n_donors": audit.get("n_donors"),
        "n_datasets": audit.get("n_datasets"),
        "n_source_blocks": audit.get("n_source_blocks"),
        "rdf": audit.get("rdf"),
        "drop_dominant_status": audit.get("drop_dominant_status"),
        "full_status": audit.get("full_status"),
        "full_reproduced": audit.get("full_reproduced"),
    }


def guarded_verdict(rows, audit):
    verdict = apply_table(rows, audit)
    if verdict['verdict'].startswith('FROZEN_PASS'):
        if audit.get('source_dependent') is True or audit.get('drop_dominant_status') == 'NOT_ESTIMABLE':
            verdict.update(verdict='SINGLE_SOURCE_DEPENDENT', reason='Source dependence overrides the ordered table')
        elif audit.get('estimand_mode') == 'joint_common_slope':
            verdict.update(verdict='INCONCLUSIVE', reason='Repeated-donor mode remains exploratory pending calibration and valid method support')
        elif not all(audit.get(k) is True for k in ('lodo_complete', 'formal_gates_pass')) or audit.get('source_dependent') is not False:
            verdict.update(verdict='INCONCLUSIVE', reason='Complete source sensitivity and formal gates required')
        if not verdict['verdict'].startswith('FROZEN_PASS'):
            verdict['matched_row'] = None
    verdict['source_evidence'] = 'INTERNAL_SENSITIVITY_ONLY'
    verdict['external_validation'] = 'NOT_ASSESSED'
    verdict['calibration_status'] = 'NOT_ESTABLISHED_FOR_THIS_ANALYSIS'
    verdict['scientifically_calibrated'] = False
    verdict['inference_scope'] = 'WITHIN_ARM_ONLY'
    verdict['project_multiplicity_control'] = 'NOT_ESTABLISHED'
    verdict['whole_pipeline_calibrated'] = False
    verdict['uncalibrated_components'] = ['QC_selection', 'annotation_selection', 'source_mapping', 'pathways']
    verdict['model_mode_status'] = audit.get('model_mode_status', 'NOT_ASSESSED')
    if audit.get('estimand_mode') == 'joint_common_slope':
        verdict['evidence_ceiling'] = 'exploratory'
        verdict['model_mode_status'] = 'EXPLORATORY_ONLY_CALIBRATION_CONCERN'
    verdict['precision_status'] = audit.get('precision_status', 'NOT_ASSESSED')
    verdict['precision_summary'] = audit.get('precision_summary')
    verdict['protocol_chronology'] = audit.get('protocol_chronology', 'NOT_VERIFIED')
    verdict['interpretation_limit'] = (
        'Protocol gates are not empirical FDR calibration, external validation, '
        'or evidence of adequate power. Read effects and intervals.'
    )
    return verdict


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("audit", type=Path)
    parser.add_argument("--table", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)

    verdict = guarded_verdict(load_table(args.table), load_audit(args.audit))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    _require_new(args.out)
    args.out.write_text(json.dumps(verdict, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"verdict={verdict['verdict']} row={verdict['matched_row']}")
    print("First matching row won. Do not paraphrase the token into a stronger claim.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

