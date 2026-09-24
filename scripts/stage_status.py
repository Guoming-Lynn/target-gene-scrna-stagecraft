#!/usr/bin/env python3
"""Show which Part 5 or Part 6 helper steps have produced artifacts and print the next command.

Informational only.

Usage:
    python scripts/stage_status.py analysis/05_arm_a
    python scripts/stage_status.py analysis/06_arm --json status.json
"""
from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
from stagecraft.io import ensure_repo_on_path as _ensure_repo_on_path  # noqa: E402

_ensure_repo_on_path(__file__)

from protocol_chronology import OUTPUT_DIRS, check_chronology  # noqa: E402
from stagecraft import EXIT_OK, EXIT_USAGE  # noqa: E402
from stagecraft.io import publish_new_files  # noqa: E402


@dataclass(frozen=True)
class Step:
    id: str
    artifacts: tuple[str, ...]
    command: str
    optional: bool = False


def _part5() -> list[Step]:
    return [
        Step(
            "protocol",
            ("00_protocol_manifest/PROTOCOL*.md",),
            "Write 00_protocol_manifest/PROTOCOL.md from references/protocol-template.md",
        ),
        Step(
            "config",
            ("00_protocol_manifest/analysis_config.yaml",),
            "Copy scripts/part5_analysis_config.example.yaml to {stage}/00_protocol_manifest/analysis_config.yaml and freeze its values",
        ),
        Step(
            "freeze",
            ("00_protocol_manifest/protocol_freeze.json",),
            "python scripts/check_protocol.py {stage}/00_protocol_manifest/PROTOCOL.md --stage-root {stage} --freeze",
        ),
        Step(
            "pseudobulk",
            ("03_pseudobulk/counts.mtx", "03_pseudobulk/metadata.csv", "03_pseudobulk/genes.csv"),
            "python scripts/part5_pseudobulk.py <locked_subtype.h5ad> --gene <TARGET_GENE> --group subtype --out {stage}/03_pseudobulk",
        ),
        Step(
            "source_blocks",
            ("{tables}/metadata_with_source_block.csv",),
            "python scripts/part5_source_blocks.py {stage}/03_pseudobulk/metadata.csv --map {stage}/00_protocol_manifest/source_block_map.yaml --out {stage}/{tables}",
        ),
        Step(
            "eligibility",
            ("{tables}/eligibility.csv",),
            "python scripts/part5_eligibility.py {stage}/{tables}/metadata_with_source_block.csv --out {stage}/{tables}/eligibility.csv",
        ),
        Step(
            "models",
            ("05_logs/model_audit.json", "{tables}/gene_evidence.csv"),
            "Rscript scripts/part5_run_models.R {stage}/00_protocol_manifest/analysis_config.yaml",
        ),
        Step(
            "pathways",
            ("{tables}/pathway_evidence.csv",),
            "Rscript scripts/part5_run_pathways.R {stage}/00_protocol_manifest/analysis_config.yaml",
            optional=True,
        ),
        Step(
            "verdict",
            ("05_logs/verdict.json",),
            "python scripts/part5_verdict.py {stage}/05_logs/model_audit.json --table {stage}/00_protocol_manifest/verdict_table.yaml --out {stage}/05_logs/verdict.json",
        ),
        Step(
            "report",
            ("06_reports/*_REPORT.md",),
            "python scripts/stage_report.py {stage}\npython scripts/claim_lint.py {stage}/06_reports",
        ),
    ]


def _part6() -> list[Step]:
    return [
        Step(
            "protocol",
            ("00_protocol_manifest/PROTOCOL*.md",),
            "Write 00_protocol_manifest/PROTOCOL.md from references/protocol-template.md",
        ),
        Step(
            "config",
            ("00_protocol_manifest/analysis_config.yaml",),
            "Copy scripts/part6_analysis_config.example.yaml to {stage}/00_protocol_manifest/analysis_config.yaml and freeze its values",
        ),
        Step(
            "freeze",
            ("00_protocol_manifest/protocol_freeze.json",),
            "python scripts/check_protocol.py {stage}/00_protocol_manifest/PROTOCOL.md --stage-root {stage} --freeze",
        ),
        Step(
            "endpoints",
            ("00_input_audit/endpoint_coverage.csv",),
            "python scripts/part6_endpoints.py --sets <endpoints.yaml> --model-genes <model_visible_genes.txt> --target <TARGET_GENE> --out {stage}/00_input_audit/endpoint_coverage.csv",
        ),
        Step(
            "tokens",
            ("{tables}/token_audit.csv",),
            "python scripts/part6_token_audit.py <token_ledger.csv> --out {stage}/{tables}/token_audit.csv",
        ),
        Step(
            "controls",
            ("05_controls/frozen_controls.csv",),
            "python scripts/part6_controls.py <gene_stats.csv> --target <TARGET_GENE> --endpoint-union <endpoint_members.txt> --out {stage}/05_controls/frozen_controls.csv",
        ),
        Step(
            "smoke",
            ("05_logs/smoke_gate.json",),
            "python scripts/part6_smoke_gate.py <parity.csv> --perturbations KO,OE --out {stage}/05_logs/smoke_gate.json",
        ),
        Step(
            "axes",
            ("03_geneformer/axes/*.npz",),
            "python scripts/part6_axes.py --cls <original_cls.npy> --cells <cell_order.csv> --scores <baseline_scores.csv> --endpoint <ENDPOINT> --out {stage}/03_geneformer/axes/<ENDPOINT>.npz",
        ),
        Step(
            "donor_eligibility",
            ("{tables}/donor_eligibility.csv", "{tables}/donor_effects_eligible.csv"),
            "python scripts/part6_eligibility.py <cell_effects.csv> --out {stage}/{tables}/donor_eligibility.csv",
        ),
        Step(
            "sign_tests",
            ("{tables}/sign_tests.csv",),
            "python scripts/part6_sign_tests.py {stage}/{tables}/donor_effects_eligible.csv --target <TARGET_GENE> --family-size <N> --out {stage}/{tables}/sign_tests.csv",
        ),
        Step(
            "verdict",
            ("05_logs/verdict.json",),
            "python scripts/part6_verdict.py <audit.json> --table {stage}/00_protocol_manifest/verdict_table.yaml --out {stage}/05_logs/verdict.json",
        ),
        Step(
            "report",
            ("06_reports/*_REPORT.md",),
            "python scripts/stage_report.py {stage}\npython scripts/claim_lint.py {stage}/06_reports",
        ),
    ]


def _fill(text: str, stage: Path, tables: str) -> str:
    return text.format(stage=stage.as_posix(), tables=tables)


def _exists(stage: Path, pattern: str) -> bool:
    if any(char in pattern for char in "*?[]"):
        if any(stage.glob(pattern)):
            return True
        if pattern.endswith("PROTOCOL*.md"):
            return any(stage.glob(pattern.replace("PROTOCOL*.md", "FROZEN_PROTOCOL*.md")))
        return False
    return (stage / pattern).is_file()


def _has_files(directory: Path) -> bool:
    return directory.is_dir() and any(path.is_file() for path in directory.rglob("*"))


def _load_verdict(stage: Path) -> dict | None:
    path = stage / "05_logs" / "verdict.json"
    if not path.is_file():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    return payload if isinstance(payload, dict) else None


def inspect_stage(stage: Path, part: int, tables: str) -> dict:
    steps = _part6() if part == 6 else _part5()
    recorded = []
    for step in steps:
        missing = []
        shown = []
        for raw in step.artifacts:
            pattern = raw.replace("{tables}", tables)
            shown.append(pattern)
            if not _exists(stage, pattern):
                missing.append(pattern)
        if not missing:
            status = "DONE"
        elif step.optional:
            status = "SKIPPED"
        else:
            status = "MISSING"
        recorded.append(
            {
                "id": step.id,
                "status": status,
                "artifacts": shown,
                "missing": missing,
                "command": _fill(step.command, stage, tables),
                "optional": step.optional,
            }
        )
    warnings = []
    for index, step in enumerate(recorded):
        if step["status"] != "MISSING" or step["optional"]:
            continue
        later = any(item["status"] == "DONE" and not item["optional"] for item in recorded[index + 1 :])
        if later:
            warnings.append(
                f"Later artifacts exist without step {step['id']}; check for skipped or deleted outputs."
            )
    verdict = _load_verdict(stage)
    receipt = stage / "00_protocol_manifest" / "protocol_freeze.json"
    manifest = stage / "00_protocol_manifest"
    protocols = sorted(manifest.glob("PROTOCOL*.md")) or sorted(manifest.glob("FROZEN_PROTOCOL*.md"))
    if receipt.is_file() and protocols:
        try:
            chronology = str(check_chronology(protocols[0], stage))
        except (ValueError, OSError, KeyError, TypeError) as exc:
            chronology = f"FAILED: {exc}"
    elif receipt.is_file():
        chronology = "FAILED: freeze receipt has no PROTOCOL.md or FROZEN_PROTOCOL.md"
    else:
        chronology = "NOT_VERIFIED"
        if any(_has_files(stage / folder) for folder in OUTPUT_DIRS):
            warnings.append(
                "Outputs exist without a freeze receipt; chronology is NOT_VERIFIED "
                "and a new freeze would be refused."
            )
    next_step = next((step for step in recorded if step["status"] == "MISSING"), None)
    return {
        "stage": stage.as_posix(),
        "part": part,
        "steps": [
            {key: step[key] for key in ("id", "status", "artifacts", "missing")}
            for step in recorded
        ],
        "next": None
        if next_step is None
        else {"id": next_step["id"], "command": next_step["command"]},
        "chronology": chronology,
        "verdict": None
        if verdict is None or "verdict" not in verdict
        else {"token": verdict.get("verdict"), "evidence_ceiling": verdict.get("evidence_ceiling")},
        "warnings": warnings,
        "_commands": {step["id"]: step["command"] for step in recorded},
    }


def _detect_part(stage: Path, requested: str) -> int:
    if requested in {"5", "6"}:
        return int(requested)
    verdict = _load_verdict(stage)
    if (stage / "00_input_audit").exists() or (stage / "03_geneformer").exists():
        return 6
    if verdict and verdict.get("evidence_class") == "exploratory_embedding_only":
        return 6
    return 5


def _tables_dir(stage: Path, requested: str) -> str:
    if requested == "02_tables" and not (stage / "02_tables").is_dir() and (stage / "03_tables").is_dir():
        return "03_tables"
    return requested


def main(argv: list[str] | None = None) -> int:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(errors="backslashreplace")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", type=Path)
    parser.add_argument("--part", choices=("auto", "5", "6"), default="auto")
    parser.add_argument("--tables-dir", default="02_tables")
    parser.add_argument("--json", dest="json_out", type=Path, default=None)
    args = parser.parse_args(argv)
    stage = args.stage
    if not stage.is_dir():
        print(f"missing dir: {stage}", file=sys.stderr)
        return EXIT_USAGE
    part = _detect_part(stage, args.part)
    tables = _tables_dir(stage, args.tables_dir)
    report = inspect_stage(stage, part, tables)
    print(f"Stage: {stage.as_posix()} (Part {part})")
    for step in report["steps"]:
        label = f"[{step['status']}]"
        shown = ", ".join(step["artifacts"])
        print(f"{label:<10} {step['id']:<16} {shown}")
    chronology = report["chronology"]
    if chronology == "NOT_VERIFIED":
        print("Chronology: NOT_VERIFIED (no freeze receipt)")
    else:
        print(f"Chronology: {chronology}")
    if report["verdict"] is None:
        print("Verdict: not written")
    else:
        token = report["verdict"]["token"]
        ceiling = report["verdict"]["evidence_ceiling"]
        print(f"Verdict: {token} (ceiling: {ceiling})")
    if report["next"] is None:
        print("Next: none")
    else:
        print(f"Next: {report['next']['id']}")
        for line in report["next"]["command"].splitlines():
            print(f"  {line}")
    if report["warnings"]:
        print("Warnings:")
        for warning in report["warnings"]:
            print(f"  - {warning}")
    print("Parts 1-4 need human-locked decisions; see references/start-here.md.")
    if args.json_out is not None:
        payload = {key: report[key] for key in ("stage", "part", "steps", "next", "chronology", "verdict", "warnings")}
        text = json.dumps(payload, indent=2) + "\n"

        def write(partials: list[Path]) -> None:
            partials[0].write_text(text, encoding="utf-8")

        publish_new_files([args.json_out], write)
    return EXIT_OK


if __name__ == "__main__":
    raise SystemExit(main())
