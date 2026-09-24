#!/usr/bin/env python3
"""Render 06_reports/PART<N>_REPORT.md from verdict.json and model_audit.json.

Formats existing numbers only.

Usage:
    python scripts/stage_report.py analysis/05_arm_a
    python scripts/stage_report.py analysis/06_arm --part 6
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
from stagecraft.io import ensure_repo_on_path as _ensure_repo_on_path  # noqa: E402

_ensure_repo_on_path(__file__)

import pandas as pd  # noqa: E402

from stagecraft import EXIT_GATE, EXIT_OK, EXIT_USAGE, __version__, stop  # noqa: E402
from stagecraft.claims import lint_text  # noqa: E402
from stagecraft.hashing import sha256_file  # noqa: E402
from stagecraft.io import publish_new_files, read_identity_csv  # noqa: E402

_GENE_COLUMNS = (
    "gene",
    "logFC",
    "CI_L",
    "CI_R",
    "q_bh",
    "dual_same_sign",
    "lodo_all_same_sign",
    "robust_primary",
)
_FOLD_COLUMNS = ("subset", "held_out_block", "status", "rdf")
_PATHWAY_COLUMNS = ("library", "pathway", "q_bh_camera", "q_bh_fgsea", "same_direction")
_SIGN_COLUMNS = (
    "endpoint",
    "perturbation",
    "analysis_population",
    "raw_p_value",
    "bh_adjusted_p_value",
    "inference_status",
)
_ELIGIBILITY_KEYS = (
    "full_status",
    "formal_gates_pass",
    "rdf_ok",
    "n_units",
    "n_donors",
    "n_datasets",
    "n_source_blocks",
    "n_blocks",
    "block_key",
    "rdf",
    "estimand_mode",
    "model_mode_status",
)
_PART6_COVERAGE = (
    "ko_eligible_cells",
    "n_eligible_donors_ko",
    "n_sign_tests_run",
    "family_size",
    "smoke_pass",
    "sham_pass",
    "ko_primary_bh_pass",
    "run_finished",
    "blocker",
)
_NOT_CLAIMED = (
    "calibration_status",
    "scientifically_calibrated",
    "external_validation",
    "source_evidence",
    "project_multiplicity_control",
    "whole_pipeline_calibrated",
    "uncalibrated_components",
    "protocol_chronology",
    "inference_scope",
)
_PRECISION_KEYS = (
    "median_ci_half_width",
    "p90_ci_half_width",
    "effect_floor",
    "fraction_ci_half_width_above_effect_floor",
)
_PART5_NO_DIRECTION = frozenset(
    {
        "NOT_ESTIMABLE",
        "INCONCLUSIVE",
        "REPRODUCTION_FAILED",
        "COLLINEAR_UNINTERPRETABLE",
    }
)
_PART6_DIRECTION = frozenset({"EMBEDDING_SHIFT_CONSISTENT", "PASS_WITH_LIMITATIONS"})


def _token(value: object) -> str:
    if isinstance(value, str) and value and set(value) <= set("ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_"):
        return f"`{value}`"
    return _show(value)


def _show(value: object) -> str:
    if value is None:
        return "not reported"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        if math.isnan(value):
            return "NA"
        return f"{value:.4g}"
    if isinstance(value, list):
        return ", ".join(_show(item) for item in value)
    return str(value)


def _cell(value: object) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "NA"
    if isinstance(value, bool):
        return "TRUE" if value else "FALSE"
    if isinstance(value, float):
        return f"{value:.4g}"
    text = str(value).replace("|", "\\|")
    if text.lower() == "nan":
        return "NA"
    return text


def _md_table(frame: pd.DataFrame, columns: tuple[str, ...]) -> str:
    present = [column for column in columns if column in frame.columns]
    if not present:
        return "No matching columns; not reported."
    header = "| " + " | ".join(present) + " |"
    rule = "|" + "|".join("---" for _ in present) + "|"
    rows = []
    for record in frame[present].itertuples(index=False):
        rows.append("| " + " | ".join(_cell(item) for item in record) + " |")
    if not rows:
        return "No rows."
    return "\n".join([header, rule, *rows])


def _lookup(verdict: dict, audit: dict | None):
    def value(key: str) -> object:
        if key in verdict:
            return verdict[key]
        if audit and key in audit:
            return audit[key]
        return None

    return value


def _kv(value, keys: tuple[str, ...]) -> str:
    lines = ["| Field | Value |", "|---|---|"]
    for key in keys:
        lines.append(f"| `{key}` | {_token(value(key))} |")
    return "\n".join(lines)


def _missing(relative: str) -> str:
    return f"`{relative}` not found; not reported."


_TRUE_TOKENS = frozenset({"true", "1", "t", "yes"})


def _true_mask(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False).astype(bool)
    return series.astype(str).str.strip().str.lower().isin(_TRUE_TOKENS)


def _read_table(path: Path) -> pd.DataFrame:
    try:
        frame = read_identity_csv(path, dtype={"gene": str})
    except (OSError, UnicodeError, ValueError) as exc:
        stop(f"{path} could not be read: {exc}", EXIT_USAGE)
    for column in ("logFC", "CI_L", "CI_R", "q_bh", "q_bh_camera", "q_bh_fgsea", "rdf"):
        if column in frame.columns:
            frame[column] = pd.to_numeric(frame[column], errors="coerce")
    return frame


def _gene_table(frame: pd.DataFrame, top: int) -> str:
    ordered = frame
    if "robust_primary" in ordered.columns:
        ordered = ordered.assign(_rank=~_true_mask(ordered["robust_primary"]))
        sort_by = ["_rank"]
        ascending = [True]
        if "q_bh" in ordered.columns:
            sort_by.append("q_bh")
            ascending.append(True)
        ordered = ordered.sort_values(sort_by, ascending=ascending, na_position="last")
    elif "q_bh" in ordered.columns:
        ordered = ordered.sort_values("q_bh", ascending=True, na_position="last")
    return _md_table(ordered.head(top), _GENE_COLUMNS)


def _wording(rows: list[tuple[str, str]]) -> str:
    lines = ["| Allowed | Not allowed |", "|---|---|"]
    for allowed, forbidden in rows:
        lines.append(f"| {allowed} | `{forbidden}` |")
    return "\n".join(lines)


def _part5_wording(value, verdict: str) -> str:
    coverage = (
        f"{_show(value('n_donors'))} donor-units across {_show(value('n_source_blocks'))} source blocks",
        "Independent replication in N donors",
    )
    if verdict in _PART5_NO_DIRECTION:
        rows = [
            (
                f"Report `{verdict}` and the stop reason only",
                "A directional association, replication claim, or equivalence",
            ),
            coverage,
        ]
        if verdict == "NOT_ESTIMABLE":
            rows.append(("`NOT_ESTIMABLE` is the result of this arm", "No association, negative result, or equivalence"))
        if verdict in {"INCONCLUSIVE", "REPRODUCTION_FAILED"}:
            rows.append(("The frozen rules did not reach a verdict", "Any directional conclusion"))
        if verdict == "COLLINEAR_UNINTERPRETABLE" or value("collinearity_review_required") is True:
            rows.append(("Two scales of one measurement", "Two independent readouts or datasets"))
        return _wording(rows)
    rows = [
        (
            "Association between exposure and outcome at the donor-unit level, within this arm",
            "Causal, regulatory, or driver language",
        ),
        coverage,
    ]
    if verdict == "SINGLE_SOURCE_DEPENDENT" or value("source_dependent") is True:
        rows.append(("Consistent within one source block", "Cross-dataset or cross-source consistency"))
    if value("evidence_ceiling") == "exploratory" or value("estimand_mode") == "joint_common_slope":
        rows.append(("Exploratory association", "Discovery or formal finding"))
    if value("collinearity_review_required") is True:
        rows.append(("Two scales of one measurement", "Two independent readouts or datasets"))
    if verdict.startswith("FROZEN_PASS"):
        rows.append(("Frozen protocol gates held", "Empirically calibrated FDR, adequate power, or validation"))
    if verdict == "PARTIALLY_CONFOUNDED":
        rows.append(
            (
                "Association survives the primary design; a frozen confound check did not clear",
                "Confound-free association",
            )
        )
    return _wording(rows)


def _part6_wording(value, tags: list, verdict: str) -> str:
    ceiling = ("The Part 5 ceiling is unchanged", "Upgrading a Part 5 source-block ceiling")
    if verdict not in _PART6_DIRECTION:
        return _wording(
            [
                (
                    f"Report `{verdict}` and the stop reason only",
                    "An embedding shift, direction consistency, or biological effect",
                ),
                ceiling,
            ]
        )
    rows = [
        (
            "Embedding shift along the frozen axis under the frozen model",
            "Predicted expression change or biological effect of the perturbation",
        ),
        (
            "Direction consistency across eligible donors in this object",
            "Replication",
        ),
    ]
    if "KO_OE_UNPAIRED" in tags or value("ko_oe_unpaired") is True:
        rows.append(("KO and OE reported as separate populations", "A paired KO/OE mechanistic mirror"))
    rows.append(ceiling)
    return _wording(rows)


def _not_claimed(value) -> str:
    lines = [f"- `{key}`: {_token(value(key))}" for key in _NOT_CLAIMED]
    limit = value("interpretation_limit")
    if isinstance(limit, str) and limit.strip():
        lines.append("")
        lines.append("> " + limit.replace("\n", "\n> "))
    lines.append("")
    lines.append("Project-wide FDR across arms is not controlled here; see scripts/project_arm_inventory.py.")
    return "\n".join(lines)


def _relative_input(stage: Path, path: Path) -> str:
    try:
        return path.resolve().relative_to(stage.resolve()).as_posix()
    except ValueError:
        return path.resolve().as_posix()


def _provenance(stage: Path, paths: list[Path]) -> str:
    lines = [
        f"Generated by stagecraft {__version__}. Numbers are copied from the files below.",
        "",
        "| Input | sha256 |",
        "|---|---|",
    ]
    for path in paths:
        lines.append(f"| `{_relative_input(stage, path)}` | `{sha256_file(path)}` |")
    return "\n".join(lines)


def render_report(stage: Path, part: int, tables: str, top: int) -> tuple[str, list[Path]]:
    verdict_path = stage / "05_logs" / "verdict.json"
    verdict = json.loads(verdict_path.read_text(encoding="utf-8"))
    if not isinstance(verdict, dict) or not isinstance(verdict.get("verdict"), str):
        stop("verdict.json must be an object with a string verdict", EXIT_USAGE)
    audit_path = stage / "05_logs" / "model_audit.json"
    read_paths = [verdict_path]
    audit = None
    if audit_path.is_file():
        try:
            audit = json.loads(audit_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as exc:
            stop(f"model_audit.json is not JSON: {exc}", EXIT_USAGE)
        except UnicodeError as exc:
            stop(f"model_audit.json is not UTF-8: {exc}", EXIT_USAGE)
        if not isinstance(audit, dict):
            stop("model_audit.json must be an object", EXIT_USAGE)
        read_paths.append(audit_path)
    value = _lookup(verdict, audit)
    token = verdict["verdict"]
    ceiling = value("evidence_ceiling")
    ceiling_text = _token(ceiling) if ceiling is not None else "not declared"
    tags = list(verdict.get("tags") or [])
    provenance_at = "@@PROVENANCE@@"
    sections = [
        f"# Part {part} report: {stage.name}",
        "",
        provenance_at,
        "",
        "## 1. Conclusion",
        "",
        f"Verdict: `{token}`. Evidence ceiling: {ceiling_text}.",
        f"Reason: {_show(value('reason'))}",
    ]
    if tags:
        sections.append("Tags: " + ", ".join(f"`{tag}`" for tag in tags) + ".")
    sections.extend(["", "## 2. Reproduction anchor", ""])
    if part == 6:
        parent = value("parent_part5_verdict")
        sections.append(
            "Parent Part 5 verdict: " + (_token(parent) if parent is not None else "not reported") + "."
        )
    else:
        sections.append(f"reproduction_status: {_token(value('reproduction_status'))}")
        sections.append(f"full_reproduced: {_show(value('full_reproduced'))}")
        if value("reproduction_status") == "NOT_APPLICABLE":
            sections.append("No parent coefficient was reproduced in this stage.")
    sections.extend(["", "## 3. Eligibility and coverage", ""])
    if audit is None:
        sections.append(_missing("05_logs/model_audit.json"))
        sections.append("")
    sections.append(_kv(value, _PART6_COVERAGE if part == 6 else _ELIGIBILITY_KEYS))
    sections.extend(["", "## 4. Primary numbers", ""])
    if part == 5:
        sections.append(f"n_robust_primary: {_show(value('n_robust_primary'))}")
        sections.append(f"q_pass: {_show(value('q_pass'))}")
        sections.append(f"dual_method_same_sign: {_show(value('dual_method_same_sign'))}")
        sections.append("")
        gene_path = stage / tables / "gene_evidence.csv"
        if gene_path.is_file():
            sections.append(_gene_table(_read_table(gene_path), top))
            read_paths.append(gene_path)
        else:
            sections.append(_missing(f"{tables}/gene_evidence.csv"))
        sections.append("")
        sections.append(f"precision_status: {_token(value('precision_status'))}")
        summary = value("precision_summary")
        if isinstance(summary, dict):
            for key in _PRECISION_KEYS:
                sections.append(f"{key}: {_show(summary.get(key))}")
            if summary.get("interpretation"):
                sections.append(str(summary["interpretation"]))
        else:
            sections.append("precision_summary: not reported")
        sections.append("")
        sections.append("q values are within-arm BH over the named gene family.")
    else:
        sign_path = stage / tables / "sign_tests.csv"
        if sign_path.is_file():
            frame = _read_table(sign_path)
            preferred = [column for column in _SIGN_COLUMNS if column in frame.columns]
            columns = tuple(preferred or list(frame.columns[:8]))
            sections.append(_md_table(frame, columns))
            read_paths.append(sign_path)
        else:
            sections.append(_missing(f"{tables}/sign_tests.csv"))
    sections.extend(["", "## 5. Robustness", ""])
    if part == 5:
        for key in (
            "source_dependent",
            "drop_dominant_status",
            "lodo_complete",
            "lodo_all_same_sign",
            "collinearity_review_required",
        ):
            sections.append(f"{key}: {_token(value(key))}")
        for key in ("confound_cleared", "sensitivity_pass"):
            if value(key) is not None:
                sections.append(f"{key}: {_show(value(key))}")
        sections.append("")
        fold_path = stage / tables / "fold_audit.csv"
        if fold_path.is_file():
            sections.append(_md_table(_read_table(fold_path), _FOLD_COLUMNS))
            read_paths.append(fold_path)
        else:
            sections.append(_missing(f"{tables}/fold_audit.csv"))
        sections.append("")
        pathway_path = stage / tables / "pathway_evidence.csv"
        if pathway_path.is_file():
            pathways = _read_table(pathway_path)
            if "dual_method_candidate" in pathways.columns:
                count = int(_true_mask(pathways["dual_method_candidate"]).sum())
                sections.append(f"dual_method_candidate rows: {count}")
            if "q_bh_camera" in pathways.columns:
                pathways = pathways.sort_values("q_bh_camera", ascending=True, na_position="last")
            sections.append(_md_table(pathways.head(5), _PATHWAY_COLUMNS))
            sections.append("")
            sections.append("Pathway results are not part of the calibrated components.")
            read_paths.append(pathway_path)
        else:
            sections.append("Pathways were not run or not found.")
    else:
        sections.append("Tags: " + (", ".join(f"`{tag}`" for tag in tags) if tags else "not reported"))
        sections.append(f"KO_OE_UNPAIRED present: {'true' if 'KO_OE_UNPAIRED' in tags else 'false'}")
        sections.append(
            "CONTROL_RANKS_NOT_ESTIMABLE present: "
            + ("true" if "CONTROL_RANKS_NOT_ESTIMABLE" in tags else "false")
        )
    sections.extend(["", "## 6. Wording consequences", ""])
    sections.append(_part6_wording(value, tags, token) if part == 6 else _part5_wording(value, token))
    sections.extend(["", "## 7. Explicitly not shown or claimed", ""])
    if part == 6:
        for key in (
            "evidence_class",
            "causal_inference",
            "expression_prediction",
            "perturbation_biological_validity",
            "evidence_ceiling",
        ):
            sections.append(f"- `{key}`: {_token(value(key))}")
    else:
        sections.append(_not_claimed(value))
    sections.append("")
    body = "\n".join(sections).replace(provenance_at, _provenance(stage, read_paths), 1)
    return body, read_paths


def main(argv: list[str] | None = None) -> int:
    import argparse

    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(errors="backslashreplace")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", type=Path)
    parser.add_argument("--part", choices=("auto", "5", "6"), default="auto")
    parser.add_argument("--tables-dir", default="02_tables")
    parser.add_argument("--top", type=int, default=10)
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args(argv)
    stage = args.stage
    verdict_path = stage / "05_logs" / "verdict.json"
    if not verdict_path.is_file():
        stop("No 05_logs/verdict.json; a report cannot precede the verdict.", EXIT_GATE)
    try:
        verdict = json.loads(verdict_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        stop(f"verdict.json is not JSON: {exc}", EXIT_USAGE)
    if not isinstance(verdict, dict) or not isinstance(verdict.get("verdict"), str):
        stop("verdict.json must be an object with a string verdict", EXIT_USAGE)
    if args.part == "auto":
        part = 6 if verdict.get("evidence_class") == "exploratory_embedding_only" else 5
    else:
        part = int(args.part)
    out = args.out or (stage / "06_reports" / f"PART{part}_REPORT.md")
    text, _paths = render_report(stage, part, args.tables_dir, args.top)

    def write(partials: list[Path]) -> None:
        partials[0].write_text(text, encoding="utf-8")

    publish_new_files([out], write)
    for finding in lint_text(out.read_text(encoding="utf-8")):
        print(
            f'WARNING claim-lint: {out}:{finding.line}:{finding.column}: '
            f'{finding.rule_id}: "{finding.matched}" -> {finding.allowed}'
        )
    print(f"wrote {out}")
    return EXIT_OK


if __name__ == "__main__":
    raise SystemExit(main())
