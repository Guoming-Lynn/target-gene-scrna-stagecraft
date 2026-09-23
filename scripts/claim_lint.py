#!/usr/bin/env python3
"""Flag report sentences that cross the claim boundaries.

Review aid, not a gate unless --strict.

Usage:
    python scripts/claim_lint.py analysis/05_arm_a/06_reports
    python scripts/claim_lint.py draft.md --strict --json findings.json
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
from stagecraft.io import ensure_repo_on_path as _ensure_repo_on_path  # noqa: E402

_ensure_repo_on_path(__file__)

from stagecraft import EXIT_GATE, EXIT_OK, EXIT_USAGE, stop  # noqa: E402
from stagecraft.claims import load_rules, lint_text  # noqa: E402
from stagecraft.io import publish_new_files  # noqa: E402


def _collect(paths: list[Path]) -> list[Path]:
    files: list[Path] = []
    for path in paths:
        if not path.exists():
            stop(f"missing path: {path}", EXIT_USAGE)
        if path.is_dir():
            files.extend(
                item
                for item in path.rglob("*")
                if item.is_file() and item.suffix.lower() in {".md", ".txt"}
            )
        else:
            files.append(path)
    return sorted(set(files))


def main(argv: list[str] | None = None) -> int:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(errors="backslashreplace")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="+", type=Path)
    parser.add_argument("--rules", type=Path, default=None)
    parser.add_argument("--json", dest="json_out", type=Path, default=None)
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args(argv)
    try:
        rules, settings = load_rules(args.rules)
    except (OSError, ValueError, UnicodeError) as exc:
        stop(str(exc), EXIT_USAGE)
    files = _collect(args.paths)
    findings = []
    for path in files:
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeError as exc:
            stop(f"{path} is not UTF-8: {exc}", EXIT_USAGE)
        for finding in lint_text(text, rules, settings):
            print(
                f'{path}:{finding.line}:{finding.column}: {finding.rule_id}: '
                f'"{finding.matched}" -> {finding.allowed}'
            )
            findings.append(
                {
                    "path": path.as_posix(),
                    "line": finding.line,
                    "column": finding.column,
                    "rule_id": finding.rule_id,
                    "category": finding.category,
                    "matched": finding.matched,
                    "allowed": finding.allowed,
                }
            )
    print(
        f"claim-lint: {len(findings)} finding(s) in {len(files)} file(s). "
        "Review flags only; rewrite, do not suppress."
    )
    if args.json_out is not None:
        payload = {
            "files": [path.as_posix() for path in files],
            "findings": findings,
            "n_findings": len(findings),
        }
        text = json.dumps(payload, indent=2, ensure_ascii=False) + "\n"

        def write(partials: list[Path]) -> None:
            partials[0].write_text(text, encoding="utf-8")

        publish_new_files([args.json_out], write)
    if findings and args.strict:
        return EXIT_GATE
    return EXIT_OK


if __name__ == "__main__":
    raise SystemExit(main())
