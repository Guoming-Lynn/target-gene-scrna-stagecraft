#!/usr/bin/env python3
"""Check that a frozen protocol has the required headings.

Usage:
    python scripts/check_protocol.py path/to/PROTOCOL.md
"""
from __future__ import annotations

import re
import sys
import argparse
from pathlib import Path
from protocol_chronology import check_chronology

REQUIRED_EN = [
    r"question|本阶段回答什么|0\.\s",
    r"hard boundar|硬边界|^#{1,3}\s+边界",
    r"input|输入|上游",
    r"unit|统计单位|eligibility|入选|生物学重复",
    r"model|模型|主模型",
    r"verdict|判定|措辞后果",
    r"not done|明确不做|明确不写|explicitly not",
]

REQUIRED_PHRASES = [
    ("evidence ceiling or 证据上限", r"evidence ceiling|证据上限|exploratory|formal"),
    ("forbidden sentences or 不能写", r"forbidden|不能写|do not write|正文不能写|不写因果|措辞后果|禁止写"),
    ("read-only upstream", r"read-only|只读"),
]


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("protocol", type=Path)
    parser.add_argument("--stage-root", type=Path, help="Verify local freeze hashes and output chronology")
    parser.add_argument("--freeze", action="store_true", help="Record a local freeze before stage artifacts exist")
    args = parser.parse_args(argv[1:])
    if args.freeze and args.stage_root is None:
        parser.error("--freeze requires --stage-root")
    path = args.protocol
    if not path.is_file():
        print(f"missing file: {path}", file=sys.stderr)
        return 1
    text = path.read_text(encoding="utf-8")
    if len(text.strip()) < 1200:
        print("FAIL", path)
        print(" - protocol is too short for a meaningful frozen protocol; keyword presence is not sufficient")
        return 1
    heads = re.findall(r"^#{1,3}\s+(.+)$", text, flags=re.M)
    blob = "\n".join("## " + head for head in heads)
    failures = []
    for i, pat in enumerate(REQUIRED_EN, 1):
        if not re.search(pat, blob, flags=re.I | re.M):
            failures.append(f"missing section pattern {i}: {pat}")
    for label, pat in REQUIRED_PHRASES:
        if not re.search(pat, text, flags=re.I | re.M):
            failures.append(f"missing {label}")
    if "NOT_ESTIMABLE" not in text and "not_estimable" not in text.lower():
        failures.append("verdict table should mention NOT_ESTIMABLE")
    # Require actual content after key headings, not headings/keywords alone.
    for label, pattern in (("question", r"(?im)^#{1,3}.*(question|本阶段回答什么)"),
                           ("input", r"(?im)^#{1,3}.*(input|输入|上游)"),
                           ("model", r"(?im)^#{1,3}.*(model|模型|主模型)"),
                           ("verdict", r"(?im)^#{1,3}.*(verdict|判定|措辞后果)")):
        match = re.search(pattern, text)
        if match:
            body = text[match.end():]
            body = re.split(r"\n#{1,3}\s+", body, maxsplit=1)[0]
            if len(re.sub(r"[`#|\-\s]", "", body)) < 80:
                failures.append(f"{label} heading lacks substantive content")
    if failures:
        print("FAIL", path)
        for f in failures:
            print(" -", f)
        return 1
    if args.stage_root:
        try:
            print(check_chronology(path, args.stage_root, args.freeze))
        except (ValueError, OSError, KeyError, TypeError) as exc:
            print("FAIL chronology:", exc)
            return 1
    else:
        print("CHRONOLOGY_NOT_VERIFIED: content checks only")
    print("OK", path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))

