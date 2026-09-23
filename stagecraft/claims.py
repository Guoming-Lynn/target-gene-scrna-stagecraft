"""Review flags for sentences that cross the claim boundaries."""

from __future__ import annotations

import re
from dataclasses import dataclass
from importlib.resources import files
from pathlib import Path

_REQUIRED = ("id", "category", "pattern", "allowed", "example")
_SEPARATORS = set(".。;；!?！？")


@dataclass(frozen=True)
class Rule:
    id: str
    category: str
    pattern: re.Pattern[str]
    allowed: str
    example: str
    ignore_negation: bool = False


@dataclass(frozen=True)
class Finding:
    line: int
    column: int
    rule_id: str
    category: str
    matched: str
    allowed: str


def load_rules(path: Path | None = None) -> tuple[list[Rule], dict]:
    """Load claim rules. ``path`` overrides the packaged yaml."""
    import yaml

    if path is None:
        text = files("stagecraft").joinpath("claim_rules.yaml").read_text(encoding="utf-8")
    else:
        text = Path(path).read_text(encoding="utf-8")
    raw = yaml.safe_load(text) or {}
    if not isinstance(raw, dict) or not isinstance(raw.get("rules"), list) or not raw["rules"]:
        raise ValueError("claim rules need a non-empty rules list")
    settings = {
        "negation_window": int(raw.get("negation_window") or 40),
        "negations_en": list(raw.get("negations_en") or []),
        "negations_zh": list(raw.get("negations_zh") or []),
    }
    rules: list[Rule] = []
    seen: set[str] = set()
    for index, item in enumerate(raw["rules"], 1):
        if not isinstance(item, dict):
            raise ValueError(f"rule {index} is not a mapping")
        rule_id = str(item.get("id") or index)
        missing = [key for key in _REQUIRED if not str(item.get(key) or "").strip()]
        if missing:
            raise ValueError(f"rule {rule_id} missing {', '.join(missing)}")
        if rule_id in seen:
            raise ValueError(f"duplicate rule id {rule_id}")
        seen.add(rule_id)
        flags = 0 if item.get("case_sensitive") else re.IGNORECASE
        try:
            compiled = re.compile(str(item["pattern"]), flags)
        except re.error as exc:
            raise ValueError(f"rule {rule_id} has an invalid pattern: {exc}") from exc
        rules.append(
            Rule(
                id=rule_id,
                category=str(item["category"]),
                pattern=compiled,
                allowed=str(item["allowed"]),
                example=str(item["example"]),
                ignore_negation=bool(item.get("ignore_negation")),
            )
        )
    return rules, settings


def _mask_inline(line: str) -> str:
    chars = list(line)
    index = 0
    while index < len(chars):
        if chars[index] != "`":
            index += 1
            continue
        end = index + 1
        while end < len(chars) and chars[end] != "`":
            end += 1
        if end >= len(chars):
            break
        for cursor in range(index, end + 1):
            chars[cursor] = " "
        index = end + 1
    return "".join(chars)


def _negated(window: str, settings: dict) -> bool:
    last = -1
    for index, char in enumerate(window):
        if char in _SEPARATORS:
            last = index
    if last >= 0:
        window = window[last + 1 :]
    english = [re.escape(word) for word in settings.get("negations_en") or []]
    if english and re.search(r"\b(" + "|".join(english) + r")\b", window, re.IGNORECASE):
        return True
    if "n't" in window.lower():
        return True
    return any(token in window for token in settings.get("negations_zh") or [])


def lint_text(
    text: str,
    rules: list[Rule] | None = None,
    settings: dict | None = None,
) -> list[Finding]:
    """Return review findings. Code blocks, inline code, and allow-lines are ignored."""
    if rules is None or settings is None:
        loaded_rules, loaded_settings = load_rules()
        rules = loaded_rules if rules is None else rules
        settings = loaded_settings if settings is None else settings
    window_size = int(settings.get("negation_window") or 40)
    found: dict[tuple[int, int, str], Finding] = {}
    in_code = False
    for line_number, line in enumerate(text.splitlines(), 1):
        stripped = line.lstrip()
        if stripped.startswith("```") or stripped.startswith("~~~"):
            in_code = not in_code
            continue
        if in_code or "claim-lint: allow" in line:
            continue
        masked = _mask_inline(line)
        for rule in rules:
            for match in rule.pattern.finditer(masked):
                if not rule.ignore_negation:
                    start = match.start()
                    window = masked[max(0, start - window_size) : start]
                    if _negated(window, settings):
                        continue
                column = match.start() + 1
                key = (line_number, column, rule.id)
                found[key] = Finding(
                    line=line_number,
                    column=column,
                    rule_id=rule.id,
                    category=rule.category,
                    matched=match.group(0),
                    allowed=rule.allowed,
                )
    return [found[key] for key in sorted(found)]
