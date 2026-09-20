"""Shared path, YAML, and CSV helpers."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

# Human-facing review worksheets opened in Excel. Machine tables use UTF-8.
CSV_EXCEL = "utf-8-sig"
CSV_MACHINE = "utf-8"


def repo_root_from_script(script_file: str) -> Path:
    return Path(script_file).resolve().parents[1]


def ensure_repo_on_path(script_file: str) -> Path:
    root = repo_root_from_script(script_file)
    text = str(root)
    if text not in sys.path:
        sys.path.insert(0, text)
    return root


def require_new(path: Path) -> Path:
    path = Path(path)
    if path.exists():
        raise SystemExit(f"Refusing to overwrite: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def require_yaml():
    try:
        import yaml
    except ImportError as exc:  # pragma: no cover
        raise SystemExit("PyYAML required: pip install pyyaml") from exc
    return yaml


def load_yaml(path: Path) -> Any:
    yaml = require_yaml()
    return yaml.safe_load(Path(path).read_text(encoding="utf-8"))


def load_yaml_rows(path: Path) -> list[dict[str, Any]]:
    raw = load_yaml(path) or {}
    rows = raw.get("rows", raw) if isinstance(raw, dict) else raw
    if not isinstance(rows, list) or not rows:
        raise SystemExit("verdict table needs a non-empty 'rows' list")
    return rows


def exit_reason(exc: BaseException) -> str:
    if getattr(exc, "args", None):
        return str(exc.args[0])
    return str(exc)
