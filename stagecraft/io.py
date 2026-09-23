"""Shared path, YAML, and CSV helpers."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

# Human-facing review worksheets opened in Excel. Machine tables use UTF-8.
CSV_EXCEL = "utf-8-sig"
CSV_MACHINE = "utf-8"
IDENTITY_COLUMNS = (
    "dataset",
    "donor_id",
    "dataset_donor_id",
    "unit_id",
    "cell_id",
    "omitted_donor",
)
_BOOL_TRUE = {"true", "1"}
_BOOL_FALSE = {"false", "0"}


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


def parse_bool_column(series, name: str = "flag"):
    """Map true/false/1/0. Any other token, including yes/no, is a hard failure."""
    import pandas as pd

    if series.isna().any():
        raise SystemExit(f"{name} has missing values")
    if pd.api.types.is_bool_dtype(series):
        return series.astype(bool)
    mapped = series.astype(str).str.strip().str.lower().map(
        {**{token: True for token in _BOOL_TRUE}, **{token: False for token in _BOOL_FALSE}}
    )
    if mapped.isna().any():
        bad = sorted(set(series.astype(str)[mapped.isna()]))[:8]
        raise SystemExit(f"{name} must be true/false or 1/0; got {bad}")
    return mapped.astype(bool)


def read_identity_csv(path: Path, **kwargs):
    """Read a CSV without coercing donor or cell identifiers to numbers."""
    import pandas as pd

    path = Path(path)
    header = pd.read_csv(path, nrows=0)
    dtype = {column: str for column in IDENTITY_COLUMNS if column in header.columns}
    dtype.update(kwargs.pop("dtype", {}) or {})
    kwargs.setdefault("encoding", "utf-8-sig")
    return pd.read_csv(path, dtype=dtype, **kwargs)
