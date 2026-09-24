"""Shared path, YAML, and CSV helpers."""

from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Any, Callable, Sequence

# UTF-8 with BOM. Part 5 tables are read by R with fileEncoding = "UTF-8-BOM",
# and the same encoding opens cleanly in Excel. CSV_MACHINE is UTF-8 without
# a BOM for text that is neither an R input nor a review sheet.
CSV_EXCEL = "utf-8-sig"
CSV_MACHINE = "utf-8"
IDENTITY_COLUMNS = (
    "dataset",
    "donor_id",
    "dataset_donor_id",
    "unit_id",
    "cell_id",
    "omitted_donor",
    "sample_id",
    "library_id",
    "source_block",
)
_BOOL_TRUE = {"true", "1"}
_BOOL_FALSE = {"false", "0"}


def repo_root_from_script(script_file: str) -> Path:
    return Path(script_file).resolve().parents[1]


def ensure_repo_on_path(script_file: str) -> Path:
    """Add the repository root and the script directory to ``sys.path``.

    The script directory stays ahead of the repository root so sibling
    modules resolve before an unrelated installed package of the same name.
    """
    script = Path(script_file).resolve()
    root = script.parents[1]
    for path in (root, script.parent):
        text = str(path)
        if text not in sys.path:
            sys.path.insert(0, text)
    return root


def require_new(path: Path) -> Path:
    """Reserve a new file. A concurrent writer loses instead of overwriting."""
    path = Path(path)
    if path.exists():
        raise SystemExit(f"Refusing to overwrite: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        fd = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        raise SystemExit(f"Refusing to overwrite: {path}") from None
    os.close(fd)
    return path


def _partial_path(path: Path) -> Path:
    """Sibling partial that keeps the original suffix.

    Matrix Market writers only emit a banner when the path ends in ``.mtx``.
    """
    if path.suffix:
        return path.with_name(path.stem + ".partial" + path.suffix)
    return path.with_name(path.name + ".partial")


def publish_new_files(paths: Sequence[Path], write: Callable[[list[Path]], None]) -> None:
    """Write every file to a sibling partial, then rename the set into place.

    Existence is checked before any final path is created. A failure deletes
    partials and any final already renamed by this call.
    """
    finals = [Path(path) for path in paths]
    if len(finals) != len(set(finals)):
        raise SystemExit("Refusing to publish duplicate output paths")
    partials = [_partial_path(path) for path in finals]
    for path in (*finals, *partials):
        if path.exists():
            raise SystemExit(f"Refusing to overwrite: {path}")
    for path in finals:
        path.parent.mkdir(parents=True, exist_ok=True)
    opened: list[Path] = []
    published: list[Path] = []
    try:
        for partial in partials:
            fd = os.open(partial, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            os.close(fd)
            opened.append(partial)
        write(partials)
        missing = [path.name for path in partials if not path.is_file()]
        if missing:
            raise SystemExit(f"Partial output was not written: {missing}")
        for partial, final in zip(partials, finals):
            os.rename(partial, final)
            opened.remove(partial)
            published.append(final)
    except BaseException:
        for partial in opened:
            partial.unlink(missing_ok=True)
        for final in published:
            final.unlink(missing_ok=True)
        raise


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
