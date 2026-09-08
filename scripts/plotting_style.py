#!/usr/bin/env python3
"""Publication theme and PNG/PDF/SVG export with a parameter sidecar.

Used by Part 1 lock figures, the Part 2/4 TARGET_GENE catalogs, Part 3
compartment-round figures, Part 5 donor-unit association panels, and
Part 6 virtual-knockout Δaxis panels.
"""
from __future__ import annotations

import hashlib
import json
import platform
from datetime import datetime, timezone
from functools import lru_cache
from pathlib import Path
from textwrap import fill
from typing import Any, Mapping, Sequence

import matplotlib as mpl
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

try:
    import yaml
except ImportError:  # pragma: no cover
    yaml = None

PathLike = str | Path

DEFAULT_CONFIG = Path(__file__).with_name("plotting_config.yaml")


def load_plotting_config(path: PathLike | None = None) -> dict[str, Any]:
    cfg_path = Path(path or DEFAULT_CONFIG).resolve()
    if not cfg_path.is_file():
        raise FileNotFoundError(f"Plotting configuration not found: {cfg_path}")
    if yaml is None:
        raise SystemExit("PyYAML required: pip install pyyaml")
    with cfg_path.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    required = {"fonts", "figure_sizes_inches", "output", "theme", "umap"}
    missing = required.difference(config or {})
    if missing:
        raise ValueError(f"Plotting configuration lacks sections: {sorted(missing)}")
    return config


def load_color_map(path: PathLike) -> dict[str, str]:
    """Accept `{cell_type: {name: hex}}` or a flat `{name: hex}` YAML."""
    color_path = Path(path).resolve()
    if yaml is None:
        raise SystemExit("PyYAML required: pip install pyyaml")
    with color_path.open("r", encoding="utf-8") as handle:
        raw = yaml.safe_load(handle)
    if not isinstance(raw, dict):
        raise ValueError("Color YAML must be a mapping.")
    if "cell_type" in raw and isinstance(raw["cell_type"], dict):
        return {str(k): str(v) for k, v in raw["cell_type"].items()}
    if all(isinstance(v, str) and str(v).startswith("#") for v in raw.values()):
        return {str(k): str(v) for k, v in raw.items()}
    for key in ("cell_type_machine", "colors"):
        if key in raw and isinstance(raw[key], dict):
            return {str(k): str(v) for k, v in raw[key].items()}
    raise ValueError("Color YAML needs 'cell_type' or a flat name→#hex map.")


def select_font_stack(config: Mapping[str, Any]) -> list[str]:
    installed = {font.name for font in fm.fontManager.ttflist}
    candidates = list(config["fonts"]["latin"]) + list(config["fonts"]["cjk"])
    selected = list(dict.fromkeys(name for name in candidates if name in installed))
    if "DejaVu Sans" not in selected:
        selected.append("DejaVu Sans")
    return selected


def apply_publication_style(config: Mapping[str, Any]) -> list[str]:
    fonts = select_font_stack(config)
    theme = config["theme"]
    font_cfg = config["fonts"]
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": fonts,
            "font.size": theme["base_fontsize"],
            "axes.titlesize": theme["title_fontsize"],
            "axes.labelsize": theme["label_fontsize"],
            "xtick.labelsize": theme["tick_fontsize"],
            "ytick.labelsize": theme["tick_fontsize"],
            "legend.fontsize": theme["legend_fontsize"],
            "axes.linewidth": theme["axis_linewidth"],
            "axes.edgecolor": theme["foreground"],
            "axes.labelcolor": theme["foreground"],
            "text.color": theme["foreground"],
            "xtick.color": theme["foreground"],
            "ytick.color": theme["foreground"],
            "figure.facecolor": theme["background"],
            "axes.facecolor": theme["background"],
            "savefig.facecolor": theme["background"],
            "grid.color": theme["grid_color"],
            "grid.alpha": theme["grid_alpha"],
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": int(font_cfg["pdf_fonttype"]),
            "ps.fonttype": int(font_cfg["ps_fonttype"]),
            "svg.fonttype": font_cfg["svg_fonttype"],
            "mathtext.fontset": font_cfg["mathtext"],
        }
    )
    return fonts


def figure_size(config: Mapping[str, Any], preset: str) -> tuple[float, float]:
    width, height = config["figure_sizes_inches"][preset]
    return float(width), float(height)


def auto_point_style(n_cells: int, config: Mapping[str, Any]) -> tuple[float, float]:
    size_cfg = config["umap"]["point_size"]
    alpha_cfg = config["umap"]["alpha"]
    if n_cells < 5_000:
        key = "under_5000"
    elif n_cells < 20_000:
        key = "under_20000"
    elif n_cells < 75_000:
        key = "under_75000"
    else:
        key = "otherwise"
    return float(size_cfg[key]), float(alpha_cfg[key])


def wrap_labels(labels: Sequence[str], width: int = 28) -> list[str]:
    return [fill(str(label), width=width, break_long_words=False) for label in labels]


def _json_safe(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value.resolve())
    if isinstance(value, Mapping):
        return {str(k): _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_json_safe(v) for v in value]
    if hasattr(value, "item"):
        try:
            return value.item()
        except (ValueError, TypeError):
            pass
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return str(value)


@lru_cache(maxsize=64)
def file_sha256(path: PathLike, block_size: int = 2**20) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def save_figure(
    fig: Figure,
    output_stem: PathLike,
    config: Mapping[str, Any],
    *,
    parameters: Mapping[str, Any] | None = None,
    input_files: Sequence[PathLike] = (),
    formats: Sequence[str] | None = None,
    close: bool = True,
) -> dict[str, Path]:
    stem = Path(output_stem).resolve()
    stem.parent.mkdir(parents=True, exist_ok=True)
    output_cfg = config["output"]
    selected = list(formats or output_cfg["formats"])
    written: dict[str, Path] = {}
    for extension in selected:
        extension = extension.lower().lstrip(".")
        if extension not in {"png", "pdf", "svg"}:
            raise ValueError(f"Unsupported figure format: {extension}")
        path = stem.with_suffix(f".{extension}")
        kwargs: dict[str, Any] = {
            "bbox_inches": output_cfg["bbox_inches"],
            "pad_inches": output_cfg["pad_inches"],
            "transparent": bool(output_cfg["transparent"]),
        }
        if extension == "png":
            kwargs["dpi"] = int(output_cfg["png_dpi"])
        fig.savefig(path, **kwargs)
        written[extension] = path

    if output_cfg.get("write_parameter_sidecar", True):
        inputs = []
        for raw_path in input_files:
            path = Path(raw_path).resolve()
            inputs.append(
                {
                    "path": str(path),
                    "size_bytes": path.stat().st_size if path.exists() else None,
                    "sha256": file_sha256(path) if path.is_file() else None,
                }
            )
        sidecar = stem.with_suffix(".parameters.json")
        payload = {
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "python": platform.python_version(),
            "matplotlib": mpl.__version__,
            "outputs": {key: str(value) for key, value in written.items()},
            "parameters": _json_safe(parameters or {}),
            "inputs": inputs,
        }
        with sidecar.open("w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, ensure_ascii=False)
        written["parameters"] = sidecar
    if close:
        plt.close(fig)
    return written

