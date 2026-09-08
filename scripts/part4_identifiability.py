#!/usr/bin/env python3
"""Part 4 identifiability forecast for a later donor-unit slope.

Does not fit limma. Does not reopen labels. Reads a unit×subtype table
(from part2/part4 figure source_data) and writes flags.

Usage:
    python scripts/part4_identifiability.py unit_summary.csv \\
        --group subtype --out tables/part4
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Any, Mapping

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from plotting_style import (  # noqa: E402
    apply_publication_style,
    figure_size,
    load_plotting_config,
    save_figure,
)

FLAG_COLORS = {
    "RANGE_OK": "#009E73",
    "EXPLORATORY_N": "#E69F00",
    "WITHIN_SOURCE_RANGE": "#56B4E9",
    "NO_EXPOSURE_RANGE": "#D55E00",
    "LOW_N": "#808080",
}

UNLIKELY_PATTERN = re.compile(
    r"unresolved|stressed|doublet|debris|contaminant|low[\s_-]?qc",
    re.I,
)


def _iqr(values: np.ndarray) -> float:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return float("nan")
    q75, q25 = np.percentile(finite, [75, 25])
    return float(q75 - q25)


def _spearman(x: np.ndarray, y: np.ndarray) -> float:
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 3:
        return float("nan")
    return float(pd.Series(x[mask]).corr(pd.Series(y[mask]), method="spearman"))


def flag_slice(
    n_units: int,
    n_sources: int,
    sd_detection: float,
    sd_log2cpm: float,
    *,
    n_formal: int = 12,
    n_exploratory: int = 8,
    sd_detection_min: float = 0.02,
    sd_log2cpm_min: float = 0.15,
    pooled: bool = False,
) -> str:
    if n_units < n_exploratory:
        return "LOW_N"
    det_dead = not np.isfinite(sd_detection) or sd_detection < sd_detection_min
    cpm_dead = not np.isfinite(sd_log2cpm) or sd_log2cpm < sd_log2cpm_min
    if det_dead and cpm_dead:
        return "NO_EXPOSURE_RANGE"
    if pooled and n_sources < 2:
        return "WITHIN_SOURCE_RANGE"
    if n_units < n_formal:
        return "EXPLORATORY_N"
    return "RANGE_OK"


def _notes(group: str, rho: float, collinear_min: float) -> str:
    bits = []
    if UNLIKELY_PATTERN.search(str(group) or ""):
        bits.append("UNLIKELY_PART5_ARM")
    if np.isfinite(rho) and abs(rho) >= collinear_min:
        bits.append("COLLINEAR_SCALES")
    return ";".join(bits)


def identifiability_table(
    units: pd.DataFrame,
    group_key: str,
    *,
    source_key: str = "dataset",
    donor_key: str = "dataset_donor_id",
    detection_col: str = "positive_fraction",
    cpm_col: str = "log2_cpm_plus_1",
    min_cells: int = 20,
    n_formal: int = 12,
    n_exploratory: int = 8,
    sd_detection_min: float = 0.02,
    sd_log2cpm_min: float = 0.15,
    collinear_min: float = 0.9,
) -> pd.DataFrame:
    frame = units.copy()
    required = {group_key, source_key, donor_key, detection_col, cpm_col, "n_cells"}
    if required - set(frame.columns):
        raise SystemExit(f"unit table missing columns: {sorted(required - set(frame.columns))}")
    for key in (group_key, source_key, donor_key):
        if frame[key].isna().any() or frame[key].astype(str).str.strip().eq("").any():
            raise SystemExit(f"Missing identity: {key}")
    if frame.duplicated([group_key, donor_key]).any():
        raise SystemExit("Require one biological donor row per subtype; aggregate technical repeats upstream")
    cell_counts = pd.to_numeric(frame["n_cells"], errors="coerce")
    if not np.isfinite(cell_counts).all() or (cell_counts < 0).any() or (cell_counts % 1 != 0).any():
        raise SystemExit("n_cells must contain nonnegative integers")
    frame["n_cells"] = cell_counts
    if "eligible" in frame.columns:
        flags = frame["eligible"].astype(str).str.strip().str.lower()
        if not flags.isin(["true", "false", "1", "0"]).all():
            raise SystemExit("eligible must contain explicit true/false or 1/0 values")
        eligible = frame.loc[flags.isin(["true", "1"]) & (cell_counts >= min_cells)].copy()
    else:
        eligible = frame.loc[frame["n_cells"] >= min_cells].copy()
    if source_key not in eligible.columns:
        raise SystemExit(f"unit table missing source column {source_key}")
    rows = []
    for key in (detection_col, cpm_col):
        values = pd.to_numeric(eligible[key], errors="coerce")
        if not np.isfinite(values).all():
            raise SystemExit(f"Non-finite eligible exposure: {key}")
        eligible[key] = values
    if not eligible[detection_col].between(0, 1).all():
        raise SystemExit("Detection fraction must be between zero and one")

    def summarize(subset: pd.DataFrame, group: str, source: str, n_sources: int) -> dict[str, Any]:
        det = subset[detection_col].to_numpy(float) if detection_col in subset else np.array([])
        cpm = subset[cpm_col].to_numpy(float) if cpm_col in subset else np.array([])
        rho = _spearman(det, cpm)
        sd_det = float(np.nanstd(det, ddof=1)) if len(det) >= 2 else float("nan")
        sd_cpm = float(np.nanstd(cpm, ddof=1)) if len(cpm) >= 2 else float("nan")
        flag = flag_slice(
            len(subset),
            n_sources,
            sd_det,
            sd_cpm,
            n_formal=n_formal,
            n_exploratory=n_exploratory,
            sd_detection_min=sd_detection_min,
            sd_log2cpm_min=sd_log2cpm_min,
            pooled=source == "ALL",
        )
        return {
            group_key: str(group),
            "source": str(source),
            "n_units_eligible": int(len(subset)),
            "n_cells": int(subset["n_cells"].sum()) if "n_cells" in subset else pd.NA,
            "n_sources_in_scope": int(n_sources),
            "detection_mean": float(np.nanmean(det)) if len(det) else float("nan"),
            "detection_sd": sd_det,
            "detection_iqr": _iqr(det),
            "detection_min": float(np.nanmin(det)) if len(det) else float("nan"),
            "detection_max": float(np.nanmax(det)) if len(det) else float("nan"),
            "log2cpm_mean": float(np.nanmean(cpm)) if len(cpm) else float("nan"),
            "log2cpm_sd": sd_cpm,
            "log2cpm_iqr": _iqr(cpm),
            "spearman_detection_log2cpm": rho,
            "flag": flag,
            "notes": _notes(group, rho, collinear_min),
            "source_key": source_key,
            "donor_key": donor_key,
            "source_independence": "NOT_ESTABLISHED_BY_FORECAST",
        }

    for group in frame[group_key].drop_duplicates():
        gframe = eligible.loc[eligible[group_key].eq(group)]
        sources = list(gframe[source_key].astype(str).unique())
        n_sources = len(sources)
        for source, subset in gframe.groupby(source_key, observed=True, sort=False):
            rows.append(summarize(subset, group, str(source), 1))
        rows.append(summarize(gframe, group, "ALL", n_sources))
    table = pd.DataFrame(rows)
    if table.empty:
        return table
    return table.sort_values([group_key, "source"])


def plot_identifiability(
    table: pd.DataFrame,
    group_key: str,
    config: Mapping[str, Any],
    output_stem: Path,
    *,
    gene: str = "",
    input_files: list[Path] | None = None,
) -> None:
    slice_table = table.loc[table["source"].astype(str).ne("ALL")].copy()
    if slice_table.empty:
        slice_table = table.copy()
    groups = list(dict.fromkeys(slice_table[group_key].astype(str)))
    sources = list(dict.fromkeys(slice_table["source"].astype(str)))
    grid = (
        slice_table.pivot_table(index=group_key, columns="source", values="flag", aggfunc="first")
        .reindex(index=groups, columns=sources)
    )
    code = {name: i for i, name in enumerate(FLAG_COLORS)}
    def _code(val: object) -> float:
        if pd.isna(val):
            return np.nan
        return float(code.get(str(val), np.nan))

    numeric = grid.apply(lambda series: series.map(_code))
    from matplotlib.colors import ListedColormap

    cmap = ListedColormap(list(FLAG_COLORS.values()))
    fig_w = max(figure_size(config, "wide_panel")[0], 0.55 * len(sources) + 3.2)
    fig_h = max(3.4, 0.38 * len(groups) + 1.8)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), layout="constrained")
    ax.imshow(
        numeric.to_numpy(float),
        aspect="auto",
        cmap=cmap,
        vmin=-0.5,
        vmax=len(FLAG_COLORS) - 0.5,
    )
    ax.set_xticks(range(len(sources)), sources, rotation=40, ha="right", fontsize=6.5)
    ax.set_yticks(range(len(groups)), groups, fontsize=7)
    ax.set_xlabel("Source")
    ax.set_ylabel("Locked subtype")
    title = "Part 5 identifiability forecast"
    if gene:
        title = f"{gene}: sd ≈ 0 cannot identify a Part 5 slope"
    ax.set_title(title, loc="left", fontweight="bold")
    counts = slice_table.pivot_table(
        index=group_key, columns="source", values="n_units_eligible", aggfunc="first"
    ).reindex(index=groups, columns=sources)
    for i, group in enumerate(groups):
        for j, source in enumerate(sources):
            n = counts.iat[i, j] if i < counts.shape[0] and j < counts.shape[1] else np.nan
            if pd.notna(n):
                ax.text(j, i, f"{int(n)}", ha="center", va="center", fontsize=6.2, color="#111111")
    handles = [Patch(facecolor=color, edgecolor="none", label=name) for name, color in FLAG_COLORS.items()]
    ax.legend(handles=handles, frameon=False, bbox_to_anchor=(1.02, 0.5), loc="center left", fontsize=6.5)
    save_figure(
        fig,
        output_stem,
        config,
        parameters={
            "plot": "identifiability_flag_heatmap",
            "cell_text": "n eligible donor-units",
            "color": "forecast flag, not a fitted verdict",
            "part5_not_run": True,
        },
        input_files=input_files or (),
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("unit_csv", type=Path)
    parser.add_argument("--group", default="subtype")
    parser.add_argument("--source-key", default="dataset")
    parser.add_argument("--donor-key", default="dataset_donor_id")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--gene", default="")
    parser.add_argument("--min-cells", type=int, default=20)
    parser.add_argument("--n-formal", type=int, default=12)
    parser.add_argument("--n-exploratory", type=int, default=8)
    parser.add_argument("--sd-detection-min", type=float, default=0.02)
    parser.add_argument("--sd-log2cpm-min", type=float, default=0.15)
    parser.add_argument("--config", default=str(SCRIPT_DIR / "plotting_config.yaml"))
    parser.add_argument("--no-plot", action="store_true")
    args = parser.parse_args(argv)

    units = pd.read_csv(args.unit_csv)
    table = identifiability_table(
        units,
        args.group,
        source_key=args.source_key,
        donor_key=args.donor_key,
        min_cells=args.min_cells,
        n_formal=args.n_formal,
        n_exploratory=args.n_exploratory,
        sd_detection_min=args.sd_detection_min,
        sd_log2cpm_min=args.sd_log2cpm_min,
    )
    args.out.mkdir(parents=True, exist_ok=True)
    csv_path = args.out / "identifiability.csv"
    if csv_path.exists():
        raise SystemExit(f"Refusing to overwrite: {csv_path}")
    table.to_csv(csv_path, index=False, encoding="utf-8-sig")
    if not args.no_plot:
        config = load_plotting_config(args.config)
        apply_publication_style(config)
        gene = args.gene or "TARGET_GENE"
        plot_identifiability(
            table,
            args.group,
            config,
            args.out / f"F04_15_{gene}_identifiability",
            gene=gene,
            input_files=[args.unit_csv],
        )
    print(f"wrote identifiability table ({len(table)} rows): {csv_path}")
    print("Flags are Part 5 forecasts, not fitted verdicts.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

