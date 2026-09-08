#!/usr/bin/env python3
"""Part 2 TARGET_GENE figure catalog.

Read-only: does not rename cells, recluster, or overwrite layers.
Detection = raw count > 0. Intensity displays use layers['normalized'].

Usage:
    python scripts/part2_figures.py atlas.h5ad \\
        --gene SYMBOL --group cell_type --colors colors.yaml \\
        --out figures/part2 --min-cells 20
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any, Mapping, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy import sparse

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from plotting_style import (  # noqa: E402
    apply_publication_style,
    auto_point_style,
    figure_size,
    load_color_map,
    load_plotting_config,
    save_figure,
    wrap_labels,
)

try:
    import anndata as sc
except ImportError:  # pragma: no cover
    sc = None


def gene_vector(adata, gene: str, layer: str) -> np.ndarray:
    if gene not in adata.var_names:
        raise KeyError(f"TARGET_GENE '{gene}' is not in var_names. Stop.")
    if layer not in adata.layers:
        raise KeyError(f"Missing layers['{layer}'].")
    matrix = adata[:, gene].layers[layer]
    if sparse.issparse(matrix):
        return np.asarray(matrix.toarray()).ravel()
    return np.asarray(matrix).ravel()


def _coords(adata, key: str = "X_umap") -> np.ndarray:
    if key not in adata.obsm:
        raise KeyError(f"Missing obsm['{key}']. Do not recompute UMAP in Part 2.")
    xy = np.asarray(adata.obsm[key])[:, :2]
    return xy


def _color(colors: Mapping[str, str], name: str, fallback: str = "#808080") -> str:
    return str(colors.get(name, fallback))


def summarize_groups(
    adata,
    gene: str,
    group_key: str,
    *,
    donor_key: str = "donor_id",
    dataset_key: str = "dataset",
    counts_layer: str = "counts",
    normalized_layer: str = "normalized",
) -> pd.DataFrame:
    frame = adata.obs[[group_key, donor_key, dataset_key]].copy()
    frame["raw_count"] = gene_vector(adata, gene, counts_layer)
    frame["normalized"] = gene_vector(adata, gene, normalized_layer)
    frame["positive"] = frame["raw_count"] > 0
    total = float(frame["raw_count"].sum())
    n_cells_all = len(frame)
    rows = []
    for group, subset in frame.groupby(group_key, observed=True, sort=False):
        pos = subset.loc[subset["positive"], "normalized"]
        raw_sum = float(subset["raw_count"].sum())
        rows.append(
            {
                group_key: str(group),
                "n_cells": int(len(subset)),
                "cell_fraction": len(subset) / n_cells_all if n_cells_all else np.nan,
                "n_donors": int(subset[[dataset_key, donor_key]].drop_duplicates().shape[0]),
                "n_datasets": int(subset[dataset_key].nunique()),
                "n_positive": int(subset["positive"].sum()),
                "positive_fraction": float(subset["positive"].mean()),
                "mean_normalized_all": float(subset["normalized"].mean()),
                "mean_normalized_positive": float(pos.mean()) if len(pos) else 0.0,
                "raw_count_sum": raw_sum,
                "count_contribution": raw_sum / total if total > 0 else np.nan,
            }
        )
    return pd.DataFrame(rows)


def summarize_units(
    adata,
    gene: str,
    group_key: str,
    *,
    donor_key: str = "donor_id",
    dataset_key: str = "dataset",
    unit_key: str = "dataset_donor_id",
    counts_layer: str = "counts",
    normalized_layer: str = "normalized",
    min_cells: int = 20,
) -> pd.DataFrame:
    obs = adata.obs
    if unit_key not in obs.columns:
        units = obs[dataset_key].astype(str) + "_" + obs[donor_key].astype(str)
    else:
        units = obs[unit_key].astype(str)
    frame = pd.DataFrame(
        {
            group_key: obs[group_key].astype(str).to_numpy(),
            dataset_key: obs[dataset_key].astype(str).to_numpy(),
            donor_key: obs[donor_key].astype(str).to_numpy(),
            unit_key: units.to_numpy() if hasattr(units, "to_numpy") else np.asarray(units),
            "raw_count": gene_vector(adata, gene, counts_layer),
            "normalized": gene_vector(adata, gene, normalized_layer),
            "total_umi": np.asarray(
                adata.layers[counts_layer].sum(axis=1)
            ).ravel(),
        }
    )
    frame["positive"] = frame["raw_count"] > 0
    rows = []
    for (unit, group), subset in frame.groupby([unit_key, group_key], observed=True, sort=False):
        n = len(subset)
        umi = float(subset["total_umi"].sum())
        gene_sum = float(subset["raw_count"].sum())
        cpm = (gene_sum / umi * 1e6) if umi > 0 else np.nan
        pos = subset.loc[subset["positive"], "normalized"]
        rows.append(
            {
                unit_key: str(unit),
                dataset_key: str(subset[dataset_key].iloc[0]),
                donor_key: str(subset[donor_key].iloc[0]),
                group_key: str(group),
                "n_cells": int(n),
                "eligible": n >= min_cells,
                "positive_fraction": float(subset["positive"].mean()),
                "mean_normalized_all": float(subset["normalized"].mean()),
                "mean_normalized_positive": float(pos.mean()) if len(pos) else np.nan,
                "cpm": cpm,
                "log2_cpm_plus_1": np.log2(cpm + 1) if np.isfinite(cpm) else np.nan,
                "mean_umi": float(subset["total_umi"].mean()),
            }
        )
    return pd.DataFrame(rows)


def leave_one_dataset_out(
    group_summary: pd.DataFrame,
    unit_table: pd.DataFrame,
    group_key: str,
    dataset_key: str = "dataset",
    metric: str = "positive_fraction",
) -> pd.DataFrame:
    """Recompute group metric after dropping each dataset. Descriptive only."""
    datasets = sorted(unit_table[dataset_key].astype(str).unique())
    rows = []
    for held in datasets:
        remain = unit_table.loc[unit_table[dataset_key].astype(str).ne(held)]
        if remain.empty:
            continue
        # weight unit means by n_cells
        for group, subset in remain.groupby(group_key, observed=True):
            weights = subset["n_cells"].to_numpy(float)
            values = subset[metric].to_numpy(float)
            mask = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
            if not mask.any():
                value = np.nan
            else:
                value = float(np.average(values[mask], weights=weights[mask]))
            rows.append(
                {
                    "held_out_dataset": held,
                    group_key: str(group),
                    metric: value,
                    "n_units_remaining": int(subset.shape[0]),
                }
            )
    return pd.DataFrame(rows)


def plot_umap_categorical(
    adata,
    group_key: str,
    colors: Mapping[str, str],
    config: Mapping[str, Any],
    *,
    order: Sequence[str] | None = None,
    output_stem: Path | None = None,
    title: str = "Locked global cell types",
) -> None:
    xy = _coords(adata)
    groups = adata.obs[group_key].astype(str)
    categories = list(order) if order is not None else list(dict.fromkeys(groups))
    fig, ax = plt.subplots(figsize=figure_size(config, "double_panel"), layout="constrained")
    size, alpha = auto_point_style(adata.n_obs, config)
    handles = []
    for category in reversed(categories):
        mask = groups.eq(category).to_numpy()
        ax.scatter(
            xy[mask, 0],
            xy[mask, 1],
            s=size,
            alpha=alpha,
            color=_color(colors, category),
            linewidths=0,
            rasterized=True,
        )
    for category in categories:
        handles.append(
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="",
                markersize=4,
                color=_color(colors, category),
                label=category,
            )
        )
    ax.legend(
        handles=handles,
        title=group_key,
        bbox_to_anchor=(1.02, 0.5),
        loc="center left",
        frameon=False,
        fontsize=6.5,
        labelspacing=0.35,
    )
    ax.set_title(title, loc="left", fontweight="bold")
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_xticks([])
    ax.set_yticks([])
    if output_stem is not None:
        save_figure(
            fig,
            output_stem,
            config,
            parameters={
                "plot": "umap_categorical",
                "group_key": group_key,
                "n_cells": int(adata.n_obs),
                "order_rule": "legend = cell-count descending; draw small groups last",
            },
        )


def plot_umap_feature(
    adata,
    gene: str,
    config: Mapping[str, Any],
    *,
    layer: str = "normalized",
    clip_percentile: float | None = 99,
    output_stem: Path | None = None,
    title: str | None = None,
) -> None:
    xy = _coords(adata)
    values = gene_vector(adata, gene, layer)
    display = values.copy()
    vmax = float(np.max(display)) if len(display) else 0.0
    if clip_percentile is not None:
        positive = display[display > 0]
        vmax = float(np.percentile(positive, clip_percentile)) if len(positive) else 0.0
        display = np.minimum(display, vmax)
    fig, ax = plt.subplots(figsize=(4.15, 3.65), layout="constrained")
    size, alpha = auto_point_style(adata.n_obs, config)
    order = np.argsort(display, kind="stable")
    artist = ax.scatter(
        xy[order, 0],
        xy[order, 1],
        c=display[order],
        cmap=config["theme"]["expression_continuous"],
        vmin=0,
        vmax=vmax if vmax > 0 else 1,
        s=size,
        alpha=alpha,
        linewidths=0,
        rasterized=True,
    )
    cbar = fig.colorbar(artist, ax=ax, orientation="horizontal", fraction=0.055, pad=0.08, aspect=35)
    suffix = f" (display clipped at p{clip_percentile:g})" if clip_percentile else ""
    cbar.set_label(f"Log-normalized {gene}{suffix}")
    ax.set_title(title or f"{gene} expression", loc="left", fontweight="bold")
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_xticks([])
    ax.set_yticks([])
    if output_stem is not None:
        save_figure(
            fig,
            output_stem,
            config,
            parameters={
                "plot": "umap_feature",
                "gene": gene,
                "layer": layer,
                "clip_percentile": clip_percentile,
                "display_vmax": vmax,
                "zeros_retained": True,
            },
        )


def plot_umap_binary(
    adata,
    gene: str,
    config: Mapping[str, Any],
    *,
    layer: str = "counts",
    output_stem: Path | None = None,
    title: str | None = None,
) -> None:
    xy = _coords(adata)
    positive = gene_vector(adata, gene, layer) > 0
    fig, ax = plt.subplots(figsize=figure_size(config, "single_panel"), layout="constrained")
    size, _ = auto_point_style(adata.n_obs, config)
    ax.scatter(
        xy[~positive, 0],
        xy[~positive, 1],
        s=size,
        color=config["theme"].get("undetected", "#D9D9D9"),
        alpha=0.28,
        linewidths=0,
        rasterized=True,
        label=f"Not detected (n={int((~positive).sum()):,})",
    )
    ax.scatter(
        xy[positive, 0],
        xy[positive, 1],
        s=max(size, 1.1),
        color=config["theme"].get("detected", "#D73027"),
        alpha=0.75,
        linewidths=0,
        rasterized=True,
        label=f"Detected (n={int(positive.sum()):,})",
    )
    ax.set_title(title or f"{gene} raw-count detection", loc="left", fontweight="bold")
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.legend(frameon=False, loc="best", fontsize=7)
    if output_stem is not None:
        save_figure(
            fig,
            output_stem,
            config,
            parameters={
                "plot": "umap_binary",
                "gene": gene,
                "layer": layer,
                "positive_definition": "raw count > 0",
                "n_positive": int(positive.sum()),
            },
        )


def plot_umap_pair(
    adata,
    gene: str,
    group_key: str,
    colors: Mapping[str, str],
    config: Mapping[str, Any],
    *,
    order: Sequence[str],
    output_stem: Path | None = None,
) -> None:
    xy = _coords(adata)
    groups = adata.obs[group_key].astype(str)
    values = gene_vector(adata, gene, "normalized")
    positive = values[values > 0]
    vmax = float(np.percentile(positive, 99)) if len(positive) else 1.0
    display = np.minimum(values, vmax)
    size, alpha = auto_point_style(adata.n_obs, config)
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.7), layout="constrained")
    for category in reversed(list(order)):
        mask = groups.eq(category).to_numpy()
        axes[0].scatter(
            xy[mask, 0],
            xy[mask, 1],
            s=size,
            alpha=alpha,
            color=_color(colors, category),
            linewidths=0,
            rasterized=True,
        )
    axes[0].set_title("a  Locked cell types", loc="left", fontweight="bold")
    sort = np.argsort(display, kind="stable")
    artist = axes[1].scatter(
        xy[sort, 0],
        xy[sort, 1],
        c=display[sort],
        cmap=config["theme"]["expression_continuous"],
        vmin=0,
        vmax=vmax if vmax > 0 else 1,
        s=size,
        alpha=alpha,
        linewidths=0,
        rasterized=True,
    )
    fig.colorbar(artist, ax=axes[1], orientation="horizontal", fraction=0.05, pad=0.08, aspect=30)
    axes[1].set_title(f"b  {gene} intensity", loc="left", fontweight="bold")
    for ax in axes:
        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        ax.set_xticks([])
        ax.set_yticks([])
    if output_stem is not None:
        save_figure(
            fig,
            output_stem,
            config,
            parameters={"plot": "umap_pair", "gene": gene, "same_coordinates": True},
        )


def plot_lollipop(
    summary: pd.DataFrame,
    group_key: str,
    gene: str,
    colors: Mapping[str, str],
    config: Mapping[str, Any],
    *,
    output_stem: Path | None = None,
) -> None:
    ordered = summary.sort_values("mean_normalized_all", ascending=True).reset_index(drop=True)
    height = max(3.6, 0.33 * len(ordered) + 0.9)
    fig, ax = plt.subplots(figsize=(6.4, height), layout="constrained")
    y = np.arange(len(ordered))
    means = ordered["mean_normalized_all"].to_numpy(float)
    fractions = ordered["positive_fraction"].to_numpy(float)
    ax.hlines(y, 0, means, color="#D9D9D9", linewidth=0.9, zorder=1)
    ax.scatter(
        means,
        y,
        s=36 + fractions * 680,
        color=[_color(colors, str(v)) for v in ordered[group_key]],
        edgecolor="white",
        linewidth=0.7,
        zorder=2,
    )
    ax.set_yticks(y)
    ax.set_yticklabels(wrap_labels(ordered[group_key].astype(str), 28))
    ax.set_xlim(0, max(0.1, float(np.nanmax(means)) * 1.23))
    ax.set_xlabel(f"Mean log-normalized {gene} (zeros retained)")
    ax.set_title(f"{gene} detection and intensity", loc="left", fontweight="bold")
    handles = [
        Line2D(
            [],
            [],
            marker="o",
            linestyle="None",
            markersize=float(np.sqrt(30 + value * 520)),
            markerfacecolor="white",
            markeredgecolor="#333333",
        )
        for value in (0.1, 0.3, 0.5)
    ]
    ax.legend(handles, ["10%", "30%", "50%"], title="Detected fraction", frameon=False, loc="lower right")
    for y_pos, donors in zip(y, ordered["n_donors"]):
        ax.text(1.01, y_pos, f"n={int(donors)}", transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=6.3, color="#555555")
    if output_stem is not None:
        save_figure(
            fig,
            output_stem,
            config,
            parameters={
                "plot": "lollipop",
                "dot_size": "positive_fraction",
                "dot_position": "mean_normalized_all",
                "not": "a one-column Scanpy dotplot",
            },
        )


def plot_violin_box(
    adata,
    gene: str,
    group_key: str,
    colors: Mapping[str, str],
    config: Mapping[str, Any],
    *,
    positive_only: bool,
    output_stem: Path | None = None,
) -> None:
    values = gene_vector(adata, gene, "normalized")
    counts = gene_vector(adata, gene, "counts")
    groups = adata.obs[group_key].astype(str).to_numpy()
    mask = np.isfinite(values)
    if positive_only:
        mask &= counts > 0
    frame = pd.DataFrame({group_key: groups[mask], "value": values[mask]})
    medians = frame.groupby(group_key, observed=True)["value"].median().sort_values(ascending=False, kind="stable")
    order = medians.index.astype(str).tolist()
    arrays = [frame.loc[frame[group_key].eq(g), "value"].to_numpy(float) for g in order]
    positions = np.arange(len(order))
    fig, ax = plt.subplots(figsize=(7.4, max(3.7, 0.53 * len(order) + 1.2)), layout="constrained")
    violin = ax.violinplot(arrays, positions=positions, vert=False, showmeans=False, showmedians=False, showextrema=False, widths=0.78)
    for body, group in zip(violin["bodies"], order):
        body.set_facecolor(_color(colors, group))
        body.set_edgecolor("none")
        body.set_alpha(0.42)
    box = ax.boxplot(
        arrays,
        positions=positions,
        vert=False,
        widths=0.24,
        showfliers=False,
        patch_artist=True,
        whis=(5, 95),
        medianprops={"color": config["theme"].get("median", "#B2182B"), "linewidth": 1.1},
        whiskerprops={"color": "#555555", "linewidth": 0.75},
        capprops={"color": "#555555", "linewidth": 0.75},
    )
    for patch, group in zip(box["boxes"], order):
        patch.set_facecolor(_color(colors, group))
        patch.set_alpha(0.18)
        patch.set_edgecolor(_color(colors, group))
    ax.set_yticks(positions)
    ax.set_yticklabels(wrap_labels(order, 26))
    ax.invert_yaxis()
    ax.set_xlabel(f"Log-normalized {gene}")
    title = f"{gene} among detected cells" if positive_only else f"{gene} across all cells, including zeros"
    ax.set_title(title, loc="left", fontweight="bold")
    ax.grid(axis="x")
    for position, arr in enumerate(arrays):
        ax.text(1.005, position, f"cells={len(arr):,}", transform=ax.get_yaxis_transform(), va="center", fontsize=6.2, color="#555555")
    if output_stem is not None:
        save_figure(
            fig,
            output_stem,
            config,
            parameters={
                "plot": "violin_box",
                "positive_only": positive_only,
                "order_rule": "descending median of this panel",
                "box_whiskers": "5th–95th percentiles",
                "cells_are_not_n": True,
            },
        )


def plot_detection_vs_intensity(
    summary: pd.DataFrame,
    group_key: str,
    gene: str,
    colors: Mapping[str, str],
    config: Mapping[str, Any],
    *,
    output_stem: Path | None = None,
) -> None:
    fig, ax = plt.subplots(figsize=(7.3, 5.1), layout="constrained")
    for _, row in summary.iterrows():
        name = str(row[group_key])
        ax.scatter(
            row["positive_fraction"],
            row["mean_normalized_positive"],
            s=35 + 0.035 * float(row["n_cells"]),
            color=_color(colors, name),
            edgecolor="white",
            linewidth=0.8,
            zorder=3,
        )
        ax.annotate(name, (row["positive_fraction"], row["mean_normalized_positive"]), xytext=(6, 5), textcoords="offset points", fontsize=6.2)
    ax.set_xlabel(f"{gene} detected fraction (raw count > 0)")
    ax.set_ylabel("Mean normalized expression among detected cells")
    ax.set_title("Detection and positive-cell intensity are distinct", loc="left", fontweight="bold")
    ax.grid(alpha=0.25)
    if output_stem is not None:
        save_figure(
            fig,
            output_stem,
            config,
            parameters={"plot": "detection_vs_intensity", "size": "n_cells"},
        )


def plot_composition_contribution(
    summary: pd.DataFrame,
    group_key: str,
    gene: str,
    colors: Mapping[str, str],
    config: Mapping[str, Any],
    *,
    output_stem: Path | None = None,
) -> None:
    order = summary.sort_values("count_contribution", ascending=False)[group_key].astype(str).tolist()
    data = summary.set_index(group_key).loc[order].reset_index()
    y = np.arange(len(data))
    height = 0.32
    fig, ax = plt.subplots(figsize=(7.5, max(4.0, 0.48 * len(data) + 1.0)), layout="constrained")
    ax.barh(y - height / 2, data["cell_fraction"], height=height, color="#C8CDD3", label="Cell fraction")
    ax.barh(
        y + height / 2,
        data["count_contribution"],
        height=height,
        color=[_color(colors, str(v)) for v in data[group_key]],
        label=f"{gene} raw-count contribution",
    )
    ax.set_yticks(y)
    ax.set_yticklabels(wrap_labels(data[group_key].astype(str), 28))
    ax.invert_yaxis()
    ax.set_xlabel("Fraction of cells or of all TARGET_GENE UMIs")
    ax.set_title("Composition vs count contribution", loc="left", fontweight="bold")
    ax.legend(frameon=False, loc="lower right")
    ax.grid(axis="x")
    if output_stem is not None:
        save_figure(
            fig,
            output_stem,
            config,
            parameters={"plot": "composition_vs_contribution", "order_rule": "descending count contribution"},
        )


def _box_strip(ax, data: pd.DataFrame, metric: str, group_key: str, colors: Mapping[str, str], title: str, xlabel: str, seed: int, median_color: str) -> None:
    order = (
        data.groupby(group_key, observed=True)[metric]
        .median()
        .sort_values(ascending=False, kind="stable")
        .index.astype(str)
        .tolist()
    )
    arrays = [data.loc[data[group_key].astype(str).eq(g), metric].dropna().to_numpy(float) for g in order]
    positions = np.arange(len(order))
    box = ax.boxplot(
        arrays,
        positions=positions,
        widths=0.55,
        vert=False,
        showfliers=False,
        patch_artist=True,
        medianprops={"color": median_color, "linewidth": 1.0},
        whiskerprops={"color": "#555555", "linewidth": 0.7},
        capprops={"color": "#555555", "linewidth": 0.7},
    )
    for patch, group in zip(box["boxes"], order):
        patch.set_facecolor(_color(colors, group))
        patch.set_alpha(0.23)
        patch.set_edgecolor(_color(colors, group))
    rng = np.random.default_rng(seed)
    for position, group in enumerate(order):
        subset = data.loc[data[group_key].astype(str).eq(group), metric].dropna()
        jitter = rng.uniform(-0.12, 0.12, len(subset))
        ax.scatter(subset, position + jitter, s=15, color=_color(colors, group), edgecolor="white", linewidth=0.35, alpha=0.9, zorder=3)
        ax.text(1.01, position, f"n={len(subset)}", transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=5.8, color="#555555")
    ax.set_yticks(positions)
    ax.set_yticklabels(wrap_labels(order, 22))
    ax.invert_yaxis()
    ax.set_xlabel(xlabel)
    ax.set_title(title, loc="left", fontweight="bold")
    ax.grid(axis="x")


def plot_donor_endpoint_grid(
    units: pd.DataFrame,
    group_key: str,
    gene: str,
    colors: Mapping[str, str],
    config: Mapping[str, Any],
    *,
    min_cells: int,
    output_stem: Path | None = None,
    seed: int = 42,
) -> None:
    eligible = units.loc[units["n_cells"] >= min_cells].copy()
    if eligible.empty:
        raise ValueError(f"No unit × {group_key} rows with ≥{min_cells} cells.")
    metrics = [
        ("positive_fraction", "a  Detection fraction", "Fraction detected"),
        ("mean_normalized_positive", "b  Positive-cell intensity", "Mean normalized expression"),
        ("mean_normalized_all", "c  All-cell mean", "Mean normalized expression"),
        ("log2_cpm_plus_1", "d  Pseudobulk burden", f"log2({gene} CPM + 1)"),
    ]
    n_groups = eligible[group_key].astype(str).nunique()
    fig, axes = plt.subplots(2, 2, figsize=(8.4, max(6.4, 0.78 * n_groups + 1.0)), layout="constrained")
    median = config["theme"].get("median", "#B2182B")
    for i, (metric, title, xlabel) in enumerate(metrics):
        _box_strip(axes.flat[i], eligible, metric, group_key, colors, title, xlabel, seed + i, median)
    fig.suptitle(
        f"{gene} across dataset × donor units (≥{min_cells} cells per unit × type)",
        fontsize=11,
        fontweight="bold",
    )
    if output_stem is not None:
        save_figure(
            fig,
            output_stem,
            config,
            parameters={
                "plot": "donor_endpoint_grid",
                "biological_unit": "dataset × donor_id",
                "one_point": "one eligible unit × cell type",
                "panel_order": "each endpoint independently by median",
                "minimum_cells": min_cells,
            },
        )


def plot_unit_heatmap(
    units: pd.DataFrame,
    group_key: str,
    unit_key: str,
    config: Mapping[str, Any],
    *,
    metric: str = "positive_fraction",
    output_stem: Path | None = None,
) -> None:
    pivot = units.pivot_table(index=unit_key, columns=group_key, values=metric, aggfunc="mean")
    fig, ax = plt.subplots(figsize=(max(6.2, 0.42 * pivot.shape[1] + 3.2), max(4.2, 0.22 * pivot.shape[0] + 1.4)), layout="constrained")
    im = ax.imshow(pivot.to_numpy(float), aspect="auto", cmap="magma", vmin=0, interpolation="nearest")
    ax.set_xticks(np.arange(pivot.shape[1]))
    ax.set_xticklabels(wrap_labels(pivot.columns.astype(str), 16), rotation=45, ha="right", fontsize=6.5)
    ax.set_yticks(np.arange(pivot.shape[0]))
    ax.set_yticklabels(pivot.index.astype(str), fontsize=5.8)
    ax.set_title(f"Unit × type {metric}", loc="left", fontweight="bold")
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    if output_stem is not None:
        save_figure(fig, output_stem, config, parameters={"plot": "unit_heatmap", "metric": metric, "missing": "NaN = type absent in that unit"})


def plot_lodo_heatmap(
    lodo: pd.DataFrame,
    group_key: str,
    config: Mapping[str, Any],
    *,
    metric: str = "positive_fraction",
    output_stem: Path | None = None,
) -> None:
    pivot = lodo.pivot_table(index="held_out_dataset", columns=group_key, values=metric, aggfunc="mean")
    fig, ax = plt.subplots(figsize=(max(6.0, 0.45 * pivot.shape[1] + 2.8), max(2.8, 0.45 * pivot.shape[0] + 1.2)), layout="constrained")
    im = ax.imshow(pivot.to_numpy(float), aspect="auto", cmap="magma", vmin=0, interpolation="nearest")
    ax.set_xticks(np.arange(pivot.shape[1]))
    ax.set_xticklabels(wrap_labels(pivot.columns.astype(str), 16), rotation=45, ha="right", fontsize=6.5)
    ax.set_yticks(np.arange(pivot.shape[0]))
    ax.set_yticklabels(pivot.index.astype(str), fontsize=7)
    ax.set_ylabel("Held-out dataset")
    ax.set_title("Leave-one-dataset-out (descriptive)", loc="left", fontweight="bold")
    fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    if output_stem is not None:
        save_figure(
            fig,
            output_stem,
            config,
            parameters={"plot": "lodo_heatmap", "claim": "description, not independent replication"},
        )


def plot_depth(
    units: pd.DataFrame,
    group_key: str,
    gene: str,
    colors: Mapping[str, str],
    config: Mapping[str, Any],
    *,
    output_stem: Path | None = None,
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.6), layout="constrained")
    for _, row in units.iterrows():
        c = _color(colors, str(row[group_key]))
        axes[0].scatter(np.log1p(row["mean_umi"]), row["positive_fraction"], s=18, color=c, alpha=0.85, edgecolor="white", linewidth=0.3)
        axes[1].scatter(row["n_cells"], row["positive_fraction"], s=18, color=c, alpha=0.85, edgecolor="white", linewidth=0.3)
    axes[0].set_xlabel("log1p mean UMI")
    axes[1].set_xlabel("Cells in unit × type")
    for ax in axes:
        ax.set_ylabel(f"{gene} detection fraction")
        ax.grid(alpha=0.25)
    axes[0].set_title("a  Depth", loc="left", fontweight="bold")
    axes[1].set_title("b  Cell number", loc="left", fontweight="bold")
    if output_stem is not None:
        save_figure(fig, output_stem, config, parameters={"plot": "depth_diagnostics", "unit": "dataset × donor × type"})


def render_catalog(
    adata,
    *,
    gene: str,
    group_key: str,
    colors: Mapping[str, str],
    config: Mapping[str, Any],
    out_dir: Path,
    min_cells: int = 20,
    donor_key: str = "donor_id",
    dataset_key: str = "dataset",
    stem_prefix: str = "F02",
    umap_title: str = "Locked global cell types",
) -> dict[str, Path]:
    apply_publication_style(config)
    out_dir.mkdir(parents=True, exist_ok=True)
    tables = out_dir / "source_data"
    tables.mkdir(exist_ok=True)
    group_summary = summarize_groups(adata, gene, group_key, donor_key=donor_key, dataset_key=dataset_key)
    units = summarize_units(adata, gene, group_key, donor_key=donor_key, dataset_key=dataset_key, min_cells=min_cells)
    group_csv = tables / "group_summary.csv"
    unit_csv = tables / "unit_summary.csv"
    group_summary.to_csv(group_csv, index=False)
    units.to_csv(unit_csv, index=False)
    cell_count_order = group_summary.sort_values("n_cells", ascending=False)[group_key].astype(str).tolist()
    stem = {
        "01": out_dir / f"{stem_prefix}_01_{gene}_celltype_UMAP",
        "02": out_dir / f"{stem_prefix}_02_{gene}_feature_full",
        "03": out_dir / f"{stem_prefix}_03_{gene}_feature_p99",
        "04": out_dir / f"{stem_prefix}_04_{gene}_detection_UMAP",
        "05": out_dir / f"{stem_prefix}_05_{gene}_celltype_and_expression",
        "06": out_dir / f"{stem_prefix}_06_{gene}_lollipop",
        "07": out_dir / f"{stem_prefix}_07_{gene}_violin_all",
        "08": out_dir / f"{stem_prefix}_08_{gene}_violin_positive",
        "09": out_dir / f"{stem_prefix}_09_{gene}_detection_vs_intensity",
        "10": out_dir / f"{stem_prefix}_10_{gene}_composition_vs_contribution",
        "11": out_dir / f"{stem_prefix}_11_{gene}_donor_endpoints",
        "12": out_dir / f"{stem_prefix}_12_{gene}_unit_heatmap",
        "13": out_dir / f"{stem_prefix}_13_{gene}_lodo",
        "14": out_dir / f"{stem_prefix}_14_{gene}_depth",
    }
    plot_umap_categorical(
        adata,
        group_key,
        colors,
        config,
        order=cell_count_order,
        output_stem=stem["01"],
        title=umap_title,
    )
    plot_umap_feature(adata, gene, config, clip_percentile=None, output_stem=stem["02"], title=f"{gene} expression (full scale)")
    plot_umap_feature(adata, gene, config, clip_percentile=99, output_stem=stem["03"], title=f"{gene} expression (positive-cell p99 clip)")
    plot_umap_binary(adata, gene, config, output_stem=stem["04"])
    plot_umap_pair(adata, gene, group_key, colors, config, order=cell_count_order, output_stem=stem["05"])
    plot_lollipop(group_summary, group_key, gene, colors, config, output_stem=stem["06"])
    plot_violin_box(adata, gene, group_key, colors, config, positive_only=False, output_stem=stem["07"])
    plot_violin_box(adata, gene, group_key, colors, config, positive_only=True, output_stem=stem["08"])
    plot_detection_vs_intensity(group_summary, group_key, gene, colors, config, output_stem=stem["09"])
    plot_composition_contribution(group_summary, group_key, gene, colors, config, output_stem=stem["10"])
    plot_donor_endpoint_grid(units, group_key, gene, colors, config, min_cells=min_cells, output_stem=stem["11"])
    unit_col = "dataset_donor_id" if "dataset_donor_id" in units.columns else str(units.columns[0])
    plot_unit_heatmap(units, group_key, unit_col, config, output_stem=stem["12"])
    lodo = leave_one_dataset_out(group_summary, units, group_key, dataset_key=dataset_key)
    lodo.to_csv(tables / "lodo.csv", index=False)
    if lodo["held_out_dataset"].nunique() >= 2:
        plot_lodo_heatmap(lodo, group_key, config, output_stem=stem["13"])
    plot_depth(units.loc[units["eligible"]], group_key, gene, colors, config, output_stem=stem["14"])
    return stem


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Render the Part 2 TARGET_GENE figure catalog.")
    parser.add_argument("h5ad")
    parser.add_argument("--gene", required=True)
    parser.add_argument("--group", default="cell_type")
    parser.add_argument("--colors", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--config", default=str(SCRIPT_DIR / "plotting_config.yaml"))
    parser.add_argument("--min-cells", type=int, default=20)
    parser.add_argument("--donor-key", default="donor_id")
    parser.add_argument("--dataset-key", default="dataset")
    parser.add_argument("--stem-prefix", default="F02")
    parser.add_argument("--umap-title", default="Locked global cell types")
    args = parser.parse_args(argv)
    if sc is None:
        raise SystemExit("anndata required")
    adata = sc.read_h5ad(args.h5ad)
    colors = load_color_map(args.colors)
    config = load_plotting_config(args.config)
    render_catalog(
        adata,
        gene=args.gene,
        group_key=args.group,
        colors=colors,
        config=config,
        out_dir=Path(args.out),
        min_cells=args.min_cells,
        donor_key=args.donor_key,
        dataset_key=args.dataset_key,
        stem_prefix=args.stem_prefix,
        umap_title=args.umap_title,
    )
    print(f"wrote {args.stem_prefix} catalog under {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

