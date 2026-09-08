#!/usr/bin/env python3
"""Part 3 compartment-round figures.

Formal panels only. Does not name clusters, delete cells, or recompute HVG.
STOP POINT 3 lock is the only command that draws names.

Usage:
    python scripts/part3_figures.py stop1 clustered.h5ad --out figures/round_00_initial
    python scripts/part3_figures.py stop2 clustered.h5ad --leiden-key leiden_r0_5 \\
        --panels config/marker_panels.yaml --out figures/round_00_initial
    python scripts/part3_figures.py exclusion parent.h5ad --partition manifest.csv --out figures/round_00
    python scripts/part3_figures.py membership membership.csv --out figures/round_01
    python scripts/part3_figures.py lock annotated.h5ad --group subtype --colors colors.yaml --out figures/final
"""
from __future__ import annotations

import argparse
import re
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
)

try:
    import anndata as sc
except ImportError:  # pragma: no cover
    sc = None

try:
    import yaml
except ImportError:  # pragma: no cover
    yaml = None


def _coords(adata, key: str = "X_umap") -> np.ndarray:
    if key not in adata.obsm:
        raise KeyError(f"Missing obsm['{key}'].")
    return np.asarray(adata.obsm[key])[:, :2]


def _layer_vector(adata, gene: str, layer: str) -> np.ndarray:
    if gene not in adata.var_names:
        raise KeyError(gene)
    matrix = adata[:, gene].layers[layer]
    if sparse.issparse(matrix):
        return np.asarray(matrix.toarray()).ravel()
    return np.asarray(matrix).ravel()


def _sort_ids(ids: Sequence[str]) -> list[str]:
    return sorted(ids, key=lambda x: (0, int(x)) if str(x).isdigit() else (1, str(x)))


def parse_resolution(column: str) -> float:
    match = re.search(r"r(\d+)_(\d+)$", str(column))
    if not match:
        return float("nan")
    return float(f"{match.group(1)}.{match.group(2)}")


def leiden_columns(adata, keys: Sequence[str] | None = None) -> list[str]:
    if keys:
        missing = [k for k in keys if k not in adata.obs]
        if missing:
            raise SystemExit(f"Missing Leiden columns: {missing}")
        return list(keys)
    cols = [c for c in adata.obs.columns if re.search(r"leiden", str(c), re.I)]
    if not cols:
        raise SystemExit("No leiden* columns in obs.")
    return sorted(cols, key=lambda c: (np.nan_to_num(parse_resolution(c), nan=99), str(c)))


def cluster_color_map(cluster_ids: Sequence[str], config: Mapping[str, Any]) -> dict[str, str]:
    cycle = list(config["theme"]["cluster_cycle"])
    mapping: dict[str, str] = {}
    for cid in cluster_ids:
        text = str(cid)
        try:
            index = int(text)
        except ValueError:
            index = sum(ord(ch) for ch in text)
        mapping[text] = cycle[index % len(cycle)]
    return mapping


def unit_labels(obs: pd.DataFrame, dataset_key: str, donor_key: str) -> pd.Series:
    if "dataset_donor_id" in obs.columns:
        return obs["dataset_donor_id"].astype(str)
    return obs[dataset_key].astype(str) + "_" + obs[donor_key].astype(str)


def _source_dir(out_dir: Path) -> Path:
    path = out_dir / "source_data"
    path.mkdir(parents=True, exist_ok=True)
    return path


def write_source(frame: pd.DataFrame, out_dir: Path, stem: str) -> Path:
    path = _source_dir(out_dir) / f"{stem}.csv"
    frame.to_csv(path, index=False, encoding="utf-8-sig")
    return path


def draw_cluster_numbers(ax, xy: np.ndarray, labels: pd.Series) -> None:
    frame = pd.DataFrame({"x": xy[:, 0], "y": xy[:, 1], "c": labels.astype(str).to_numpy()})
    for cluster, subset in frame.groupby("c", observed=True):
        ax.text(
            float(subset["x"].median()),
            float(subset["y"].median()),
            str(cluster),
            ha="center",
            va="center",
            fontweight="bold",
            fontsize=7.0,
            color="#111111",
            bbox={
                "boxstyle": "round,pad=0.15",
                "facecolor": "white",
                "edgecolor": "#202124",
                "linewidth": 0.45,
                "alpha": 0.92,
            },
            zorder=4,
        )


def plot_umap_clusters(
    adata,
    leiden_key: str,
    config: Mapping[str, Any],
    output_stem: Path,
    *,
    title: str,
    names: Mapping[str, str] | None = None,
    named_colors: Mapping[str, str] | None = None,
    input_files: Sequence[Path] = (),
) -> None:
    xy = _coords(adata)
    labels = adata.obs[leiden_key].astype(str)
    counts = labels.value_counts()
    order = list(counts.index)
    colors = named_colors or cluster_color_map(order, config)
    unknown = config["theme"].get("unknown", "#808080")
    size, alpha = auto_point_style(adata.n_obs, config)
    fig, ax = plt.subplots(figsize=figure_size(config, "double_panel"), layout="constrained")
    for category in reversed(order):
        mask = labels.eq(category).to_numpy()
        ax.scatter(
            xy[mask, 0],
            xy[mask, 1],
            s=size,
            alpha=alpha,
            color=colors.get(category, unknown),
            linewidths=0,
            rasterized=True,
        )
    draw_cluster_numbers(ax, xy, labels)
    handles = []
    for category in order:
        n = int(counts.loc[category])
        if names:
            pretty = names.get(category, category)
            label = f"{category}: {pretty} (n={n:,})" if pretty != category else f"{pretty} (n={n:,})"
        else:
            label = f"{category} (n={n:,})"
        handles.append(
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="",
                markersize=4,
                color=colors.get(category, unknown),
                label=label,
            )
        )
    ax.legend(
        handles=handles,
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
    save_figure(
        fig,
        output_stem,
        config,
        parameters={
            "plot": "umap_clusters",
            "leiden_key": leiden_key,
            "n_cells": int(adata.n_obs),
            "legend_order": "cell-count descending; small groups drawn last",
            "names_on_data": False,
            "numbers_on_data": True,
            "cluster_0_not_identity_across_rounds": True,
        },
        input_files=input_files,
    )


def plot_umap_categorical(
    adata,
    column: str,
    config: Mapping[str, Any],
    output_stem: Path,
    *,
    title: str,
    max_legend: int = 12,
    input_files: Sequence[Path] = (),
) -> None:
    xy = _coords(adata)
    labels = adata.obs[column].astype(str)
    counts = labels.value_counts()
    order = list(counts.index)
    cycle = list(config["theme"]["cluster_cycle"])
    colors = {name: cycle[i % len(cycle)] for i, name in enumerate(order)}
    size, alpha = auto_point_style(adata.n_obs, config)
    fig, ax = plt.subplots(figsize=figure_size(config, "double_panel"), layout="constrained")
    for category in reversed(order):
        mask = labels.eq(category).to_numpy()
        ax.scatter(
            xy[mask, 0],
            xy[mask, 1],
            s=size,
            alpha=alpha,
            color=colors[category],
            linewidths=0,
            rasterized=True,
        )
    ax.set_title(title, loc="left", fontweight="bold")
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_xticks([])
    ax.set_yticks([])
    if len(order) <= max_legend:
        handles = [
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="",
                markersize=4,
                color=colors[name],
                label=f"{name} (n={int(counts.loc[name]):,})",
            )
            for name in order
        ]
        ax.legend(
            handles=handles,
            bbox_to_anchor=(1.02, 0.5),
            loc="center left",
            frameon=False,
            fontsize=6.5,
        )
    save_figure(
        fig,
        output_stem,
        config,
        parameters={
            "plot": "umap_categorical",
            "column": column,
            "n_categories": len(order),
            "legend": "drawn" if len(order) <= max_legend else f"omitted (n>{max_legend})",
        },
        input_files=input_files,
    )


def plot_umap_continuous(
    adata,
    values: np.ndarray,
    config: Mapping[str, Any],
    output_stem: Path,
    *,
    title: str,
    colorbar_label: str,
    input_files: Sequence[Path] = (),
) -> None:
    xy = _coords(adata)
    display = np.asarray(values, dtype=float)
    finite = np.isfinite(display)
    vmin = float(np.nanmin(display[finite])) if finite.any() else 0.0
    vmax = float(np.nanmax(display[finite])) if finite.any() else 1.0
    if vmin == vmax:
        vmax = vmin + 1.0
    size, alpha = auto_point_style(adata.n_obs, config)
    order = np.argsort(np.nan_to_num(display, nan=vmin), kind="stable")
    fig, ax = plt.subplots(figsize=(4.15, 3.65), layout="constrained")
    artist = ax.scatter(
        xy[order, 0],
        xy[order, 1],
        c=display[order],
        cmap=config["theme"]["colorblind_continuous"],
        vmin=vmin,
        vmax=vmax,
        s=size,
        alpha=alpha,
        linewidths=0,
        rasterized=True,
    )
    cbar = fig.colorbar(artist, ax=ax, orientation="horizontal", fraction=0.055, pad=0.08, aspect=35)
    cbar.set_label(colorbar_label)
    ax.set_title(title, loc="left", fontweight="bold")
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_xticks([])
    ax.set_yticks([])
    save_figure(
        fig,
        output_stem,
        config,
        parameters={"plot": "umap_continuous", "title": title, "vmin": vmin, "vmax": vmax},
        input_files=input_files,
    )


def plot_hvg(
    adata,
    config: Mapping[str, Any],
    output_stem: Path,
    *,
    gene: str = "",
    input_files: Sequence[Path] = (),
) -> None:
    var = adata.var
    mean_col = next((c for c in ("means", "mean", "mean_counts") if c in var.columns), None)
    disp_col = next(
        (c for c in ("variances_norm", "dispersions_norm", "dispersions") if c in var.columns),
        None,
    )
    if mean_col is None or disp_col is None or "highly_variable" not in var.columns:
        return
    frame = pd.DataFrame(
        {
            "gene": var.index.astype(str),
            "mean": var[mean_col].to_numpy(float),
            "dispersion": var[disp_col].to_numpy(float),
            "highly_variable": var["highly_variable"].to_numpy(bool),
        }
    )
    fig, ax = plt.subplots(figsize=figure_size(config, "single_panel"), layout="constrained")
    other = ~frame["highly_variable"]
    ax.scatter(
        frame.loc[other, "mean"],
        frame.loc[other, "dispersion"],
        s=6,
        color="#B0B0B0",
        alpha=0.35,
        linewidths=0,
        rasterized=True,
        label="not HVG",
    )
    ax.scatter(
        frame.loc[frame["highly_variable"], "mean"],
        frame.loc[frame["highly_variable"], "dispersion"],
        s=8,
        color="#0072B2",
        alpha=0.7,
        linewidths=0,
        rasterized=True,
        label="HVG",
    )
    in_hvg = False
    if gene and gene in set(frame["gene"]):
        row = frame.loc[frame["gene"].eq(gene)].iloc[0]
        in_hvg = bool(row["highly_variable"])
        ax.scatter(
            [row["mean"]],
            [row["dispersion"]],
            s=36,
            facecolors="none",
            edgecolors="#D73027",
            linewidths=1.1,
            zorder=3,
            label=gene,
        )
    ax.set_title(f"{gene or 'TARGET_GENE'} excluded from PCA input", loc="left", fontweight="bold")
    ax.set_xlabel("Mean expression")
    ax.set_ylabel("Normalized variance")
    ax.legend(frameon=False, fontsize=6.5)
    save_figure(
        fig,
        output_stem,
        config,
        parameters={
            "plot": "hvg_mean_dispersion",
            "target_gene": gene,
            "target_gene_in_highly_variable": in_hvg,
            "failed_assert_if_in_hvg": bool(gene) and in_hvg,
        },
        input_files=input_files,
    )


def plot_pca_variance(
    adata,
    config: Mapping[str, Any],
    output_stem: Path,
    *,
    n_used: int = 30,
    input_files: Sequence[Path] = (),
) -> None:
    if "pca" not in adata.uns or "variance_ratio" not in adata.uns["pca"]:
        return
    ratios = np.asarray(adata.uns["pca"]["variance_ratio"], dtype=float)
    frame = pd.DataFrame({"pc": np.arange(1, len(ratios) + 1), "variance_ratio": ratios})
    fig, ax = plt.subplots(figsize=figure_size(config, "single_panel"), layout="constrained")
    ax.bar(frame["pc"], frame["variance_ratio"], color="#0072B2", width=0.85)
    ax.axvline(n_used + 0.5, color="#D55E00", linewidth=0.9, linestyle="--", label=f"used {n_used} PCs")
    ax.set_title("PCA variance after TARGET_GENE exclusion", loc="left", fontweight="bold")
    ax.set_xlabel("PC (intrinsic order)")
    ax.set_ylabel("Variance ratio")
    ax.legend(frameon=False, fontsize=6.5)
    save_figure(
        fig,
        output_stem,
        config,
        parameters={"plot": "pca_variance", "n_pcs_used": n_used, "n_pcs_computed": int(len(ratios))},
        input_files=input_files,
    )


def plot_harmony_before_after(
    adata,
    config: Mapping[str, Any],
    output_stem: Path,
    *,
    dataset_key: str = "dataset",
    input_files: Sequence[Path] = (),
) -> None:
    after_key = "X_pca" if "X_pca" in adata.obsm else None
    before_key = "X_pca_unintegrated" if "X_pca_unintegrated" in adata.obsm else None
    if after_key is None:
        return
    labels = adata.obs[dataset_key].astype(str)
    order = list(labels.value_counts().index)
    cycle = list(config["theme"]["cluster_cycle"])
    colors = {name: cycle[i % len(cycle)] for i, name in enumerate(order)}
    size, alpha = auto_point_style(adata.n_obs, config)
    n_panels = 2 if before_key else 1
    fig, axes = plt.subplots(1, n_panels, figsize=(3.6 * n_panels + 1.4, 3.4), layout="constrained")
    if n_panels == 1:
        axes = [axes]
    panels = []
    if before_key:
        panels.append((axes[0], before_key, "Unintegrated PCs"))
        panels.append((axes[1], after_key, "Harmony-corrected PCs"))
    else:
        panels.append((axes[0], after_key, "Harmony-corrected PCs (unintegrated not stashed)"))
    for ax, key, title in panels:
        xy = np.asarray(adata.obsm[key])[:, :2]
        for category in reversed(order):
            mask = labels.eq(category).to_numpy()
            ax.scatter(
                xy[mask, 0],
                xy[mask, 1],
                s=size,
                alpha=alpha,
                color=colors[category],
                linewidths=0,
                rasterized=True,
            )
        ax.set_title(title, loc="left", fontweight="bold")
        ax.set_xlabel("PC 1")
        ax.set_ylabel("PC 2")
        ax.set_xticks([])
        ax.set_yticks([])
    handles = [
        Line2D([0], [0], marker="o", linestyle="", markersize=4, color=colors[name], label=name)
        for name in order
    ]
    fig.legend(handles=handles, loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False, fontsize=6.5)
    save_figure(
        fig,
        output_stem,
        config,
        parameters={
            "plot": "harmony_before_after",
            "dataset_key": dataset_key,
            "unintegrated_stashed": before_key is not None,
        },
        input_files=input_files,
    )


def coverage_table(
    adata,
    leiden_key: str,
    *,
    dataset_key: str = "dataset",
    donor_key: str = "donor_id",
) -> pd.DataFrame:
    units = unit_labels(adata.obs, dataset_key, donor_key)
    clusters = adata.obs[leiden_key].astype(str)
    frame = pd.crosstab(units, clusters)
    frame = frame.reindex(columns=_sort_ids(frame.columns.astype(str)))
    return frame


def plot_coverage_heatmap(
    counts: pd.DataFrame,
    config: Mapping[str, Any],
    output_stem: Path,
    *,
    title: str,
    log1p: bool = True,
    input_files: Sequence[Path] = (),
) -> None:
    values = np.log1p(counts.to_numpy(float)) if log1p else counts.to_numpy(float)
    fig_w = max(6.4, 0.32 * counts.shape[1] + 2.8)
    fig_h = max(4.2, 0.22 * counts.shape[0] + 1.5)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), layout="constrained")
    image = ax.imshow(values, aspect="auto", cmap=config["theme"]["expression_continuous"], interpolation="nearest")
    ax.set_xticks(range(counts.shape[1]), counts.columns.astype(str), fontsize=6.5)
    ax.set_yticks(range(counts.shape[0]), counts.index.astype(str), fontsize=5.8)
    ax.set_xlabel("Leiden cluster")
    ax.set_ylabel("dataset × donor")
    ax.set_title(title, loc="left", fontweight="bold")
    cbar = fig.colorbar(image, ax=ax, pad=0.01)
    cbar.set_label("log1p cell count" if log1p else "cell count")
    save_figure(
        fig,
        output_stem,
        config,
        parameters={
            "plot": "coverage_heatmap",
            "log1p_display": log1p,
            "n_units": int(counts.shape[0]),
            "n_clusters": int(counts.shape[1]),
        },
        input_files=input_files,
    )


def resolution_overview(
    adata,
    columns: Sequence[str],
    *,
    dataset_key: str = "dataset",
    donor_key: str = "donor_id",
) -> pd.DataFrame:
    units = unit_labels(adata.obs, dataset_key, donor_key)
    rows = []
    for column in columns:
        labels = adata.obs[column].astype(str)
        sizes = labels.value_counts()
        composition = pd.crosstab(labels, adata.obs[dataset_key].astype(str), normalize="index")
        rows.append(
            {
                "leiden_key": column,
                "resolution": parse_resolution(column),
                "n_clusters": int(labels.nunique()),
                "min_cluster_size": int(sizes.min()),
                "median_cluster_size": float(sizes.median()),
                "n_units": int(units.nunique()),
                "max_dataset_fraction": float(composition.max(axis=1).max()) if not composition.empty else np.nan,
            }
        )
    return pd.DataFrame(rows)


def plot_resolution_landscape(
    summary: pd.DataFrame,
    config: Mapping[str, Any],
    output_stem: Path,
    *,
    input_files: Sequence[Path] = (),
) -> None:
    frame = summary.sort_values("resolution")
    fig, axes = plt.subplots(1, 3, figsize=(8.8, 2.9), layout="constrained")
    specs = [
        ("n_clusters", "N clusters"),
        ("min_cluster_size", "Min cluster size"),
        ("max_dataset_fraction", "Max dataset fraction"),
    ]
    for ax, (column, ylabel) in zip(axes, specs):
        ax.plot(frame["resolution"], frame[column], marker="o", color="#0072B2", linewidth=1.1)
        ax.set_xlabel("Leiden resolution")
        ax.set_ylabel(ylabel)
        ax.set_title(ylabel, loc="left", fontweight="bold")
    save_figure(
        fig,
        output_stem,
        config,
        parameters={"plot": "resolution_landscape", "resolution_axis": "intrinsic numeric order"},
        input_files=input_files,
    )


def plot_graph_sensitivity(
    table: pd.DataFrame,
    config: Mapping[str, Any],
    output_stem: Path,
    *,
    input_files: Sequence[Path] = (),
) -> None:
    def _col(options: tuple[str, ...]) -> str:
        for name in options:
            if name in table.columns:
                return name
        raise SystemExit(f"graph-sensitivity CSV needs one of {options}")

    pcs_col = _col(("n_pcs", "pcs", "n_pc"))
    nn_col = _col(("n_neighbors", "neighbors", "n_neighbour"))
    val_col = _col(("n_clusters", "metric", "value"))
    pivot = table.pivot_table(index=nn_col, columns=pcs_col, values=val_col, aggfunc="mean")
    pivot = pivot.sort_index().reindex(sorted(pivot.columns), axis=1)
    fig, ax = plt.subplots(figsize=figure_size(config, "single_panel"), layout="constrained")
    image = ax.imshow(pivot.to_numpy(float), cmap="cividis", aspect="auto")
    ax.set_xticks(range(pivot.shape[1]), [str(c) for c in pivot.columns])
    ax.set_yticks(range(pivot.shape[0]), [str(i) for i in pivot.index])
    ax.set_xlabel("PCs")
    ax.set_ylabel("Neighbors")
    ax.set_title("Graph sensitivity (n clusters)", loc="left", fontweight="bold")
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            value = pivot.to_numpy(float)[i, j]
            if np.isfinite(value):
                ax.text(j, i, f"{value:.0f}", ha="center", va="center", fontsize=7, color="white")
    fig.colorbar(image, ax=ax, pad=0.01)
    save_figure(
        fig,
        output_stem,
        config,
        parameters={"plot": "graph_sensitivity", "value": val_col},
        input_files=input_files,
    )


def plot_cluster_size_bar(
    adata,
    leiden_key: str,
    config: Mapping[str, Any],
    output_stem: Path,
    *,
    donor_key: str = "donor_id",
    input_files: Sequence[Path] = (),
) -> None:
    labels = adata.obs[leiden_key].astype(str)
    counts = labels.value_counts()
    donors = adata.obs.groupby(labels, observed=True)[donor_key].nunique() if donor_key in adata.obs else None
    colors = cluster_color_map(counts.index.astype(str), config)
    fig, ax = plt.subplots(figsize=figure_size(config, "wide_panel"), layout="constrained")
    x = np.arange(len(counts))
    ax.bar(x, counts.to_numpy(), color=[colors[str(i)] for i in counts.index], width=0.72)
    ax.set_xticks(x, [str(i) for i in counts.index])
    ax.set_ylabel("N cells")
    ax.set_xlabel("Leiden cluster (cell-count descending)")
    ax.set_title("Cluster size", loc="left", fontweight="bold")
    if donors is not None:
        for idx, cluster in enumerate(counts.index.astype(str)):
            ax.text(idx, counts.iloc[idx], f" d={int(donors.loc[cluster])}", ha="center", va="bottom", fontsize=6)
    save_figure(
        fig,
        output_stem,
        config,
        parameters={"plot": "cluster_size_bar", "order": "cell-count descending"},
        input_files=input_files,
    )


def _metric_values(adata, column: str) -> np.ndarray:
    if column == "log10_umi":
        umi = adata.obs["total_counts"].to_numpy(float) if "total_counts" in adata.obs else np.asarray(
            adata.layers["counts"].sum(axis=1)
        ).ravel()
        return np.log10(umi + 1)
    return adata.obs[column].to_numpy(float)


def plot_qc_violins(
    adata,
    leiden_key: str,
    config: Mapping[str, Any],
    output_stem: Path,
    *,
    input_files: Sequence[Path] = (),
) -> None:
    labels = adata.obs[leiden_key].astype(str)
    cluster_ids = _sort_ids(labels.unique())
    colors = cluster_color_map(cluster_ids, config)
    panels = []
    if "pct_counts_mt" in adata.obs:
        panels.append(("pct_counts_mt", "% mitochondrial counts", False))
    if "total_counts" in adata.obs or "counts" in adata.layers:
        panels.append(("log10_umi", "log10(UMI + 1)", False))
    if "n_genes_by_counts" in adata.obs:
        panels.append(("n_genes_by_counts", "Detected genes", False))
    if not panels:
        return
    fig, axes = plt.subplots(1, len(panels), figsize=(4.1 * len(panels), 3.6), layout="constrained")
    if len(panels) == 1:
        axes = [axes]
    panel_orders: dict[str, list[str]] = {}
    median_hex = config["theme"].get("median", "#B2182B")
    for ax, (column, ylabel, _) in zip(axes, panels):
        values = {cid: _metric_values(adata, column)[labels.eq(cid).to_numpy()] for cid in cluster_ids}
        order = sorted(
            cluster_ids,
            key=lambda cid: (-float(np.nanmedian(values[cid])) if np.isfinite(np.nanmedian(values[cid])) else 0, cid),
        )
        panel_orders[column] = order
        bodies = []
        for cid in order:
            arr = np.asarray(values[cid], dtype=float)
            arr = arr[np.isfinite(arr)]
            bodies.append(arr if len(arr) else np.array([np.nan]))
        positions = np.arange(1, len(order) + 1)
        violin_pos = [pos for pos, arr in zip(positions, bodies) if len(arr) >= 2]
        violin_data = [arr for arr in bodies if len(arr) >= 2]
        if violin_data:
            violin = ax.violinplot(violin_data, positions=violin_pos, showextrema=False, showmedians=False)
            body_ids = [cid for cid, arr in zip(order, bodies) if len(arr) >= 2]
            for body, cid in zip(violin["bodies"], body_ids):
                body.set_facecolor(colors[cid])
                body.set_edgecolor(colors[cid])
                body.set_alpha(0.45)
        ax.boxplot(
            bodies,
            positions=positions,
            widths=0.18,
            patch_artist=True,
            showfliers=False,
            whis=(5, 95),
            boxprops={"facecolor": "white", "edgecolor": "#222222", "linewidth": 0.65, "alpha": 0.85},
            medianprops={"color": median_hex, "linewidth": 1.25},
            whiskerprops={"color": "#222222", "linewidth": 0.65},
            capprops={"color": "#222222", "linewidth": 0.65},
        )
        ax.set_xticks(positions, order)
        ax.set_xlabel("Leiden cluster (median descending)")
        ax.set_ylabel(ylabel)
        ax.set_title(ylabel, loc="left", fontweight="bold")
    save_figure(
        fig,
        output_stem,
        config,
        parameters={
            "plot": "qc_violin_box",
            "ordering": "each panel independently ordered by its displayed metric median, descending",
            "panel_cluster_orders": panel_orders,
            "strip_points": False,
            "color_follows_cluster_id": True,
        },
        input_files=input_files,
    )


def load_panels(path: Path) -> dict[str, Any]:
    if yaml is None:
        raise SystemExit("PyYAML required: pip install pyyaml")
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def marker_matrix(
    adata,
    leiden_key: str,
    genes: Sequence[str],
    *,
    exclude: set[str] | None = None,
) -> tuple[pd.DataFrame, list[str]]:
    if "counts" not in adata.layers or "normalized" not in adata.layers:
        raise ValueError("Marker evidence requires counts and normalized layers.")
    wanted = [g for g in dict.fromkeys(genes) if not exclude or g.upper() not in exclude]
    missing = [g for g in wanted if g not in adata.var_names]
    present = [g for g in wanted if g in adata.var_names]
    labels = adata.obs[leiden_key].astype(str)
    rows = []
    for gene in present:
        counts = _layer_vector(adata, gene, "counts") if "counts" in adata.layers else np.zeros(adata.n_obs)
        norm = _layer_vector(adata, gene, "normalized") if "normalized" in adata.layers else counts
        for cluster, mask in labels.groupby(labels, observed=True):
            index = labels.eq(cluster).to_numpy()
            rows.append(
                {
                    "gene": gene,
                    "cluster": str(cluster),
                    "detection": float((counts[index] > 0).mean()),
                    "mean_normalized": float(norm[index].mean()),
                }
            )
    return pd.DataFrame(rows), missing


def plot_marker_matrix(
    evidence: pd.DataFrame,
    cluster_order: Sequence[str],
    gene_order: Sequence[str],
    config: Mapping[str, Any],
    output_stem: Path,
    *,
    title: str,
    input_files: Sequence[Path] = (),
) -> None:
    if evidence.empty or not gene_order:
        return
    subset = evidence.loc[evidence.gene.isin(gene_order) & evidence.cluster.isin(cluster_order)]
    x_map = {gene: i for i, gene in enumerate(gene_order)}
    y_map = {cluster: i for i, cluster in enumerate(cluster_order)}
    fig, ax = plt.subplots(
        figsize=(max(7.2, 0.28 * len(gene_order) + 2.6), max(3.4, 0.38 * len(cluster_order) + 1.6)),
        layout="constrained",
    )
    plotted = ax.scatter(
        subset["gene"].map(x_map),
        subset["cluster"].map(y_map),
        s=18 + 260 * subset["detection"],
        c=subset["mean_normalized"],
        cmap="cividis",
        vmin=0,
        linewidths=0.25,
        edgecolors="white",
    )
    ax.set_xticks(range(len(gene_order)), gene_order, rotation=60, ha="right", fontsize=6.2)
    ax.set_yticks(range(len(cluster_order)), [f"cluster {c}" for c in cluster_order])
    ax.invert_yaxis()
    ax.set_xlim(-0.7, len(gene_order) - 0.3)
    ax.set_title(title, loc="left", fontweight="bold")
    ax.set_xlabel("Genes; dot area = raw-count detection")
    cbar = fig.colorbar(plotted, ax=ax, pad=0.01)
    cbar.set_label("Mean log-normalized expression")
    save_figure(
        fig,
        output_stem,
        config,
        parameters={
            "plot": "marker_dot_matrix",
            "row_order": "cluster cell-count descending",
            "column_order": "max mean log-normalized descending",
            "p_values_shown": False,
        },
        input_files=input_files,
    )


def plot_exclusion_umap(
    adata,
    status: pd.Series,
    config: Mapping[str, Any],
    output_stem: Path,
    *,
    input_files: Sequence[Path] = (),
) -> None:
    xy = _coords(adata)
    aligned = status.reindex(adata.obs_names).fillna("unknown").astype(str).str.lower()
    size, alpha = auto_point_style(adata.n_obs, config)
    fig, ax = plt.subplots(figsize=figure_size(config, "double_panel"), layout="constrained")
    palette = {
        "retained": config["theme"].get("retained", "#0072B2"),
        "removed": config["theme"].get("removed", "#D55E00"),
        "unknown": config["theme"].get("unknown", "#808080"),
    }
    for key in ("retained", "unknown", "removed"):
        mask = aligned.eq(key).to_numpy()
        if not mask.any():
            continue
        ax.scatter(
            xy[mask, 0],
            xy[mask, 1],
            s=size,
            alpha=alpha,
            color=palette[key],
            linewidths=0,
            rasterized=True,
            label=f"{key} (n={int(mask.sum()):,})",
        )
    ax.legend(frameon=False, loc="best", fontsize=7)
    ax.set_title("Parent-round UMAP: retained vs removed", loc="left", fontweight="bold")
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_xticks([])
    ax.set_yticks([])
    save_figure(
        fig,
        output_stem,
        config,
        parameters={
            "plot": "exclusion_umap",
            "coordinates": "parent round",
            "child_manifold_recomputed": True,
        },
        input_files=input_files,
    )


def plot_membership_heatmap(
    counts: pd.DataFrame,
    config: Mapping[str, Any],
    output_stem: Path,
    *,
    input_files: Sequence[Path] = (),
) -> None:
    counts = counts.copy()
    for key in ("parent_cluster", "child_cluster"):
        if counts[key].isna().any():
            raise ValueError("Membership cluster IDs must not be missing.")
        counts[key] = counts[key].astype(str)
    parent_ids = _sort_ids(counts["parent_cluster"].astype(str).unique())
    child_ids = _sort_ids(counts["child_cluster"].astype(str).unique())
    matrix = (
        counts.pivot_table(index="parent_cluster", columns="child_cluster", values="n_cells", aggfunc="sum")
        .reindex(index=parent_ids, columns=child_ids)
        .fillna(0)
    )
    fig_w = max(5.6, 0.42 * matrix.shape[1] + 2.4)
    fig_h = max(4.2, 0.38 * matrix.shape[0] + 1.8)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), layout="constrained")
    image = ax.imshow(matrix.to_numpy(float), cmap="cividis", aspect="auto")
    ax.set_xticks(range(matrix.shape[1]), matrix.columns.astype(str))
    ax.set_yticks(range(matrix.shape[0]), matrix.index.astype(str))
    ax.set_xlabel("Child Leiden ID (new)")
    ax.set_ylabel("Parent Leiden ID (old)")
    ax.set_title("Old IDs are not new IDs", loc="left", fontweight="bold")
    if matrix.size <= 400:
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                value = matrix.to_numpy(float)[i, j]
                if value > 0:
                    ax.text(j, i, f"{int(value)}", ha="center", va="center", fontsize=6.2, color="white")
    fig.colorbar(image, ax=ax, pad=0.01, label="shared barcodes")
    save_figure(
        fig,
        output_stem,
        config,
        parameters={"plot": "membership_heatmap", "ids_comparable": False},
        input_files=input_files,
    )


def plot_partition_counts(
    manifest: pd.DataFrame,
    config: Mapping[str, Any],
    output_stem: Path,
    *,
    input_files: Sequence[Path] = (),
) -> None:
    status = manifest["status"].astype(str).str.lower()
    counts = status.value_counts().reindex(["retained", "removed"]).fillna(0)
    fig, ax = plt.subplots(figsize=figure_size(config, "single_panel"), layout="constrained")
    ax.bar(
        ["retained", "removed"],
        counts.to_numpy(),
        color=[config["theme"].get("retained", "#0072B2"), config["theme"].get("removed", "#D55E00")],
        width=0.55,
    )
    ax.set_ylabel("N cells")
    ax.set_title("Partition counts", loc="left", fontweight="bold")
    save_figure(
        fig,
        output_stem,
        config,
        parameters={"plot": "partition_counts", "n_retained": int(counts.get("retained", 0)), "n_removed": int(counts.get("removed", 0))},
        input_files=input_files,
    )


def _mkdir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def render_stop1(
    adata,
    config: Mapping[str, Any],
    out_dir: Path,
    *,
    gene: str = "",
    dataset_key: str = "dataset",
    donor_key: str = "donor_id",
    library_key: str = "library_id",
    coverage_key: str | None = None,
    graph_sensitivity: Path | None = None,
    n_pcs_used: int = 30,
    h5ad_path: Path | None = None,
) -> None:
    apply_publication_style(config)
    inputs = [h5ad_path] if h5ad_path else []
    cols = leiden_columns(adata)
    coverage_key = coverage_key or cols[min(len(cols) // 2, len(cols) - 1)]
    folders = {
        "hvg": _mkdir(out_dir / "01_HVG"),
        "pca": _mkdir(out_dir / "02_PCA"),
        "harmony": _mkdir(out_dir / "03_Harmony"),
        "umap": _mkdir(out_dir / "04_UMAP"),
        "meta": _mkdir(out_dir / "05_metadata"),
        "leiden": _mkdir(out_dir / "06_Leiden"),
    }
    plot_hvg(adata, config, folders["hvg"] / "F03_S1_01_hvg", gene=gene, input_files=inputs)
    plot_pca_variance(adata, config, folders["pca"] / "F03_S1_02_pca_variance", n_used=n_pcs_used, input_files=inputs)
    plot_harmony_before_after(
        adata, config, folders["harmony"] / "F03_S1_03_harmony_dataset", dataset_key=dataset_key, input_files=inputs
    )
    plot_umap_categorical(
        adata, dataset_key, config, folders["umap"] / "F03_S1_04_umap_dataset", title="UMAP by dataset", input_files=inputs
    )
    if donor_key in adata.obs:
        plot_umap_categorical(
            adata,
            donor_key,
            config,
            folders["meta"] / "F03_S1_05_umap_donor",
            title="UMAP by donor_id",
            input_files=inputs,
        )
    lib = library_key if library_key in adata.obs else ("sample_id" if "sample_id" in adata.obs else "")
    if lib:
        plot_umap_categorical(
            adata, lib, config, folders["meta"] / f"F03_S1_05_umap_{lib}", title=f"UMAP by {lib}", input_files=inputs
        )
    if "total_counts" in adata.obs:
        plot_umap_continuous(
            adata,
            np.log10(adata.obs["total_counts"].to_numpy(float) + 1),
            config,
            folders["meta"] / "F03_S1_05_umap_log_umi",
            title="UMAP by log10 UMI",
            colorbar_label="log10(UMI + 1)",
            input_files=inputs,
        )
    if "n_genes_by_counts" in adata.obs:
        plot_umap_continuous(
            adata,
            adata.obs["n_genes_by_counts"].to_numpy(float),
            config,
            folders["meta"] / "F03_S1_05_umap_n_genes",
            title="UMAP by n genes",
            colorbar_label="n genes",
            input_files=inputs,
        )
    if "scrublet_doublet_score" in adata.obs:
        plot_umap_continuous(
            adata,
            adata.obs["scrublet_doublet_score"].to_numpy(float),
            config,
            folders["meta"] / "F03_S1_05_umap_scrublet",
            title="UMAP by Scrublet score",
            colorbar_label="Scrublet score",
            input_files=inputs,
        )
    if gene and gene in adata.var_names and "normalized" in adata.layers:
        plot_umap_continuous(
            adata,
            _layer_vector(adata, gene, "normalized"),
            config,
            folders["meta"] / f"F03_S1_05_umap_{gene}_diagnostic",
            title=f"Diagnostic {gene} expression (not a cluster definition)",
            colorbar_label=f"Log-normalized {gene}",
            input_files=inputs,
        )
    for column in cols:
        res = parse_resolution(column)
        slug = f"r{str(res).replace('.', '_')}" if np.isfinite(res) else column
        plot_umap_clusters(
            adata,
            column,
            config,
            folders["leiden"] / f"F03_S1_06_leiden_{slug}",
            title=f"{column} (numbers only)",
            input_files=inputs,
        )
    summary = resolution_overview(adata, cols, dataset_key=dataset_key, donor_key=donor_key)
    write_source(summary, out_dir, "F03_S1_07_resolution_landscape")
    plot_resolution_landscape(
        summary, config, out_dir / "06_Leiden" / "F03_S1_07_resolution_landscape", input_files=inputs
    )
    if graph_sensitivity and Path(graph_sensitivity).is_file():
        table = pd.read_csv(graph_sensitivity)
        write_source(table, out_dir, "F03_S1_08_graph_sensitivity")
        plot_graph_sensitivity(
            table, config, out_dir / "03_Harmony" / "F03_S1_08_graph_sensitivity", input_files=list(inputs) + [Path(graph_sensitivity)]
        )
    coverage = coverage_table(adata, coverage_key, dataset_key=dataset_key, donor_key=donor_key)
    write_source(coverage.reset_index(names="dataset_donor_id"), out_dir, "F03_S1_09_coverage")
    plot_coverage_heatmap(
        coverage,
        config,
        out_dir / "06_Leiden" / "F03_S1_09_coverage",
        title=f"dataset × donor coverage ({coverage_key})",
        input_files=inputs,
    )


def render_stop2(
    adata,
    config: Mapping[str, Any],
    out_dir: Path,
    *,
    leiden_key: str,
    panels_path: Path | None = None,
    gene: str = "",
    dataset_key: str = "dataset",
    donor_key: str = "donor_id",
    h5ad_path: Path | None = None,
) -> None:
    apply_publication_style(config)
    inputs = [h5ad_path] if h5ad_path else []
    exclude = {gene.upper()} if gene else set()
    folders = {
        "leiden": _mkdir(out_dir / "06_Leiden"),
        "markers": _mkdir(out_dir / "07_markers"),
        "contam": _mkdir(out_dir / "08_contamination"),
        "qc": _mkdir(out_dir / "09_QC"),
    }
    plot_umap_clusters(
        adata,
        leiden_key,
        config,
        folders["leiden"] / "F03_S2_01_selected_umap",
        title=f"Selected {leiden_key} (numbers only)",
        input_files=inputs,
    )
    plot_cluster_size_bar(
        adata, leiden_key, config, folders["qc"] / "F03_S2_02_cluster_size", donor_key=donor_key, input_files=inputs
    )
    plot_qc_violins(adata, leiden_key, config, folders["qc"] / "F03_S2_03_qc_violins", input_files=inputs)
    coverage = coverage_table(adata, leiden_key, dataset_key=dataset_key, donor_key=donor_key)
    write_source(coverage.reset_index(names="dataset_donor_id"), out_dir, "F03_S2_04_coverage")
    plot_coverage_heatmap(
        coverage,
        config,
        folders["qc"] / "F03_S2_04_coverage",
        title=f"dataset × donor coverage ({leiden_key})",
        input_files=inputs,
    )
    cluster_order = list(adata.obs[leiden_key].astype(str).value_counts().index)
    if panels_path:
        panels = load_panels(panels_path)
        dropped: list[str] = []
        for block_name, folder, stem in (
            ("lineage", folders["markers"], "F03_S2_05_lineage_markers"),
            ("contamination", folders["contam"], "F03_S2_06_contamination_markers"),
        ):
            block = panels.get(block_name) or {}
            genes = []
            for group in (block.get("genes") or {}).values():
                genes.extend(list(group))
            if gene:
                kept = []
                for item in genes:
                    if str(item).upper() == gene.upper():
                        dropped.append(str(item))
                    else:
                        kept.append(item)
                genes = kept
            evidence, missing = marker_matrix(adata, leiden_key, genes, exclude=exclude)
            if not evidence.empty:
                gene_order = (
                    evidence.groupby("gene")["mean_normalized"].max().sort_values(ascending=False).index.tolist()
                )
                write_source(evidence, out_dir, stem)
                plot_marker_matrix(
                    evidence,
                    cluster_order,
                    gene_order,
                    config,
                    folder / stem,
                    title=block.get("title") or block_name,
                    input_files=list(inputs) + [panels_path],
                )
            if missing:
                (folder / f"{stem}_missing_genes.txt").write_text("\n".join(missing) + "\n", encoding="utf-8")
        if dropped:
            (_source_dir(out_dir) / "dropped_target_gene_from_panels.txt").write_text(
                "\n".join(dropped) + "\n", encoding="utf-8"
            )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    stop1 = sub.add_parser("stop1")
    stop1.add_argument("h5ad")
    stop1.add_argument("--out", required=True)
    stop1.add_argument("--gene", default="")
    stop1.add_argument("--dataset-key", default="dataset")
    stop1.add_argument("--donor-key", default="donor_id")
    stop1.add_argument("--library-key", default="library_id")
    stop1.add_argument("--coverage-key", default="")
    stop1.add_argument("--graph-sensitivity", default="")
    stop1.add_argument("--n-pcs-used", type=int, default=30)
    stop1.add_argument("--config", default=str(SCRIPT_DIR / "plotting_config.yaml"))

    stop2 = sub.add_parser("stop2")
    stop2.add_argument("h5ad")
    stop2.add_argument("--leiden-key", required=True)
    stop2.add_argument("--out", required=True)
    stop2.add_argument("--panels", default="")
    stop2.add_argument("--gene", default="")
    stop2.add_argument("--dataset-key", default="dataset")
    stop2.add_argument("--donor-key", default="donor_id")
    stop2.add_argument("--config", default=str(SCRIPT_DIR / "plotting_config.yaml"))

    exclusion = sub.add_parser("exclusion")
    exclusion.add_argument("h5ad")
    exclusion.add_argument("--partition", required=True)
    exclusion.add_argument("--out", required=True)
    exclusion.add_argument("--status-column", default="status")
    exclusion.add_argument("--barcode-column", default="parent_cell_id")
    exclusion.add_argument("--config", default=str(SCRIPT_DIR / "plotting_config.yaml"))

    membership = sub.add_parser("membership")
    membership.add_argument("csv")
    membership.add_argument("--out", required=True)
    membership.add_argument("--config", default=str(SCRIPT_DIR / "plotting_config.yaml"))

    lock = sub.add_parser("lock")
    lock.add_argument("h5ad")
    lock.add_argument("--group", required=True)
    lock.add_argument("--colors", required=True)
    lock.add_argument("--out", required=True)
    lock.add_argument("--id-key", default="", help="Numeric Leiden key for numbers-on-data; names come from --group")
    lock.add_argument("--config", default=str(SCRIPT_DIR / "plotting_config.yaml"))

    args = parser.parse_args(argv)
    if args.command != "membership" and sc is None:
        raise SystemExit("anndata required")
    config = load_plotting_config(args.config)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    if args.command == "stop1":
        adata = sc.read_h5ad(args.h5ad)
        render_stop1(
            adata,
            config,
            out,
            gene=args.gene,
            dataset_key=args.dataset_key,
            donor_key=args.donor_key,
            library_key=args.library_key,
            coverage_key=args.coverage_key or None,
            graph_sensitivity=Path(args.graph_sensitivity) if args.graph_sensitivity else None,
            n_pcs_used=args.n_pcs_used,
            h5ad_path=Path(args.h5ad),
        )
    elif args.command == "stop2":
        adata = sc.read_h5ad(args.h5ad)
        render_stop2(
            adata,
            config,
            out,
            leiden_key=args.leiden_key,
            panels_path=Path(args.panels) if args.panels else None,
            gene=args.gene,
            dataset_key=args.dataset_key,
            donor_key=args.donor_key,
            h5ad_path=Path(args.h5ad),
        )
    elif args.command == "exclusion":
        adata = sc.read_h5ad(args.h5ad)
        apply_publication_style(config)
        manifest = pd.read_csv(args.partition)
        status = pd.Series(
            manifest[args.status_column].astype(str).to_numpy(),
            index=manifest[args.barcode_column].astype(str),
        )
        write_source(manifest, out, "F03_R_01_exclusion")
        plot_exclusion_umap(
            adata,
            status,
            config,
            _mkdir(out / "10_removal") / "F03_R_01_exclusion_umap",
            input_files=[Path(args.h5ad), Path(args.partition)],
        )
        plot_partition_counts(
            manifest,
            config,
            out / "10_removal" / "F03_R_03_partition_counts",
            input_files=[Path(args.partition)],
        )
    elif args.command == "membership":
        apply_publication_style(config)
        table = pd.read_csv(args.csv)
        write_source(table, out, "F03_R_02_membership")
        plot_membership_heatmap(
            table,
            config,
            _mkdir(out / "10_removal") / "F03_R_02_membership",
            input_files=[Path(args.csv)],
        )
    else:
        adata = sc.read_h5ad(args.h5ad)
        apply_publication_style(config)
        colors = load_color_map(args.colors)
        id_key = args.id_key or args.group
        names = None
        named_colors = colors
        if args.id_key and args.id_key != args.group:
            names = {
                str(i): str(n)
                for i, n in zip(adata.obs[args.id_key].astype(str), adata.obs[args.group].astype(str))
            }
            named_colors = {
                str(i): colors.get(str(n), config["theme"].get("unknown", "#808080"))
                for i, n in names.items()
            }
        plot_umap_clusters(
            adata,
            id_key,
            config,
            _mkdir(out / "08_Final_annotation") / "F03_L_01_named_umap",
            title="Locked subtypes",
            names=names,
            named_colors=named_colors if args.id_key else colors,
            input_files=[Path(args.h5ad), Path(args.colors)],
        )
    print(f"wrote Part 3 figures under {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

