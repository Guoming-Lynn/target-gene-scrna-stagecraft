#!/usr/bin/env python3
"""Part 5 association figures. Reads tables, not cells.

Usage:
    python scripts/part5_figures.py coverage eligibility.csv --out 03_figures
    python scripts/part5_figures.py exposure metadata.csv --out 03_figures
    python scripts/part5_figures.py collinear metadata.csv --out 03_figures
    python scripts/part5_figures.py volcano genes.csv --out 03_figures --gene SYMBOL
    python scripts/part5_figures.py forest gene_effects.csv --out 03_figures
    python scripts/part5_figures.py holdout holdout.csv --out 03_figures
    python scripts/part5_figures.py loo loo_with_source.csv --out 03_figures
    python scripts/part5_figures.py pathways pathway_evidence.csv --out 03_figures
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any, Mapping

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from plotting_style import apply_publication_style, figure_size, load_plotting_config, save_figure  # noqa: E402

STATUS_COLORS = {
    "formal": "#009E73",
    "exploratory": "#E69F00",
    "NOT_ESTIMABLE": "#808080",
    "SUCCESS": "#0072B2",
    "RANGE_OK": "#009E73",
    "NO_EXPOSURE_RANGE": "#D55E00",
}
TECHNICAL_PREFIXES = ("MT-", "RPS", "RPL")
TECHNICAL_EXACT = {"MALAT1", "XIST"}


def _is_technical(name: str) -> bool:
    upper = str(name).upper()
    if upper in TECHNICAL_EXACT or upper.startswith(("HBA", "HBB", "HBD", "HBG", "HBM", "HBQ", "HBZ")):
        return True
    return upper.startswith(TECHNICAL_PREFIXES)


def _config(path: str | None) -> dict[str, Any]:
    config = load_plotting_config(path or SCRIPT_DIR / "plotting_config.yaml")
    apply_publication_style(config)
    return config


def plot_coverage(table: pd.DataFrame, config: Mapping[str, Any], stem: Path, input_files: list[Path]) -> None:
    frame = table.copy()
    if "exposure_sd" in frame.columns:
        frame = frame.sort_values("exposure_sd", ascending=False, na_position="last")
    labels = frame["arm"].astype(str).tolist()
    rdf = frame["rdf"].to_numpy(float) if "rdf" in frame.columns else np.zeros(len(frame))
    status = frame["status"].astype(str) if "status" in frame.columns else pd.Series(["formal"] * len(frame))
    colors = [STATUS_COLORS.get(str(s), "#808080") for s in status]
    fig, ax = plt.subplots(figsize=figure_size(config, "wide_panel"), layout="constrained")
    ax.barh(range(len(labels)), np.nan_to_num(rdf), color=colors, height=0.7)
    ax.set_yticks(range(len(labels)), labels, fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("Residual df of the declared design")
    ax.set_title("n is dataset × donor_id; rdf is of the declared design", loc="left", fontweight="bold")
    save_figure(fig, stem, config, parameters={"plot": "eligibility_rdf"}, input_files=input_files)


def plot_exposure(
    meta: pd.DataFrame,
    config: Mapping[str, Any],
    stem: Path,
    *,
    exposure: str,
    source_key: str,
    input_files: list[Path],
) -> None:
    work = meta.loc[meta["eligible"].astype(bool)].copy() if "eligible" in meta.columns else meta.copy()
    sources = list(dict.fromkeys(work[source_key].astype(str)))
    fig, ax = plt.subplots(figsize=figure_size(config, "wide_panel"), layout="constrained")
    data = [work.loc[work[source_key].astype(str).eq(s), exposure].to_numpy(float) for s in sources]
    ax.boxplot(data, showfliers=False, widths=0.55)
    ax.set_xticklabels(sources, rotation=30, ha="right")
    rng = np.random.default_rng(0)
    for i, values in enumerate(data, start=1):
        finite = values[np.isfinite(values)]
        if finite.size == 0:
            continue
        x = rng.uniform(-0.12, 0.12, size=finite.size) + i
        ax.scatter(x, finite, s=12, c="#202124", alpha=0.75, zorder=3, linewidths=0)
    ax.set_ylabel(exposure)
    ax.set_xlabel("Source")
    ax.set_title("Unit-level exposure; sd ≈ 0 cannot identify a slope", loc="left", fontweight="bold")
    save_figure(
        fig,
        stem,
        config,
        parameters={"plot": "exposure_by_source", "points": "eligible donor-units"},
        input_files=input_files,
    )


def plot_collinear(
    meta: pd.DataFrame,
    config: Mapping[str, Any],
    stem: Path,
    *,
    x: str,
    y: str,
    input_files: list[Path],
) -> None:
    work = meta.loc[meta["eligible"].astype(bool)].copy() if "eligible" in meta.columns else meta.copy()
    xv = work[x].to_numpy(float)
    yv = work[y].to_numpy(float)
    mask = np.isfinite(xv) & np.isfinite(yv)
    rho = float(pd.Series(xv[mask]).corr(pd.Series(yv[mask]), method="spearman")) if mask.sum() >= 3 else float("nan")
    fig, ax = plt.subplots(figsize=figure_size(config, "single_panel"), layout="constrained")
    ax.scatter(xv[mask], yv[mask], s=16, c="#0072B2", alpha=0.8, linewidths=0)
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    ax.set_title("these two scales may be collinear", loc="left", fontweight="bold")
    ax.text(0.04, 0.96, f"Spearman ρ = {rho:.2f}", transform=ax.transAxes, va="top", fontsize=7)
    save_figure(fig, stem, config, parameters={"plot": "scale_collinearity", "spearman": rho}, input_files=input_files)


def plot_volcano(
    genes: pd.DataFrame,
    config: Mapping[str, Any],
    stem: Path,
    *,
    effect: str,
    pcol: str,
    qcol: str,
    logfc_floor: float,
    q_cut: float,
    input_files: list[Path],
) -> None:
    frame = genes.copy()
    name_col = "gene" if "gene" in frame.columns else frame.columns[0]
    x = frame[effect].to_numpy(float)
    p = frame[pcol].to_numpy(float)
    y = -np.log10(np.clip(p, 1e-300, 1))
    q = frame[qcol].to_numpy(float) if qcol in frame.columns else np.full(len(frame), np.nan)
    robust = frame["robust_primary"].astype(bool).to_numpy() if "robust_primary" in frame.columns else (
        (np.abs(x) >= logfc_floor) & (q < q_cut)
    )
    technical = frame[name_col].astype(str).map(_is_technical).to_numpy()
    fig, ax = plt.subplots(figsize=figure_size(config, "single_panel"), layout="constrained")
    ax.scatter(x[~robust & ~technical], y[~robust & ~technical], s=6, c="#B0B0B0", linewidths=0, rasterized=True)
    ax.scatter(x[technical], y[technical], s=6, c="#808080", linewidths=0, rasterized=True, label="technical")
    ax.scatter(x[robust & ~technical], y[robust & ~technical], s=10, c="#D55E00", linewidths=0, label="robust_primary")
    ax.axvline(logfc_floor, color="#D9D9D9", lw=0.6)
    ax.axvline(-logfc_floor, color="#D9D9D9", lw=0.6)
    labelled = frame.loc[robust & ~technical].copy()
    labelled["_abs"] = np.abs(labelled[effect].to_numpy(float))
    labelled = labelled.sort_values("_abs", ascending=False).head(20)
    for _, row in labelled.iterrows():
        ax.text(row[effect], -np.log10(max(float(row[pcol]), 1e-300)), str(row[name_col]), fontsize=5.5, ha="left")
    ax.set_xlabel(effect)
    ax.set_ylabel(f"−log10 {pcol}")
    ax.set_title("Donor-unit gene model; TARGET_GENE is not an outcome point", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=6.5)
    save_figure(fig, stem, config, parameters={"plot": "donor_unit_volcano", "stars": "genes not cells"}, input_files=input_files)


def plot_forest(
    table: pd.DataFrame,
    config: Mapping[str, Any],
    stem: Path,
    *,
    gene_col: str,
    effect: str,
    lo: str,
    hi: str,
    title: str,
    input_files: list[Path],
) -> None:
    frame = table.copy()
    labels = frame[gene_col].astype(str).tolist()
    fig_h = max(figure_size(config, "wide_panel")[1], 0.28 * len(labels) + 1.4)
    fig, ax = plt.subplots(figsize=(figure_size(config, "wide_panel")[0], fig_h), layout="constrained")
    y = np.arange(len(labels))
    for i, row in frame.reset_index(drop=True).iterrows():
        status = str(row.get("status", "SUCCESS"))
        if status == "NOT_ESTIMABLE" or not np.isfinite(row.get(effect, np.nan)):
            ax.text(0, i, "NOT_ESTIMABLE", ha="center", va="center", fontsize=6.5, color="#808080")
            continue
        ax.plot([row[lo], row[hi]], [i, i], color="#0072B2", lw=1.2)
        ax.plot(row[effect], i, "o", color="#0072B2", ms=4)
    ax.axvline(0, color="#D9D9D9", lw=0.6)
    ax.set_yticks(y, labels, fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel(effect)
    ax.set_title(title, loc="left", fontweight="bold")
    save_figure(fig, stem, config, parameters={"plot": "forest"}, input_files=input_files)


def plot_holdout(
    table: pd.DataFrame,
    config: Mapping[str, Any],
    stem: Path,
    *,
    subset_col: str,
    effect: str,
    lo: str,
    hi: str,
    model_col: str,
    input_files: list[Path],
) -> None:
    frame = table.copy()
    models = list(dict.fromkeys(frame[model_col].astype(str))) if model_col in frame.columns else ["primary"]
    limma_models = [name for name in models if "limma" in str(name).lower()]
    if limma_models:
        models = limma_models[:2]
    subsets = list(dict.fromkeys(frame[subset_col].astype(str)))
    n_col = max(1, len(models))
    fig, axes = plt.subplots(
        1,
        n_col,
        figsize=figure_size(config, "double_panel" if n_col > 1 else "single_panel"),
        layout="constrained",
        sharey=True,
    )
    if n_col == 1:
        axes = [axes]
    for ax, model in zip(axes, models):
        sub = frame.loc[frame[model_col].astype(str).eq(model)] if model_col in frame.columns else frame
        sub = sub.set_index(subset_col).reindex(subsets)
        y = np.arange(len(subsets))
        for i, name in enumerate(subsets):
            row = sub.iloc[i] if i < len(sub) else None
            status = str(row.get("status", "")) if row is not None else "NOT_ESTIMABLE"
            value = row[effect] if row is not None else np.nan
            if status == "NOT_ESTIMABLE" or not np.isfinite(value):
                ax.text(0, i, "NOT_ESTIMABLE", ha="center", va="center", fontsize=6.5, color="#808080")
                continue
            ax.plot([row[lo], row[hi]], [i, i], color="#0072B2", lw=1.2)
            ax.plot(value, i, "o", color="#0072B2", ms=4.5)
        ax.axvline(0, color="#D9D9D9", lw=0.6)
        ax.set_yticks(y, subsets, fontsize=7)
        ax.set_title(str(model), loc="left", fontweight="bold")
        ax.set_xlabel(effect)
    axes[0].invert_yaxis()
    fig.suptitle("holdout of a source block, not of one donor", fontsize=9, fontweight="bold")
    save_figure(
        fig,
        stem,
        config,
        parameters={"plot": "source_block_holdout", "empty_slots": "NOT_ESTIMABLE kept"},
        input_files=input_files,
    )


def plot_loo(
    table: pd.DataFrame,
    config: Mapping[str, Any],
    stem: Path,
    *,
    effect: str,
    source_key: str,
    input_files: list[Path],
) -> None:
    frame = table.copy()
    sources = list(dict.fromkeys(frame[source_key].astype(str)))
    palette = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"]
    color = {name: palette[i % len(palette)] for i, name in enumerate(sources)}
    fig, ax = plt.subplots(figsize=figure_size(config, "wide_panel"), layout="constrained")
    x = np.arange(len(frame))
    y = frame[effect].to_numpy(float)
    c = frame[source_key].astype(str).map(color)
    ax.scatter(x, y, c=c, s=18, linewidths=0)
    ax.axhline(0, color="#D9D9D9", lw=0.6)
    ax.set_xlabel("Omitted donor (iteration)")
    ax.set_ylabel(effect)
    ax.set_title("donor LOO (leverage of one person), not source-block LODO", loc="left", fontweight="bold")
    handles = [plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=color[s], label=s) for s in sources]
    ax.legend(handles=handles, frameon=False, fontsize=6.5)
    save_figure(fig, stem, config, parameters={"plot": "donor_loo"}, input_files=input_files)


def plot_pathways(
    table: pd.DataFrame,
    config: Mapping[str, Any],
    stem: Path,
    input_files: list[Path],
) -> None:
    frame = table.copy()
    xcol = "camera_q" if "camera_q" in frame.columns else frame.columns[1]
    ycol = "fgsea_q" if "fgsea_q" in frame.columns else frame.columns[2]
    fig, ax = plt.subplots(figsize=figure_size(config, "single_panel"), layout="constrained")
    ax.scatter(
        -np.log10(np.clip(frame[xcol].to_numpy(float), 1e-300, 1)),
        -np.log10(np.clip(frame[ycol].to_numpy(float), 1e-300, 1)),
        s=10,
        c="#0072B2",
        linewidths=0,
    )
    ax.axvline(-np.log10(0.05), color="#D9D9D9", lw=0.6)
    ax.axhline(-np.log10(0.05), color="#D9D9D9", lw=0.6)
    ax.set_xlabel("−log10 CAMERA q")
    ax.set_ylabel("−log10 fgsea q")
    ax.set_title("robust_primary needs both; empty CAMERA is allowed", loc="left", fontweight="bold")
    save_figure(fig, stem, config, parameters={"plot": "pathway_concordance"}, input_files=input_files)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["coverage", "exposure", "collinear", "volcano", "forest", "holdout", "loo", "pathways", "within"])
    parser.add_argument("table", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--gene", default="TARGET_GENE")
    parser.add_argument("--config", default=str(SCRIPT_DIR / "plotting_config.yaml"))
    parser.add_argument("--source-key", default="source_block")
    parser.add_argument("--exposure", default="jeffreys_per_10pct")
    parser.add_argument("--effect", default="logFC")
    parser.add_argument("--lo", default="CI_L")
    parser.add_argument("--hi", default="CI_R")
    parser.add_argument("--pcol", default="P")
    parser.add_argument("--qcol", default="q_bh")
    parser.add_argument("--model-col", default="model")
    parser.add_argument("--subset-col", default="subset")
    parser.add_argument("--gene-col", default="gene")
    parser.add_argument("--x", default="jeffreys_per_10pct")
    parser.add_argument("--y", default="log2cpm")
    args = parser.parse_args(argv)

    config = _config(args.config)
    table = pd.read_csv(args.table)
    if args.command in {"holdout", "loo", "forest", "within"} and args.gene and args.gene != "TARGET_GENE":
        if "gene" in table.columns:
            table = table.loc[table["gene"].astype(str).eq(args.gene) | table["gene"].isna()].copy()
            if table.empty:
                raise SystemExit(f"no rows for focus gene {args.gene}")
    args.out.mkdir(parents=True, exist_ok=True)
    inputs = [args.table]
    if args.command == "coverage":
        plot_coverage(table, config, args.out / "F05_01_eligibility", inputs)
    elif args.command == "exposure":
        key = args.source_key if args.source_key in table.columns else "dataset"
        plot_exposure(table, config, args.out / "F05_02_exposure_by_source", exposure=args.exposure, source_key=key, input_files=inputs)
    elif args.command == "collinear":
        plot_collinear(table, config, args.out / "F05_03_scale_collinearity", x=args.x, y=args.y, input_files=inputs)
    elif args.command == "volcano":
        plot_volcano(
            table,
            config,
            args.out / f"F05_04_{args.gene}_volcano",
            effect=args.effect,
            pcol=args.pcol,
            qcol=args.qcol,
            logfc_floor=0.10,
            q_cut=0.05,
            input_files=inputs,
        )
    elif args.command == "forest":
        plot_forest(
            table,
            config,
            args.out / "F05_05_target_gene_forest",
            gene_col=args.gene_col,
            effect=args.effect,
            lo=args.lo,
            hi=args.hi,
            title="Pre-declared panel; blank is not zero",
            input_files=inputs,
        )
    elif args.command == "holdout":
        plot_holdout(
            table,
            config,
            args.out / "F05_08_source_holdout_forest",
            subset_col=args.subset_col,
            effect=args.effect,
            lo=args.lo,
            hi=args.hi,
            model_col=args.model_col,
            input_files=inputs,
        )
    elif args.command == "loo":
        key = args.source_key if args.source_key in table.columns else "dataset"
        plot_loo(table, config, args.out / "F05_07_donor_loo", effect=args.effect, source_key=key, input_files=inputs)
    elif args.command == "pathways":
        plot_pathways(table, config, args.out / "F05_06_pathway_camera_fgsea", inputs)
    elif args.command == "within":
        if "subset" in table.columns:
            table = table.loc[table["subset"].astype(str).str.startswith("WITHIN_")].copy()
        plot_forest(
            table,
            config,
            args.out / "F05_09_within_source",
            gene_col=args.subset_col if "subset" in table.columns else args.gene_col,
            effect=args.effect,
            lo=args.lo,
            hi=args.hi,
            title="within-source agreement is not cross-source replication",
            input_files=inputs,
        )
    print(f"wrote F05 panel under {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


