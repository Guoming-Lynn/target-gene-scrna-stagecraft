#!/usr/bin/env python3
"""Part 6 virtual-KO figures. Donor medians of Δaxis, not cells.

Usage:
    python scripts/part6_figures.py observability token_audit.csv --out 03_figures
    python scripts/part6_figures.py eligibility donor_eligibility.csv --out 03_figures
    python scripts/part6_figures.py donors donor_effects.csv --out 03_figures
    python scripts/part6_figures.py forest donor_effects.csv --out 03_figures
    python scripts/part6_figures.py ranks control_ranks.csv --out 03_figures
    python scripts/part6_figures.py support support_summary.csv --out 03_figures
    python scripts/part6_figures.py symmetry donor_effects.csv --out 03_figures
    python scripts/part6_figures.py sham sham.json --out 03_figures
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Mapping

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from plotting_style import apply_publication_style, figure_size, load_plotting_config, save_figure as _save_figure  # noqa: E402
from part6_sign_tests import median_order_statistic_interval  # noqa: E402


def save_figure(fig, stem, config, *, parameters=None, input_files=None):
    parameters = dict(parameters or {})
    parameters.update(evidence_class='exploratory_embedding_only', causal_inference=False,
                      expression_prediction=False)
    fig.text(0.99, 0.005, 'Exploratory embedding shift; not expression or causal evidence',
             ha='right', va='bottom', fontsize=5)
    return _save_figure(fig, stem, config, parameters=parameters, input_files=input_files)

KO_COLOR = "#0072B2"
OE_COLOR = "#D55E00"
PASS_COLOR = "#009E73"
FAIL_COLOR = "#808080"
ZERO_LINE = "#202124"


def _config(path: str | None) -> dict[str, Any]:
    config = load_plotting_config(path or SCRIPT_DIR / "plotting_config.yaml")
    apply_publication_style(config)
    return config


def plot_observability(table: pd.DataFrame, config: Mapping[str, Any], stem: Path, input_files: list[Path]) -> None:
    col = "token_class" if "token_class" in table.columns else "class"
    counts = table[col].astype(str).value_counts()
    labels = counts.index.astype(str).tolist()
    fig, ax = plt.subplots(figsize=figure_size(config, "single_panel"), layout="constrained")
    colors = [FAIL_COLOR if "ILLEGAL" in lab else KO_COLOR for lab in labels]
    ax.barh(range(len(labels)), counts.to_numpy(), color=colors, height=0.7)
    ax.set_yticks(range(len(labels)), labels, fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("Cells")
    ax.set_title("KO uses final-token-present cells only", loc="left", fontweight="bold")
    save_figure(fig, stem, config, parameters={"plot": "token_observability"}, input_files=input_files)


def plot_eligibility(table: pd.DataFrame, config: Mapping[str, Any], stem: Path, input_files: list[Path]) -> None:
    frame = table.copy()
    if "perturbation" in frame.columns:
        frame = frame.loc[frame["perturbation"].astype(str).str.upper().eq("KO")].copy()
    if frame.empty:
        frame = table.copy()
    donors = frame["dataset_donor_id"].astype(str).tolist()
    n = frame["n_success"].to_numpy(float) if "n_success" in frame.columns else np.zeros(len(frame))
    ok = frame["eligible"].astype(bool) if "eligible" in frame.columns else np.ones(len(frame), dtype=bool)
    colors = [PASS_COLOR if flag else FAIL_COLOR for flag in ok]
    fig, ax = plt.subplots(figsize=figure_size(config, "wide_panel"), layout="constrained")
    ax.barh(range(len(donors)), np.nan_to_num(n), color=colors, height=0.7)
    ax.set_yticks(range(len(donors)), donors, fontsize=6)
    ax.invert_yaxis()
    ax.set_xlabel("Successful KO cells / donor-unit")
    ax.set_title("eligibility before summaries", loc="left", fontweight="bold")
    save_figure(fig, stem, config, parameters={"plot": "donor_eligibility"}, input_files=input_files)


def plot_donors(
    table: pd.DataFrame,
    config: Mapping[str, Any],
    stem: Path,
    *,
    unpaired: bool,
    input_files: list[Path],
) -> None:
    frame = table.loc[table["eligible"].astype(bool)].copy() if "eligible" in table.columns else table.copy()
    endpoints = list(dict.fromkeys(frame["endpoint"].astype(str)))
    perts = [p for p in ("KO", "OE") if p in set(frame["perturbation"].astype(str).str.upper())]
    if not perts:
        perts = list(dict.fromkeys(frame["perturbation"].astype(str)))
    n_rows, n_cols = max(len(endpoints), 1), max(len(perts), 1)
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=figure_size(config, "multi_panel" if n_rows * n_cols > 2 else "double_panel"),
        layout="constrained",
        squeeze=False,
    )
    ko_order = None
    for i, endpoint in enumerate(endpoints):
        for j, pert in enumerate(perts):
            ax = axes[i][j]
            block = frame.loc[
                frame["endpoint"].astype(str).eq(endpoint)
                & frame["perturbation"].astype(str).str.upper().eq(str(pert).upper())
            ]
            if ko_order is None and str(pert).upper() == "KO":
                ko_order = (
                    block.sort_values("median_delta_axis")["dataset_donor_id"].astype(str).tolist()
                )
            order = ko_order or block["dataset_donor_id"].astype(str).tolist()
            y = np.arange(len(order))
            lookup = {
                str(d): float(v)
                for d, v in zip(block["dataset_donor_id"], block["median_delta_axis"])
            }
            values = np.array([lookup.get(d, np.nan) for d in order])
            color = KO_COLOR if str(pert).upper() == "KO" else OE_COLOR
            ax.scatter(values, y, s=18, color=color, zorder=3)
            ax.axvline(0, color=ZERO_LINE, linewidth=0.6)
            ax.set_yticks(y, order if j == 0 else [], fontsize=6)
            ax.set_xlabel("embedding-axis shift (Δaxis)")
            ax.set_title(f"{endpoint} · {pert}", fontsize=8, loc="left")
            ax.set_ylim(-0.5, len(order) - 0.5)
            ax.invert_yaxis()
    title = "no KO–OE pairing lines" if unpaired else "same analysis_population only"
    fig.suptitle(title, fontsize=8)
    save_figure(
        fig,
        stem,
        config,
        parameters={"plot": "donor_delta_axis", "ko_oe_unpaired": unpaired, "n_is": "donor-units"},
        input_files=input_files,
    )


def plot_forest(table: pd.DataFrame, config: Mapping[str, Any], stem: Path, input_files: list[Path]) -> None:
    frame = table.copy()
    if "target_symbol" in frame and frame["target_symbol"].nunique() > 1:
        raise SystemExit("Forest requires one target; separate matched controls")
    if "across_donor_median" not in frame.columns:
        summaries = []
        for (endpoint, pert), group in frame.groupby(["endpoint", "perturbation"], sort=True):
            if "eligible" in group:
                group = group.loc[group["eligible"].astype(bool)]
            if "dataset_donor_id" not in group or group["dataset_donor_id"].duplicated().any():
                raise SystemExit("Forest requires one row per biological donor; separate populations")
            values = group["median_delta_axis"].to_numpy(float)
            low, high, _, status = median_order_statistic_interval(values)
            summaries.append(dict(endpoint=endpoint, perturbation=pert,
                                  across_donor_median=np.median(values) if len(values) else np.nan,
                                  median_ci_lower=low, median_ci_upper=high, median_ci_status=status))
        frame = pd.DataFrame(summaries)
    if not {"median_ci_lower", "median_ci_upper", "median_ci_status"}.issubset(frame.columns):
        raise SystemExit("Forest requires donor effects or sign_tests.csv with median CI columns")
    labels = [f"{e} · {p}" for e, p in zip(frame["endpoint"], frame["perturbation"])]
    y = np.arange(len(frame))
    mid = frame["across_donor_median"].to_numpy(float)
    lo = frame["median_ci_lower"].to_numpy(float)
    hi = frame["median_ci_upper"].to_numpy(float)
    finite_ci = np.isfinite(lo) & np.isfinite(hi)
    labels = [label if finite else f"{label} [{status}]"
              for label, finite, status in zip(labels, finite_ci, frame["median_ci_status"])]
    width = figure_size(config, "wide_panel")[0]
    fig, ax = plt.subplots(figsize=(width, max(2.3, 1.6 + .32 * len(frame))), layout="constrained")
    fig.set_layout_engine("constrained", rect=(0, .10, 1, .90))
    if np.any(finite_ci):
        ax.hlines(np.asarray(y[finite_ci], dtype=float), np.asarray(lo[finite_ci], dtype=float),
                  np.asarray(hi[finite_ci], dtype=float), color=KO_COLOR, linewidth=1.2)
    finite_mid = np.isfinite(mid)
    ax.scatter(mid[finite_mid], y[finite_mid], s=24, color=KO_COLOR, zorder=3)
    ax.axvline(0, color=ZERO_LINE, linewidth=0.6)
    ax.set_yticks(y, labels, fontsize=7)
    ax.set_ylim(-.6, len(frame) - .4)
    ax.invert_yaxis()
    ax.set_xlabel("across-donor median embedding-axis shift (Δaxis)")
    ax.set_title("Nominal 95% median interval\nShared-axis dependence not calibrated", loc="left", fontsize=8)
    save_figure(fig, stem, config, parameters={"plot": "primary_forest", "interval": "central_order_statistics",
                "confidence": .95, "inference_status": "NOMINAL_SHARED_AXIS_DEPENDENCE_NOT_CALIBRATED",
                "unbounded_or_missing_intervals": int((~finite_ci).sum())}, input_files=input_files)


def plot_ranks(table: pd.DataFrame, config: Mapping[str, Any], stem: Path, input_files: list[Path]) -> None:
    if table.empty or "directional_rank" not in table.columns:
        fig, ax = plt.subplots(figsize=figure_size(config, "double_panel"), layout="constrained")
        ax.text(0.5, 0.5, "NOT_ESTIMABLE", ha="center", va="center")
        ax.set_axis_off()
        save_figure(fig, stem, config, parameters={"plot": "control_ranks", "status": "NOT_ESTIMABLE"}, input_files=input_files)
        return
    fig, axes = plt.subplots(1, 2, figsize=figure_size(config, "double_panel"), layout="constrained")
    for ax, col, title in (
        (axes[0], "directional_rank", "directional rank"),
        (axes[1], "absolute_rank", "absolute rank"),
    ):
        frame = table.sort_values(col)
        names = frame["target_symbol"].astype(str).tolist()
        vals = frame[col].to_numpy(float)
        is_t = frame["is_target"].astype(bool) if "is_target" in frame.columns else [False] * len(frame)
        colors = [KO_COLOR if flag else FAIL_COLOR for flag in is_t]
        ax.barh(range(len(names)), vals, color=colors, height=0.7)
        ax.set_yticks(range(len(names)), names, fontsize=6)
        ax.invert_yaxis()
        ax.set_xlabel(title)
        ax.set_title(title, loc="left", fontsize=8)
    fig.suptitle("≤10 matched genes; not a permutation p-value", fontsize=8)
    save_figure(fig, stem, config, parameters={"plot": "control_ranks"}, input_files=input_files)


def plot_support(table: pd.DataFrame, config: Mapping[str, Any], stem: Path, input_files: list[Path]) -> None:
    endpoints = list(dict.fromkeys(table["endpoint"].astype(str)))
    perts = list(dict.fromkeys(table["perturbation"].astype(str)))
    matrix = np.full((len(endpoints), len(perts)), np.nan)
    value_col = "across_donor_median" if "across_donor_median" in table.columns else "median_delta_axis"
    for rec in table.itertuples(index=False):
        i = endpoints.index(str(rec.endpoint))
        j = perts.index(str(rec.perturbation))
        matrix[i, j] = float(getattr(rec, value_col))
    fig, ax = plt.subplots(figsize=figure_size(config, "double_panel"), layout="constrained")
    vmax = np.nanmax(np.abs(matrix)) if np.isfinite(matrix).any() else 1.0
    im = ax.imshow(matrix, cmap="coolwarm", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(perts)), perts, fontsize=8)
    ax.set_yticks(range(len(endpoints)), endpoints, fontsize=7)
    ax.set_title("support endpoints are descriptive; no p-value stars", loc="left", fontweight="bold")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Δaxis")
    save_figure(fig, stem, config, parameters={"plot": "support_heatmap"}, input_files=input_files)


def plot_symmetry(table: pd.DataFrame, config: Mapping[str, Any], stem: Path, input_files: list[Path]) -> None:
    frame = table.copy()
    ko = frame.loc[frame["perturbation"].astype(str).str.upper().eq("KO")]
    oe = frame.loc[
        frame["perturbation"].astype(str).str.upper().eq("OE")
        & frame["analysis_population"].astype(str).str.contains("SYMMETRY", case=False, na=False)
    ]
    if oe.empty:
        raise SystemExit("F06_07 needs OE symmetry population; do not plot OE primary here")
    merged = ko.merge(
        oe,
        on=["dataset_donor_id", "endpoint"],
        suffixes=("_ko", "_oe"),
        how="inner",
    )
    endpoints = list(dict.fromkeys(merged["endpoint"].astype(str)))
    fig, axes = plt.subplots(1, max(len(endpoints), 1), figsize=figure_size(config, "double_panel"), layout="constrained", squeeze=False)
    for ax, endpoint in zip(axes[0], endpoints):
        block = merged.loc[merged["endpoint"].astype(str).eq(endpoint)]
        ax.axhline(0, color=ZERO_LINE, linewidth=0.5)
        ax.axvline(0, color=ZERO_LINE, linewidth=0.5)
        ax.scatter(block["median_delta_axis_ko"], block["median_delta_axis_oe"], s=18, color=KO_COLOR)
        ax.set_xlabel("KO median Δaxis")
        ax.set_ylabel("OE-symmetry median Δaxis")
        ax.set_title(str(endpoint), loc="left", fontsize=8)
    fig.suptitle("opposite sign is descriptive; not an extra FDR test", fontsize=8)
    save_figure(fig, stem, config, parameters={"plot": "ko_oe_symmetry"}, input_files=input_files)


def plot_sham(payload: dict[str, Any], config: Mapping[str, Any], stem: Path, input_files: list[Path]) -> None:
    freeze = float(payload.get("sham_max", 1e-6))
    emb = float(payload.get("max_abs_embedding_diff", payload.get("max_abs_embedding", 0.0)))
    cos = float(payload.get("max_cosine", payload.get("max_cosine_distance", 0.0)))
    passed = bool(payload.get("pass", emb <= freeze and cos <= freeze))
    fig, ax = plt.subplots(figsize=figure_size(config, "single_panel"), layout="constrained")
    ax.bar(["|Δemb|", "cosine"], [emb, cos], color=PASS_COLOR if passed else "#D55E00")
    ax.axhline(freeze, color=ZERO_LINE, linewidth=0.8, linestyle="--")
    ax.set_ylabel("max difference vs original CLS")
    ax.set_title("sham re-inference; not a biological result", loc="left", fontweight="bold")
    save_figure(fig, stem, config, parameters={"plot": "sham", "pass": passed, "freeze": freeze}, input_files=input_files)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["observability", "eligibility", "donors", "forest", "ranks", "support", "symmetry", "sham"])
    parser.add_argument("table", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--config", default=None)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--unpaired", dest="paired", action="store_false")
    group.add_argument("--paired", dest="paired", action="store_true")
    parser.set_defaults(paired=False)
    args = parser.parse_args(argv)
    config = _config(args.config)
    args.out.mkdir(parents=True, exist_ok=True)
    inputs = [args.table]
    unpaired = not args.paired
    if args.command == "sham":
        payload = json.loads(args.table.read_text(encoding="utf-8"))
        plot_sham(payload, config, args.out / "F06_08_sham", inputs)
    else:
        table = pd.read_csv(args.table)
        if args.command == "observability":
            plot_observability(table, config, args.out / "F06_01_token_observability", inputs)
        elif args.command == "eligibility":
            plot_eligibility(table, config, args.out / "F06_02_donor_eligibility", inputs)
        elif args.command == "donors":
            plot_donors(table, config, args.out / "F06_03_donor_delta_axis", unpaired=unpaired, input_files=inputs)
        elif args.command == "forest":
            plot_forest(table, config, args.out / "F06_04_primary_forest", inputs)
        elif args.command == "ranks":
            plot_ranks(table, config, args.out / "F06_05_control_ranks", inputs)
        elif args.command == "support":
            plot_support(table, config, args.out / "F06_06_support_heatmap", inputs)
        elif args.command == "symmetry":
            plot_symmetry(table, config, args.out / "F06_07_ko_oe_symmetry", inputs)
    print(f"wrote F06 panel under {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


