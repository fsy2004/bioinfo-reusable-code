"""Small shared helpers for source-audited output modules."""
from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def require_columns(columns: Sequence[str], required: Sequence[str]) -> None:
    missing = [name for name in required if name not in columns]
    if missing:
        raise ValueError("missing required column(s): " + ", ".join(missing))


def write_summary(outdir: Path, module_id: str, summary: Mapping[str, object]) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    with (outdir / f"{module_id}_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(dict(summary), handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
    with (outdir / f"{module_id}_summary.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["metric", "value"])
        writer.writerows(summary.items())


def save_bar(outdir: Path, module_id: str, labels: Sequence[str], values: Sequence[float], ylabel: str) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "Arial", "font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(figsize=(5.2, 3.4))
    ax.bar(labels, values, color="#3C5488", width=0.68)
    ax.set_ylabel(ylabel)
    ax.tick_params(axis="x", rotation=25)
    fig.tight_layout()
    for suffix, kwargs in (("png", {"dpi": 300}), ("pdf", {})):
        fig.savefig(outdir / f"{module_id}_audit.{suffix}", bbox_inches="tight", **kwargs)
    plt.close(fig)


def save_volcano(outdir: Path, module_id: str, effects, padj, significant, xlabel: str) -> None:
    import numpy as np

    outdir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "Arial", "font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    y = -np.log10(np.clip(padj, 1e-300, 1.0))
    colors = np.where(significant, "#E64B35", "#8491B4")
    fig, ax = plt.subplots(figsize=(4.8, 3.8))
    ax.scatter(effects, y, c=colors, s=28, alpha=0.85, edgecolors="none")
    ax.axhline(-np.log10(0.05), color="#4D4D4D", linestyle="--", linewidth=0.8)
    ax.axvline(0, color="#4D4D4D", linewidth=0.7)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("-log10(adjusted P)")
    fig.tight_layout()
    for suffix, kwargs in (("png", {"dpi": 300}), ("pdf", {})):
        fig.savefig(outdir / f"{module_id}_audit.{suffix}", bbox_inches="tight", **kwargs)
    plt.close(fig)
