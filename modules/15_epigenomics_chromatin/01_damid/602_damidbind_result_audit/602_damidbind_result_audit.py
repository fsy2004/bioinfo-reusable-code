#!/usr/bin/env python3
"""Audit an exported damidBind differential-binding table."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
MODULES = ROOT.parents[2]
sys.path.insert(0, str(MODULES / "_framework"))
from output_audit import require_columns, save_volcano, write_summary


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=ROOT / "example_data" / "damidbind_results.csv")
    parser.add_argument("--outdir", type=Path, default=ROOT / "results")
    parser.add_argument("--fdr", type=float, default=0.05)
    return parser.parse_args()


def audit(path: Path, fdr: float):
    data = pd.read_csv(path)
    require_columns(data.columns, ["locus", "log2FoldChange", "pvalue", "padj"])
    if data.empty:
        raise ValueError("differential-binding table is empty")
    for column in ("log2FoldChange", "pvalue", "padj"):
        data[column] = pd.to_numeric(data[column], errors="raise")
    if not data["padj"].between(0, 1).all():
        raise ValueError("padj must be between 0 and 1")
    sig = data["padj"] <= fdr
    summary = {"tested_loci": int(len(data)), "significant_loci": int(sig.sum()), "gain_loci": int((sig & (data["log2FoldChange"] > 0)).sum()), "loss_loci": int((sig & (data["log2FoldChange"] < 0)).sum()), "fdr_threshold": float(fdr)}
    return summary, data, sig


def main() -> int:
    args = parse_args()
    try:
        summary, data, sig = audit(args.input, args.fdr)
        write_summary(args.outdir, "602", summary)
        save_volcano(args.outdir, "602", data["log2FoldChange"].to_numpy(), data["padj"].to_numpy(), sig.to_numpy(), "Log2 fold change")
        return 0
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
