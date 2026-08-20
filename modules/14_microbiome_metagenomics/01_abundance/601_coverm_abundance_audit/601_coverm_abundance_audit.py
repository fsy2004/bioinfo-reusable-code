#!/usr/bin/env python3
"""Audit long-form CoverM genome relative-abundance output."""
from __future__ import annotations

import argparse
import shlex
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
MODULES = ROOT.parents[2]
sys.path.insert(0, str(MODULES / "_framework"))
from output_audit import require_columns, save_bar, write_summary


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=ROOT / "example_data" / "coverm_abundance.csv")
    parser.add_argument("--outdir", type=Path, default=ROOT / "results")
    parser.add_argument("--print-upstream-command", action="store_true")
    parser.add_argument("--bam")
    parser.add_argument("--genome-definition")
    return parser.parse_args()


def audit(path: Path):
    data = pd.read_csv(path)
    require_columns(data.columns, ["sample", "genome", "relative_abundance"])
    if data.empty:
        raise ValueError("abundance table is empty")
    data["relative_abundance"] = pd.to_numeric(data["relative_abundance"], errors="raise")
    if (data["relative_abundance"] < 0).any():
        raise ValueError("relative_abundance must be nonnegative")
    totals = data.groupby("sample", sort=True)["relative_abundance"].sum()
    summary = {"total_rows": int(len(data)), "samples": int(data["sample"].nunique()), "genomes": int(data["genome"].nunique()), "samples_summing_to_one": int(((totals - 1.0).abs() <= 1e-6).sum()), "max_relative_abundance": float(data["relative_abundance"].max())}
    return summary, totals


def main() -> int:
    args = parse_args()
    try:
        if args.print_upstream_command:
            if not args.bam or not args.genome_definition:
                raise ValueError("CoverM command requires BAM and genome definition")
            print(f"coverm genome --bam-files {shlex.quote(args.bam)} --genome-definition {shlex.quote(args.genome_definition)} --methods relative_abundance")
            return 0
        summary, totals = audit(args.input)
        write_summary(args.outdir, "601", summary)
        save_bar(args.outdir, "601", list(totals.index), list(totals.values), "Relative-abundance sum")
        return 0
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
