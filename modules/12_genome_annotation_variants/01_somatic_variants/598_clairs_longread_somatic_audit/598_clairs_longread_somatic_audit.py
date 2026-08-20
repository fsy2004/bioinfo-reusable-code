#!/usr/bin/env python3
"""Audit ClairS-compatible VCF output; optionally print the verified upstream CLI."""
from __future__ import annotations

import argparse
import shlex
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
MODULES = ROOT.parents[2]
sys.path.insert(0, str(MODULES / "_framework"))
from output_audit import save_bar, write_summary


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=ROOT / "example_data" / "clairs_calls.vcf")
    parser.add_argument("--outdir", type=Path, default=ROOT / "results")
    parser.add_argument("--print-upstream-command", action="store_true")
    parser.add_argument("--tumor-bam")
    parser.add_argument("--normal-bam")
    parser.add_argument("--reference")
    parser.add_argument("--platform")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--upstream-output")
    return parser.parse_args()


def quoted(value: str) -> str:
    return shlex.quote(value)


def upstream_command(args) -> str:
    required = [args.tumor_bam, args.normal_bam, args.reference, args.platform, args.upstream_output]
    if any(value is None for value in required):
        raise ValueError("ClairS command requires tumor BAM, normal BAM, reference, platform, and output")
    return " ".join(["run_clairs", "-T", quoted(args.tumor_bam), "-N", quoted(args.normal_bam), "-R", quoted(args.reference), "-o", quoted(args.upstream_output), "-t", str(args.threads), "-p", quoted(args.platform)])


def audit(path: Path):
    counts = {"total_variants": 0, "pass_variants": 0, "filtered_variants": 0, "snv_variants": 0, "indel_variants": 0, "other_variants": 0}
    with path.open(encoding="utf-8") as handle:
        for lineno, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise ValueError(f"VCF line {lineno} has fewer than 8 fields")
            ref, alt, filt = fields[3], fields[4], fields[6]
            counts["total_variants"] += 1
            counts["pass_variants" if filt in {"PASS", "."} else "filtered_variants"] += 1
            alts = alt.split(",")
            if len(ref) == 1 and all(len(a) == 1 for a in alts):
                counts["snv_variants"] += 1
            elif any(len(a) != len(ref) for a in alts):
                counts["indel_variants"] += 1
            else:
                counts["other_variants"] += 1
    if counts["total_variants"] == 0:
        raise ValueError("no variant records found")
    return counts


def main() -> int:
    args = parse_args()
    try:
        if args.print_upstream_command:
            print(upstream_command(args))
            return 0
        summary = audit(args.input)
        write_summary(args.outdir, "598", summary)
        save_bar(args.outdir, "598", ["PASS", "Filtered", "SNV", "Indel", "Other"], [summary["pass_variants"], summary["filtered_variants"], summary["snv_variants"], summary["indel_variants"], summary["other_variants"]], "Variant count")
        return 0
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
