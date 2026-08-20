#!/usr/bin/env python3
"""Audit OrthoFinder Orthogroups.tsv output."""
from __future__ import annotations

import argparse
import csv
import shlex
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
MODULES = ROOT.parents[2]
sys.path.insert(0, str(MODULES / "_framework"))
from output_audit import save_bar, write_summary


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=ROOT / "example_data" / "Orthogroups.tsv")
    parser.add_argument("--outdir", type=Path, default=ROOT / "results")
    parser.add_argument("--print-upstream-command", action="store_true")
    parser.add_argument("--proteomes")
    parser.add_argument("--search-threads", type=int, default=8)
    parser.add_argument("--analysis-threads", type=int, default=4)
    return parser.parse_args()


def audit(path: Path):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or reader.fieldnames[0] != "Orthogroup" or len(reader.fieldnames) < 2:
            raise ValueError("Orthogroups table requires Orthogroup and at least one species column")
        species = reader.fieldnames[1:]
        groups = total_genes = complete = single = 0
        for row in reader:
            groups += 1
            counts = []
            for name in species:
                cell = (row.get(name) or "").strip()
                count = 0 if not cell else len([gene for gene in cell.split(",") if gene.strip()])
                counts.append(count)
                total_genes += count
            is_complete = all(count > 0 for count in counts)
            complete += int(is_complete)
            single += int(is_complete and all(count == 1 for count in counts))
    if groups == 0:
        raise ValueError("no orthogroups found")
    return {"orthogroups": groups, "species": len(species), "total_gene_assignments": total_genes, "complete_orthogroups": complete, "single_copy_complete_orthogroups": single}


def main() -> int:
    args = parse_args()
    try:
        if args.print_upstream_command:
            if not args.proteomes:
                raise ValueError("OrthoFinder command requires proteomes directory")
            print(f"orthofinder -f {shlex.quote(args.proteomes)} -t {args.search_threads} -a {args.analysis_threads}")
            return 0
        summary = audit(args.input)
        write_summary(args.outdir, "600", summary)
        save_bar(args.outdir, "600", ["All groups", "Complete", "Single-copy"], [summary["orthogroups"], summary["complete_orthogroups"], summary["single_copy_complete_orthogroups"]], "Orthogroup count")
        return 0
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
