#!/usr/bin/env python3
"""Audit EviAnn GFF3 output; optionally print its verified upstream CLI."""
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
    parser.add_argument("--input", type=Path, default=ROOT / "example_data" / "eviann_annotation.gff3")
    parser.add_argument("--outdir", type=Path, default=ROOT / "results")
    parser.add_argument("--print-upstream-command", action="store_true")
    parser.add_argument("--genome")
    parser.add_argument("--reads-list")
    parser.add_argument("--transcripts")
    parser.add_argument("--proteins")
    parser.add_argument("--threads", type=int, default=24)
    return parser.parse_args()


def parse_attributes(raw: str):
    return {piece.split("=", 1)[0]: piece.split("=", 1)[1] for piece in raw.split(";") if "=" in piece}


def upstream_command(args) -> str:
    if not args.genome or not (args.reads_list or args.transcripts):
        raise ValueError("EviAnn command requires genome and reads-list or transcripts")
    parts = ["eviann.sh", "-t", str(args.threads), "-g", shlex.quote(args.genome)]
    if args.reads_list:
        parts.extend(["-r", shlex.quote(args.reads_list)])
    if args.transcripts:
        parts.extend(["-e", shlex.quote(args.transcripts)])
    if args.proteins:
        parts.extend(["-p", shlex.quote(args.proteins)])
    return " ".join(parts)


def audit(path: Path):
    summary = {"gene_features": 0, "transcript_features": 0, "exon_features": 0, "cds_features": 0, "complete_transcripts": 0, "transcripts_with_evidence": 0, "pseudo_transcripts": 0}
    with path.open(encoding="utf-8") as handle:
        for lineno, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"GFF3 line {lineno} must have 9 tab-separated fields")
            feature = fields[2]
            attrs = parse_attributes(fields[8])
            if feature == "gene":
                summary["gene_features"] += 1
            elif feature in {"mRNA", "transcript", "ncRNA"}:
                summary["transcript_features"] += 1
                summary["complete_transcripts"] += int(attrs.get("Evidence", "").lower() == "complete")
                summary["transcripts_with_evidence"] += int(any(key in attrs for key in ("EvidenceProteinID", "EvidenceTranscriptID")))
                summary["pseudo_transcripts"] += int(attrs.get("pseudo", "").lower() == "true")
            elif feature == "exon":
                summary["exon_features"] += 1
            elif feature == "CDS":
                summary["cds_features"] += 1
    if summary["gene_features"] == 0:
        raise ValueError("no gene features found")
    return summary


def main() -> int:
    args = parse_args()
    try:
        if args.print_upstream_command:
            print(upstream_command(args))
            return 0
        summary = audit(args.input)
        write_summary(args.outdir, "599", summary)
        labels = ["Genes", "Transcripts", "Exons", "CDS", "Complete"]
        values = [summary["gene_features"], summary["transcript_features"], summary["exon_features"], summary["cds_features"], summary["complete_transcripts"]]
        save_bar(args.outdir, "599", labels, values, "Feature count")
        return 0
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
