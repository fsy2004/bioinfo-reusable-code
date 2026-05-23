#!/usr/bin/env python
"""pySCENIC command wrapper.

Expected external commands:
  pyscenic grn
  pyscenic ctx
  pyscenic aucell
"""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path


def run(cmd: list[str], log_path: Path) -> None:
    with log_path.open("a", encoding="utf-8") as log:
        log.write("$ " + " ".join(cmd) + "\n")
        subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--expression", required=True, help="Expression matrix or loom file.")
    p.add_argument("--tf_list", required=True)
    p.add_argument("--ranking_db", required=True)
    p.add_argument("--motif_annotations", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--prefix", default="pyscenic")
    p.add_argument("--workers", default="8")
    p.add_argument("--mode", choices=["all", "grn", "ctx", "aucell"], default="all")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    log_path = outdir / f"{args.prefix}.log"
    adj = outdir / f"{args.prefix}_adjacencies.tsv"
    regulons = outdir / f"{args.prefix}_regulons.csv"
    auc = outdir / f"{args.prefix}_aucell.loom"

    commands = {}
    if args.mode in {"all", "grn"}:
        commands["grn"] = [
            "pyscenic", "grn", args.expression, args.tf_list,
            "-o", str(adj), "--num_workers", args.workers,
        ]
    if args.mode in {"all", "ctx"}:
        commands["ctx"] = [
            "pyscenic", "ctx", str(adj), args.ranking_db,
            "--annotations_fname", args.motif_annotations,
            "--expression_mtx_fname", args.expression,
            "-o", str(regulons), "--num_workers", args.workers,
        ]
    if args.mode in {"all", "aucell"}:
        commands["aucell"] = [
            "pyscenic", "aucell", args.expression, str(regulons),
            "-o", str(auc), "--num_workers", args.workers,
        ]

    for name, cmd in commands.items():
        run(cmd, log_path)

    (outdir / f"{args.prefix}_commands.json").write_text(
        json.dumps(commands, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
