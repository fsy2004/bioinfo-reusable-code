#!/usr/bin/env python
"""Squidiff/PerturbDiff external runner.

This wrapper records a reproducible command and calls a user-specified training or
sampling script from the local external source checkout.
"""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--mode", choices=["train", "sample", "predict"], required=True)
    p.add_argument("--script", required=True, help="External Squidiff/PerturbDiff Python entry script.")
    p.add_argument("--input_h5ad")
    p.add_argument("--config")
    p.add_argument("--checkpoint")
    p.add_argument("--outdir", required=True)
    p.add_argument("--extra_args", default="", help="Extra command string appended as-is.")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    cmd = ["python", args.script]
    if args.input_h5ad:
        cmd += ["--input_h5ad", args.input_h5ad]
    if args.config:
        cmd += ["--config", args.config]
    if args.checkpoint:
        cmd += ["--checkpoint", args.checkpoint]
    cmd += ["--outdir", str(outdir)]
    if args.extra_args:
        cmd += args.extra_args.split()

    (outdir / "squidiff_runner_command.json").write_text(
        json.dumps({"mode": args.mode, "command": cmd}, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    with (outdir / "squidiff_runner.log").open("w", encoding="utf-8") as log:
        log.write("$ " + " ".join(cmd) + "\n")
        subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT)


if __name__ == "__main__":
    main()
