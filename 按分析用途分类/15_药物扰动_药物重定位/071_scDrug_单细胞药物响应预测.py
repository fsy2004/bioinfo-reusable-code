"""scDrug pipeline wrapper.

Runs selected scDrug steps from an external scDrug checkout:
1. single-cell preprocessing
2. cluster-level drug response prediction
3. optional treatment selection
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path


def run(cmd: list[str], cwd: Path, log_file: Path) -> None:
    with open(log_file, "w", encoding="utf-8") as log:
        log.write("COMMAND: " + " ".join(cmd) + "\n\n")
        subprocess.run(cmd, cwd=cwd, check=True, stdout=log, stderr=subprocess.STDOUT)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run scDrug single-cell drug response workflow.")
    parser.add_argument("--scdrug-repo", required=True, help="Path to cloned ailabstw/scDrug repository.")
    parser.add_argument("--input", required=True, help="Input 10x directory, CSV, or h5ad.")
    parser.add_argument("--format", default="h5ad", choices=["10x", "csv", "h5ad"])
    parser.add_argument("--metadata", default="", help="Optional metadata CSV.")
    parser.add_argument("--batch", default="", help="Batch column in metadata.")
    parser.add_argument("--clusters", default="All", help="Clusters for drug response prediction.")
    parser.add_argument("--model", default="GDSC", choices=["GDSC", "PRISM"])
    parser.add_argument("--output-dir", default="results/scdrug")
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--skip-preprocess", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    repo = Path(args.scdrug_repo).resolve()
    out = Path(args.output_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)

    single_cell = repo / "single_cell_analysis.py"
    drug_pred = repo / "drug_response_prediction.py"
    for entry in [drug_pred] + ([] if args.skip_preprocess else [single_cell]):
        if not entry.exists():
            raise SystemExit(f"Cannot find scDrug entry script: {entry}")

    commands = []
    if not args.skip_preprocess:
        cmd = [args.python, str(single_cell), "--input", args.input, "--format", args.format, "--output", str(out)]
        if args.metadata:
            cmd += ["--metadata", args.metadata]
        if args.batch:
            cmd += ["--batch", args.batch]
        commands.append(("single_cell_analysis", cmd))

    h5ad_input = str(out / "scanpyobj.h5ad") if not args.skip_preprocess else args.input
    commands.append(
        (
            "drug_response_prediction",
            [
                args.python,
                str(drug_pred),
                "--input",
                h5ad_input,
                "--output",
                str(out),
                "--clusters",
                args.clusters,
                "--model",
                args.model,
            ],
        )
    )

    (out / "scdrug_commands.json").write_text(
        json.dumps({name: cmd for name, cmd in commands}, indent=2), encoding="utf-8"
    )
    if args.dry_run:
        for name, cmd in commands:
            print(name + ":", " ".join(cmd))
        return

    for name, cmd in commands:
        run(cmd, repo, out / f"{name}.log")

    print(f"Done. scDrug outputs saved to: {out}")


if __name__ == "__main__":
    main()

