"""chemCPA orchestration wrapper.

This wrapper keeps chemCPA as an external dependency and records each run in
this repository's standard result folder.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description="Run chemCPA training/prediction from an external checkout.")
    parser.add_argument("--chemcpa-repo", required=True, help="Path to cloned theislab/chemCPA repository.")
    parser.add_argument("--config", required=True, help="Hydra/seml YAML config for chemCPA.")
    parser.add_argument("--output-dir", default="results/chemcpa", help="Output directory for logs and run metadata.")
    parser.add_argument("--python", default=sys.executable, help="Python executable in the chemCPA environment.")
    parser.add_argument("--dry-run", action="store_true", help="Only write the command without executing it.")
    args = parser.parse_args()

    repo = Path(args.chemcpa_repo).resolve()
    train_py = repo / "chemCPA" / "train_hydra.py"
    if not train_py.exists():
        raise SystemExit(f"Cannot find chemCPA training entry: {train_py}")

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    cmd = [args.python, str(train_py), f"--config-name={Path(args.config).name}"]

    metadata = {
        "chemcpa_repo": str(repo),
        "config": str(Path(args.config).resolve()),
        "command": cmd,
        "note": "Run this in an environment where chemCPA dependencies are installed.",
    }
    (out / "chemcpa_run_command.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    if args.dry_run:
        print("Dry run command:", " ".join(cmd))
        return

    with open(out / "chemcpa_stdout.log", "w", encoding="utf-8") as stdout, open(
        out / "chemcpa_stderr.log", "w", encoding="utf-8"
    ) as stderr:
        subprocess.run(cmd, cwd=repo, check=True, stdout=stdout, stderr=stderr)

    print(f"Done. Logs saved to: {out.resolve()}")


if __name__ == "__main__":
    main()

