#!/usr/bin/env python
"""AutoDock Vina / GROMACS / gmx_MMPBSA / MDAnalysis pipeline wrapper."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
from pathlib import Path


def run(cmd: list[str], log_path: Path) -> None:
    with log_path.open("a", encoding="utf-8") as log:
        log.write("$ " + " ".join(cmd) + "\n")
        subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--outdir", required=True)
    p.add_argument("--receptor_pdbqt")
    p.add_argument("--ligand_pdbqt")
    p.add_argument("--vina_config")
    p.add_argument("--trajectory")
    p.add_argument("--topology")
    p.add_argument("--gmx_mmpbsa_input")
    p.add_argument("--run_vina", action="store_true")
    p.add_argument("--run_mmpbsa", action="store_true")
    p.add_argument("--run_mdanalysis", action="store_true")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    log_path = outdir / "docking_md_pipeline.log"
    commands: dict[str, list[str]] = {}

    if args.run_vina:
        if not all([args.receptor_pdbqt, args.ligand_pdbqt, args.vina_config]):
            raise SystemExit("--run_vina requires --receptor_pdbqt --ligand_pdbqt --vina_config")
        if shutil.which("vina") is None:
            raise SystemExit("AutoDock Vina executable 'vina' not found in PATH.")
        out_pose = outdir / "vina_out.pdbqt"
        cmd = [
            "vina", "--receptor", args.receptor_pdbqt, "--ligand", args.ligand_pdbqt,
            "--config", args.vina_config, "--out", str(out_pose),
        ]
        commands["vina"] = cmd
        run(cmd, log_path)

    if args.run_mmpbsa:
        if not all([args.gmx_mmpbsa_input, args.trajectory, args.topology]):
            raise SystemExit("--run_mmpbsa requires --gmx_mmpbsa_input --trajectory --topology")
        if shutil.which("gmx_MMPBSA") is None:
            raise SystemExit("gmx_MMPBSA executable not found in PATH.")
        cmd = [
            "gmx_MMPBSA", "-i", args.gmx_mmpbsa_input,
            "-ct", args.trajectory, "-cp", args.topology,
            "-o", str(outdir / "FINAL_RESULTS_MMPBSA.dat"),
        ]
        commands["gmx_MMPBSA"] = cmd
        run(cmd, log_path)

    if args.run_mdanalysis:
        if not all([args.trajectory, args.topology]):
            raise SystemExit("--run_mdanalysis requires --trajectory --topology")
        import MDAnalysis as mda
        import pandas as pd
        from MDAnalysis.analysis import rms

        u = mda.Universe(args.topology, args.trajectory)
        r = rms.RMSD(u, select="backbone").run()
        pd.DataFrame(r.results.rmsd, columns=["frame", "time", "rmsd"]).to_csv(outdir / "mdanalysis_rmsd.csv", index=False)

    (outdir / "docking_md_commands.json").write_text(
        json.dumps(commands, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
