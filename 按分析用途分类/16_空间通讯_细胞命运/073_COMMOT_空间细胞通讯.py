"""COMMOT spatial cell-cell communication wrapper."""

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description="Run COMMOT ligand-receptor spatial communication.")
    parser.add_argument("--adata", required=True, help="Spatial AnnData h5ad.")
    parser.add_argument("--species", default="human", choices=["human", "mouse"])
    parser.add_argument("--database", default="CellChat", help="Ligand-receptor database name supported by COMMOT.")
    parser.add_argument("--distance-threshold", type=float, default=200)
    parser.add_argument("--output-dir", default="results/commot")
    args = parser.parse_args()

    try:
        import commot as ct
        import scanpy as sc
    except ImportError as exc:
        raise SystemExit("Install dependencies first: pip install commot scanpy") from exc

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(args.adata)

    if "spatial" not in adata.obsm:
        raise SystemExit("COMMOT requires spatial coordinates in adata.obsm['spatial'].")

    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    lr = ct.pp.ligand_receptor_database(database=args.database, species=args.species)
    ct.tl.spatial_communication(
        adata,
        database_name=args.database,
        df_ligrec=lr,
        dis_thr=args.distance_threshold,
        heteromeric=True,
    )

    adata.write_h5ad(out / "commot_spatial_communication.h5ad")
    for key in ["commot-" + args.database + "-sum-sender", "commot-" + args.database + "-sum-receiver"]:
        if key in adata.obsm:
            adata.obsm[key].to_csv(out / f"{key}.csv")
    print(f"Done. COMMOT results saved to: {out.resolve()}")


if __name__ == "__main__":
    main()

