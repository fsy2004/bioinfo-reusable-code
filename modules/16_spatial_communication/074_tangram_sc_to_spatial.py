"""Tangram single-cell to spatial mapping wrapper."""

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description="Map single-cell data to spatial transcriptomics with Tangram.")
    parser.add_argument("--single-cell", required=True, help="Single-cell h5ad.")
    parser.add_argument("--spatial", required=True, help="Spatial h5ad.")
    parser.add_argument("--annotation", default="", help="Optional scRNA obs column for cell-type projection.")
    parser.add_argument("--output-dir", default="results/tangram")
    parser.add_argument("--device", default="cpu")
    args = parser.parse_args()

    try:
        import scanpy as sc
        import tangram as tg
    except ImportError as exc:
        raise SystemExit("Install dependencies first: pip install tangram-sc scanpy") from exc

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    ad_sc = sc.read_h5ad(args.single_cell)
    ad_sp = sc.read_h5ad(args.spatial)
    tg.pp_adatas(ad_sc, ad_sp, genes=None)
    ad_map = tg.map_cells_to_space(ad_sc, ad_sp, device=args.device)
    ad_map.write_h5ad(out / "tangram_cell_to_space_map.h5ad")

    if args.annotation:
        tg.project_cell_annotations(ad_map, ad_sp, annotation=args.annotation)
        if f"tangram_ct_pred" in ad_sp.obsm:
            ad_sp.obsm["tangram_ct_pred"].to_csv(out / "tangram_celltype_projection.csv")
    ad_sp.write_h5ad(out / "tangram_spatial_projected.h5ad")
    print(f"Done. Tangram results saved to: {out.resolve()}")


if __name__ == "__main__":
    main()

