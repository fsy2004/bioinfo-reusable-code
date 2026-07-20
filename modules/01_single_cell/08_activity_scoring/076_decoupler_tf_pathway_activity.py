"""decoupler TF/pathway activity scoring wrapper."""

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description="Infer TF or pathway activity scores with decoupler.")
    parser.add_argument("--adata", required=True, help="Input h5ad.")
    parser.add_argument("--resource", default="collectri", choices=["collectri", "progeny"])
    parser.add_argument("--organism", default="human", choices=["human", "mouse"])
    parser.add_argument("--output-dir", default="results/decoupler")
    args = parser.parse_args()

    try:
        import decoupler as dc
        import scanpy as sc
    except ImportError as exc:
        raise SystemExit("Install dependencies first: pip install decoupler scanpy") from exc

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(args.adata)

    if args.resource == "collectri":
        net = dc.get_collectri(organism=args.organism, split_complexes=False)
        source, target, weight = "source", "target", "weight"
    else:
        net = dc.get_progeny(organism=args.organism, top=500)
        source, target, weight = "source", "target", "weight"

    if hasattr(dc, "run_mlm"):
        dc.run_mlm(adata, net, source=source, target=target, weight=weight, use_raw=False)
        key = "mlm_estimate"
        scores = adata.obsm[key] if key in adata.obsm else None
    else:
        result = dc.mt.mlm(adata, net, source=source, target=target, weight=weight)
        scores = result[0] if isinstance(result, tuple) else result

    adata.write_h5ad(out / f"decoupler_{args.resource}_activity.h5ad")
    if scores is not None:
        scores.to_csv(out / f"decoupler_{args.resource}_activity_scores.csv")
    print(f"Done. decoupler results saved to: {out.resolve()}")


if __name__ == "__main__":
    main()

