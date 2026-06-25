"""CellOracle GRN-based in silico gene perturbation wrapper.

This script expects a prepared CellOracle object saved by pickle.
It runs one or more gene knockdown simulations and exports the updated object.
"""

from __future__ import annotations

import argparse
import pickle
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description="Run CellOracle in silico perturbation from a prepared Oracle object.")
    parser.add_argument("--oracle-pkl", required=True, help="Prepared CellOracle Oracle object pickle.")
    parser.add_argument("--genes", required=True, help="Comma-separated genes for knockdown/knockout simulation.")
    parser.add_argument("--output-dir", default="results/celloracle", help="Output directory.")
    parser.add_argument("--knockdown-value", type=float, default=0.0, help="Expression value after perturbation.")
    parser.add_argument("--n-propagation", type=int, default=3)
    args = parser.parse_args()

    try:
        import celloracle  # noqa: F401
    except ImportError as exc:
        raise SystemExit("Missing dependency. Install with: pip install celloracle") from exc

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    genes = [g.strip() for g in args.genes.split(",") if g.strip()]

    with open(args.oracle_pkl, "rb") as fh:
        oracle = pickle.load(fh)

    for gene in genes:
        gene_dir = out / gene
        gene_dir.mkdir(parents=True, exist_ok=True)
        perturbed = oracle.copy()
        perturbed.simulate_shift(
            perturb_condition={gene: args.knockdown_value},
            n_propagation=args.n_propagation,
        )
        perturbed.estimate_transition_prob(n_neighbors=200, knn_random=True, sampled_fraction=1)
        perturbed.calculate_embedding_shift(sigma_corr=0.05)
        with open(gene_dir / f"{gene}_celloracle_perturbed.pkl", "wb") as fh:
            pickle.dump(perturbed, fh)

    print(f"Done. CellOracle perturbation objects saved to: {out.resolve()}")


if __name__ == "__main__":
    main()

