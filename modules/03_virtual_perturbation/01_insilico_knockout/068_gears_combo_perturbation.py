"""GEARS single-cell perturbation prediction wrapper.

Example:
python 068_gears_combo_perturbation.py --adata input.h5ad --dataset-name mydata \
  --perturbations STAT3,STAT3+JUN --output-dir results/gears --train
"""

from __future__ import annotations

import argparse
import json
import pickle
from pathlib import Path


def parse_perturbations(text: str) -> list[list[str]]:
    items = []
    for part in text.split(","):
        genes = [g.strip() for g in part.split("+") if g.strip()]
        if genes:
            items.append(genes)
    return items


def main() -> None:
    parser = argparse.ArgumentParser(description="Run GEARS on perturb-seq style AnnData.")
    parser.add_argument("--adata", help="Custom AnnData h5ad file. Requires condition/cell_type/gene_name fields.")
    parser.add_argument("--data-dir", default="data/gears", help="GEARS data directory.")
    parser.add_argument("--dataset-name", default="custom", help="GEARS dataset name or built-in dataset name.")
    parser.add_argument("--condition-key", default="condition", help="adata.obs column for perturbation condition.")
    parser.add_argument("--celltype-key", default="cell_type", help="adata.obs column for cell type.")
    parser.add_argument("--gene-name-key", default="gene_name", help="adata.var column for gene names.")
    parser.add_argument("--perturbations", required=True, help="Comma-separated genes or gene pairs, e.g. STAT3,JUN+FOS.")
    parser.add_argument("--output-dir", default="results/gears", help="Output directory.")
    parser.add_argument("--device", default="cuda", help="cuda, cuda:0, or cpu.")
    parser.add_argument("--epochs", type=int, default=20)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--test-batch-size", type=int, default=128)
    parser.add_argument("--train", action="store_true", help="Train a GEARS model before prediction.")
    parser.add_argument("--model-dir", default="", help="Existing model directory to load, or save path after training.")
    args = parser.parse_args()

    try:
        import scanpy as sc
        from gears import GEARS, PertData
    except ImportError as exc:
        raise SystemExit(
            "Missing dependency. Install GEARS first, for example: pip install cell-gears"
        ) from exc

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    perturbations = parse_perturbations(args.perturbations)

    pert_data = PertData(args.data_dir)
    if args.adata:
        adata = sc.read_h5ad(args.adata)
        if args.gene_name_key not in adata.var:
            adata.var[args.gene_name_key] = adata.var_names.astype(str)
        adata.var["gene_name"] = adata.var[args.gene_name_key].astype(str)
        adata.obs["condition"] = adata.obs[args.condition_key].astype(str)
        adata.obs["cell_type"] = adata.obs[args.celltype_key].astype(str)
        pert_data.new_data_process(dataset_name=args.dataset_name, adata=adata)
        pert_data.load(data_path=str(Path(args.data_dir) / args.dataset_name))
    else:
        pert_data.load(data_name=args.dataset_name)

    pert_data.prepare_split(split="simulation", seed=1)
    pert_data.get_dataloader(batch_size=args.batch_size, test_batch_size=args.test_batch_size)

    model = GEARS(pert_data, device=args.device)
    if args.train:
        model.model_initialize(hidden_size=64)
        model.train(epochs=args.epochs)
        save_dir = args.model_dir or str(out / "gears_model")
        model.save_model(save_dir)
    elif args.model_dir:
        model.load_pretrained(args.model_dir)
    else:
        raise SystemExit("Provide --train or --model-dir for prediction.")

    pred = model.predict(perturbations)
    with open(out / "gears_prediction.pkl", "wb") as fh:
        pickle.dump(pred, fh)

    summary = {
        "dataset": args.dataset_name,
        "perturbations": perturbations,
        "prediction_pickle": str(out / "gears_prediction.pkl"),
    }
    (out / "gears_run_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"Done. Results saved to: {out.resolve()}")


if __name__ == "__main__":
    main()

