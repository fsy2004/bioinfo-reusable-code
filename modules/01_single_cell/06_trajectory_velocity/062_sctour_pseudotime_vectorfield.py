# ==========================================================================
# 脚本名     : scTour拟时序向量场_Python参考.py
# 分类       : 08_singlecell_spatial_trajectory
# 项目来源   : 从压缩包 491.SCTOUR.rar 整理
# 原始文件   : 491.SCTOUR\06_run_official_tutorial_exact.py
# 用途       : Python 脚本，复现 scTour 官方教程：训练 scTour 模型，推断 pseudotime、潜在空间和向量场，并输出UMAP/向量场图。
# 结果图     : celltype UMAP；batch UMAP；pseudotime UMAP；vector field向量场2x2图；pseudotime表；训练日志；h5ad结果
# 非肿瘤消化适配: 适合但不是R脚本。可作为单细胞轨迹/向量场高级图参考；R主线可暂时只记录，不强行改写。
# 主要 Python 包  : Python: scanpy; sctour; matplotlib; argparse; json; pathlib
# 整理日期   : 2026-05-13
# 备注       : 保留bioinfo-reusable-code逻辑，仅添加统一说明头；运行前请把 workDir/setwd 和输入文件名改成当前项目路径。
# ==========================================================================
# -*- coding: utf-8 -*-
"""Reproduce the scTour official tutorial figure as closely as possible."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import scanpy as sc
import sctour as sct


ROOT = Path(__file__).resolve().parent


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the exact official scTour tutorial workflow."
    )
    parser.add_argument(
        "--input",
        default=str(ROOT / "EX_development_human_cortex_10X.h5ad"),
    )
    parser.add_argument(
        "--output-dir",
        default=str(ROOT / "sctour_results_official_tutorial_exact"),
    )
    parser.add_argument("--n-top-genes", type=int, default=1000)
    parser.add_argument("--nepoch", type=int, default=None)
    parser.add_argument("--percent", type=float, default=None)
    parser.add_argument("--force-cpu", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_path = Path(args.input)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    print(f"Reading: {input_path}", flush=True)
    adata = sc.read_h5ad(input_path)
    print(f"Data shape: {adata.shape}", flush=True)

    # The following block intentionally mirrors 教程.txt.
    print("Calculating QC metrics...", flush=True)
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

    print("Selecting highly variable genes with seurat_v3...", flush=True)
    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat_v3",
        n_top_genes=args.n_top_genes,
        subset=True,
    )

    trainer_kwargs = {
        "loss_mode": "nb",
        "alpha_recon_lec": 0.5,
        "alpha_recon_lode": 0.5,
    }
    if args.nepoch is not None:
        trainer_kwargs["nepoch"] = args.nepoch
    if args.percent is not None:
        trainer_kwargs["percent"] = args.percent
    if args.force_cpu:
        trainer_kwargs["use_gpu"] = False

    print("Training scTour model...", flush=True)
    tnode = sct.train.Trainer(adata, **trainer_kwargs)
    tnode.train()

    print("Inferring pseudotime, latent space, and vector field...", flush=True)
    adata.obs["ptime"] = tnode.get_time()
    mix_zs, zs, pred_zs = tnode.get_latentsp(alpha_z=0.5, alpha_predz=0.5)
    adata.obsm["X_TNODE"] = mix_zs
    adata.obsm["X_VF"] = tnode.get_vector_field(
        adata.obs["ptime"].values,
        adata.obsm["X_TNODE"],
    )

    print("Computing UMAP from X_TNODE...", flush=True)
    adata = adata[adata.obs["ptime"].values.argsort(), :].copy()
    sc.pp.neighbors(adata, use_rep="X_TNODE", n_neighbors=15)
    sc.tl.umap(adata, min_dist=0.1)

    print("Saving official tutorial-style figure...", flush=True)
    fig, axs = plt.subplots(ncols=2, nrows=2, figsize=(10, 10))
    sc.pl.umap(
        adata,
        color="celltype",
        ax=axs[0, 0],
        legend_loc="on data",
        show=False,
        frameon=False,
    )
    sc.pl.umap(
        adata,
        color="Sample batch",
        ax=axs[0, 1],
        show=False,
        frameon=False,
    )
    sc.pl.umap(
        adata,
        color="ptime",
        ax=axs[1, 0],
        show=False,
        frameon=False,
    )
    sct.vf.plot_vector_field(
        adata,
        zs_key="X_TNODE",
        vf_key="X_VF",
        use_rep_neigh="X_TNODE",
        color="celltype",
        show=False,
        ax=axs[1, 1],
        legend_loc="none",
        frameon=False,
        size=100,
        alpha=0.2,
    )
    fig.tight_layout()
    fig.savefig(output_dir / "official_tutorial_exact_2x2.png", dpi=300)
    plt.close(fig)

    print("Saving result files...", flush=True)
    adata.obs[["celltype", "Sample batch", "Name", "Order", "ptime"]].to_csv(
        output_dir / "official_tutorial_exact_pseudotime.csv",
        encoding="utf-8-sig",
    )
    if getattr(tnode, "log", None):
        import pandas as pd

        pd.DataFrame(tnode.log).to_csv(
            output_dir / "official_tutorial_exact_training_log.csv",
            index=False,
            encoding="utf-8-sig",
        )
    adata.write_h5ad(output_dir / "official_tutorial_exact_result.h5ad")

    with (output_dir / "official_tutorial_exact_parameters.json").open(
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(
            {
                "input": str(input_path),
                "n_top_genes": args.n_top_genes,
                "nepoch": args.nepoch,
                "percent": args.percent,
                "force_cpu": args.force_cpu,
                "note": "This script intentionally follows 教程.txt and does not apply pseudotime reversal or t_key in plot_vector_field.",
            },
            handle,
            indent=2,
            ensure_ascii=False,
        )

    print("Exact official tutorial workflow completed.", flush=True)
    print(f"Output directory: {output_dir}", flush=True)


if __name__ == "__main__":
    main()
