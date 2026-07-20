# -*- coding: utf-8 -*-
# =============================================================================
# 581 · veloAgent —— 空间信息驱动的 RNA velocity 与 in-silico 扰动
# -----------------------------------------------------------------------------
# 上游  : veloAgent (Raghavan/Yoon/Fonseca/Li/Ding, Mol Syst Biol 2026)
#         repo  https://github.com/mcgilldinglab/veloAgent
#         DOI   10.1038/s44320-026-00213-w   PMID 42092184
# 澄清  : 这里的 "Agent" 指 **agent-based model (ABM,基于智能体的空间模拟)**,
#         不是 LLM agent。上游 **不需要任何外部 LLM API / API key**(已读源码与
#         README 确认:依赖是 mesa 的 ABM + PyTorch VAE + STRING 基因-基因先验)。
# 本模块 : ① 一个**本机零依赖新增就能跑通的朴素基线**(scVelo 速度 + 空间 kNN 平滑
#            + 逐基因敲除的余弦扰动打分),作为 veloAgent 的对照下限;
#          ② 一个**守卫式的 veloAgent 封装**:装了才跑,没装就打印真实安装命令与
#            已核实的官方调用顺序,绝不静默降级、绝不假装跑了 veloAgent。
# 出图  : 空间速度场 quiver / 每细胞一致性 raincloud / 扰动打分 lollipop(无条形图)
# =============================================================================
"""581 · veloAgent — spatially-informed RNA velocity + in-silico perturbation.

Baseline (always runnable): scVelo velocity -> spatial kNN smoothing -> per-gene
knockout scored by cosine shift of the velocity field.
veloAgent path (--run-veloagent): guarded; needs `pip install .` from the repo.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
# 框架统一出图样式(pubstyle.py 在 modules/_framework/)
sys.path.insert(0, str(HERE.parents[2] / "_framework"))
from pubstyle import (set_pub_style, pal, save_fig, CMAP_CONT,  # noqa: E402
                      NATURE_W1, NATURE_W2)

import matplotlib                                              # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                # noqa: E402

SEED = 2026


# --------------------------------------------------------------------------- #
# 0. 读入
# --------------------------------------------------------------------------- #
def load_data(datadir: Path):
    """读 example_data/ 的四个 csv(# 开头为注释行),组装成 AnnData。"""
    import anndata as ad

    S = pd.read_csv(datadir / "spliced.csv", comment="#", index_col=0)
    U = pd.read_csv(datadir / "unspliced.csv", comment="#", index_col=0)
    meta = pd.read_csv(datadir / "spatial_meta.csv", comment="#").set_index("cell")
    meta = meta.loc[S.index]

    A = ad.AnnData(S.values.astype("float32"))
    A.obs_names, A.var_names = S.index.astype(str), S.columns.astype(str)
    A.layers["spliced"] = S.values.astype("float32")
    A.layers["unspliced"] = U.loc[S.index, S.columns].values.astype("float32")
    A.obsm["spatial"] = meta[["x", "y"]].values.astype("float32")
    A.obs["true_time"] = meta["true_time"].values
    A.obs["cluster"] = pd.Categorical(meta["cluster"].values)

    truth_f = datadir / "gene_truth.csv"
    if truth_f.exists():                       # 合成数据自带的地面真值(仅用于评估)
        tr = pd.read_csv(truth_f, comment="#").set_index("gene")
        A.var["true_direction"] = tr.reindex(A.var_names)["true_direction"].values
    return A


# --------------------------------------------------------------------------- #
# 1. 基线 A:非空间 scVelo 速度场
# --------------------------------------------------------------------------- #
def velocity_baseline(adata):
    """朴素下限:标准 scVelo,完全不看空间坐标。"""
    import scanpy as sc
    import scvelo as scv

    scv.pp.filter_and_normalize(adata, min_shared_counts=0)
    scv.pp.moments(adata, n_pcs=20, n_neighbors=15)
    # stochastic 在近退化数据上会解不出二阶矩,deterministic 是安全兜底
    mode_used = None
    for mode in ("stochastic", "deterministic"):
        try:
            scv.tl.velocity(adata, mode=mode)
            mode_used = mode
            break
        except Exception as e:
            print(f"       velocity mode='{mode}' failed ({type(e).__name__}), fallback")
    if mode_used is None:
        sys.exit("[581] scVelo 无法在该输入上拟合速度场")
    sc.pp.neighbors(adata, n_neighbors=15)
    sc.tl.umap(adata)
    return mode_used


# --------------------------------------------------------------------------- #
# 2. 基线 B:空间 kNN 平滑(ABM 的线性朴素替身)
# --------------------------------------------------------------------------- #
def spatial_smooth(adata, k: int = 12, tau: float = 2.0):
    """把每个细胞的速度向量替换为其空间 k 近邻的距离加权平均。

    这是 veloAgent 的 agent-based 微环境模拟(mesa CellModel)所对应的**最朴素
    线性对照**:同样只用「邻居的速度」这一条信息,但不做任何 agent 规则/迭代。
    上游若声称 ABM 更好,至少要赢过这条线。
    """
    from sklearn.neighbors import NearestNeighbors

    V = np.nan_to_num(np.asarray(adata.layers["velocity"], dtype=float))
    xy = np.asarray(adata.obsm["spatial"], dtype=float)
    nn = NearestNeighbors(n_neighbors=min(k + 1, len(xy))).fit(xy)
    dist, idx = nn.kneighbors(xy)
    w = np.exp(-dist / max(tau, 1e-8))
    w /= w.sum(axis=1, keepdims=True)
    adata.layers["velocity_spatial"] = np.einsum("ck,ckg->cg", w, V[idx]).astype("float32")
    return {"k": k, "tau": tau}


# --------------------------------------------------------------------------- #
# 3. 评估:每细胞速度方向 vs 地面真值
# --------------------------------------------------------------------------- #
def _row_cosine(A, B):
    an = np.linalg.norm(A, axis=1, keepdims=True)
    bn = np.linalg.norm(B, axis=1, keepdims=True)
    A = A / np.where(an == 0, 1, an)
    B = B / np.where(bn == 0, 1, bn)
    return np.sum(A * B, axis=1)


def consistency(adata, layer: str):
    """与合成真值方向的逐细胞余弦相似度;没有真值就返回 NaN 数组。"""
    if "true_direction" not in adata.var:
        return np.full(adata.n_obs, np.nan)
    truth = np.tile(np.asarray(adata.var["true_direction"], dtype=float),
                    (adata.n_obs, 1))
    V = np.nan_to_num(np.asarray(adata.layers[layer], dtype=float))
    return _row_cosine(V, truth)


# --------------------------------------------------------------------------- #
# 4. 基线 C:逐基因 in-silico 敲除(余弦位移打分)
# --------------------------------------------------------------------------- #
def knockout_scores(adata, layer: str = "velocity_spatial"):
    """把某基因的速度分量置零,量化整个速度场被推动了多少。

    打分口径对齐上游 `veloagent.perturbation_score(..., metric_option=2)` 的
    「扰动前后速度场余弦相似度」思路(已从源码读到);但**实现是本模块自己的
    numpy 朴素版**,不调用上游、也不假装等价于上游的 alpha/beta 动力学置零。
    """
    V = np.nan_to_num(np.asarray(adata.layers[layer], dtype=float))
    rows = []
    for j, g in enumerate(adata.var_names):
        Vp = V.copy()
        Vp[:, j] = 0.0
        cos = _row_cosine(Vp, V)
        rows.append({"gene": g, "score": float(1.0 - np.nanmean(cos)),
                     "mean_abs_velocity": float(np.nanmean(np.abs(V[:, j])))})
    df = pd.DataFrame(rows).sort_values("score", ascending=False).reset_index(drop=True)
    df["rank"] = np.arange(1, len(df) + 1)
    return df


# --------------------------------------------------------------------------- #
# 5. veloAgent 守卫式封装
# --------------------------------------------------------------------------- #
# 每一行后面的 "# <-" 是**上游源码里该符号的定义位置**(2026-07-21 逐个 grep 本地克隆
# C:\Users\fsy\Desktop\upstream-sources\581_veloAgent 核对,参数名与默认值均已对齐)。
VELOAGENT_WORKFLOW = [
    "veloagent.preprocess(adata)"
    "  # <- src/veloagent/preprocessing.py:8  preprocess(data, num_genes=2000, min_count=20, log_norm=True)",

    "vae, optimizer, loss_fn = veloagent.get_vae(adata, z_dim=5, lr=1e-2)"
    "  # <- src/veloagent/vae.py:318  get_vae(adata, z_dim, lr=1e-2) -> (vae, optimizer, loss_fn)",

    "best = veloagent.train_vae(adata, vae, optimizer, loss_fn, "
    "adata.uns['neighbors']['indices'], patience=25, num_epochs=1000, batch=256, device=device)"
    "  # <- src/veloagent/vae.py:346  train_vae(adata, vae, optimizer, loss_fn, global_nb_indices, "
    "patience, num_epochs=1000, batch=256, device='cpu', verbose=True) -> best_vae (state_dict)",

    "vae.load_state_dict(best); veloagent.get_embedding(adata, vae, device)"
    "  # <- src/veloagent/vae.py:472  get_embedding(adata, vae_model, device); 就地写 obsm['cell_embed']",

    "paths = veloagent.load_protein_paths(species='mouse', base='data/conn_mat')"
    "  # <- src/veloagent/gg_net.py:14  返回 [info, aliases, links] 三条路径;仅 mouse/human/chicken",

    "conn_mat = veloagent.create_con_mat(adata, adata.n_vars, paths[0], paths[1], paths[2], "
    "varname='index', confidence=False, conf_threshold=400)"
    "  # <- src/veloagent/gg_net.py:53  create_con_mat(data, num_genes, prot_names, prot_alias, "
    "gene_conn, varname, confidence=False, conf_threshold=400)  [paths 顺序与形参一一对应]",

    "genenet = veloagent.GeneNet(in_dim=adata.obsm['cell_embed'][0].shape[0], "
    "gene_dim=adata.n_vars, conn=conn_mat)"
    "  # <- src/veloagent/gg_net.py:334 (class) / :383 __init__(self, in_dim, gene_dim, conn=None)",

    "veloagent.train_gg(num_epochs=500, data=adata, embed_basis='cell_embed', genenet=genenet, "
    "optimizer=optimizer, patience=25, num_nbrs=30, dt=0.5, batch=256, device=device)"
    "  # <- src/veloagent/gg_net.py:629  (源码默认 patience=10, dt=0.3;这里的 25/0.5 是 tutorial 显式值)",

    "veloagent.train_nbr(num_epochs=50, data=adata, embed_basis='cell_embed', genenet=genenet, "
    "optimizer=optimizer, num_nbrs=30, dt=0.5, batch=256, device=device, tau_nbr=0.2)"
    "  # <- src/veloagent/gg_net.py:760  train_nbr(..., num_nbrs=30, dt=0.3, batch=256, "
    "device='cpu', tau_nbr=0.2)",

    "abm = veloagent.CellModel(adata, steps=100, tau=2, nbr_radius=30, sig_ratio=0.7); abm.step()"
    "  # <- src/veloagent/abm.py:157 (class) / :255 __init__(self, adata, steps, tau=2, "
    "nbr_radius=40, sig_ratio=0.7, max_gene_norm=10.0)  [★源码默认 nbr_radius=40,tutorial 传 30]",

    "scores = veloagent.perturbation_score(adata, cluster_name='clusters', cluster_edges=[...], "
    "vel_key='velocity_u', metric_option=2, pert_param='alpha', dt=0.5)"
    "  # <- src/veloagent/perturbations.py:10  perturbation_score(data, cluster_name, cluster_edges, "
    "vel_key='velocity_u', metric_option=1, pert_param='alpha', dt=0.5) -> DataFrame(index=var_names, col='score')\n"
    "  #    ★ cluster_edges 的形状随 metric_option 变:option 1 = [(src,dst), ...] 成对元组;"
    "option 2 = ['c1','c2', ...] 扁平列表(见 tutorial/perturbation_tutorial.ipynb)",
]

INSTALL_HINT = (
    "conda create -n veloagent python=3.11 && conda activate veloagent && "
    "git clone https://github.com/mcgilldinglab/veloAgent.git && cd veloAgent && pip install .\n"
    "       # ★ 上游 README 写的 `conda env create -f environment.yml` 已失效:该文件被上游"
    " commit 63df2a4 'remove environment file' 删除,当前仓库里没有 environment.yml。\n"
    "       # 依赖以 pyproject.toml 为准(requires-python >= 3.11.8);PyTorch 未列在 deps 里,"
    "需按平台单独装:https://pytorch.org/get-started/locally/\n"
    "       # ★ pyproject 把 anndata==0.9.2 / numpy==1.24.4 / scanpy==1.9.8 / scvelo==0.3.2 "
    "全部钉死,与本机环境冲突 —— 必须装在独立 conda env 里,别污染本机。\n"
    "       # 另需手动下载 STRING DB 三件套 (protein.info/aliases/links) 放到 "
    "<base>/<species>/,再 veloagent.load_protein_paths(species, base)\n"
    "       # 注意:perturbations.py 是 `from veloproj import *`(由 pyproject 的 "
    "`veloae @ git+https://github.com/qiaochen/VeloAE.git` 提供),abm.py 依赖 mesa==2.1.5"
)


def run_veloagent(adata):
    """装了 veloagent 才走这条路;没装就如实报告并给出真实安装命令。

    这里**不替上游拟合任何模型**——本模块从未在本机跑通 veloAgent 训练。
    上面 VELOAGENT_WORKFLOW 的每一行都抄自官方 tutorial/veloagent_tutorial.ipynb
    与 src/veloagent/*.py 的真实签名(见 README「已核实来源」),不是臆造。
    """
    try:
        import veloagent  # noqa: F401
    except ImportError as e:
        return {"status": "skipped",
                "reason": f"veloagent 未安装 ({e.name})",
                "install": INSTALL_HINT,
                "workflow_from_official_tutorial": VELOAGENT_WORKFLOW}
    missing = []
    for dep in ("mesa", "torch", "veloproj"):
        try:
            __import__(dep)
        except ImportError:
            missing.append(dep)
    exported = [n for n in ("preprocess", "get_vae", "train_vae", "get_embedding",
                            "load_protein_paths", "create_con_mat", "GeneNet",
                            "train_gg", "train_nbr", "CellModel",
                            "perturbation_score", "perturb", "perturb_score_plt")
                if hasattr(veloagent, n)]
    return {"status": "ready" if not missing else "partial",
            "missing_deps": missing,
            "exported_api_present": exported,
            "next": "按 VELOAGENT_WORKFLOW 执行;STRING DB 文件需自行下载",
            "workflow_from_official_tutorial": VELOAGENT_WORKFLOW}


# --------------------------------------------------------------------------- #
# 6. 出图(全部非条形图)
# --------------------------------------------------------------------------- #
def fig_spatial_quiver(adata, cons_raw, cons_sp, outstem):
    """空间速度场 quiver:左=非空间 scVelo,右=空间平滑后。用 UMAP1 作向量投影载体。"""
    xy = np.asarray(adata.obsm["spatial"], dtype=float)
    # 把高维速度投影到 2D 空间坐标:用 PCA 载荷把 gene-space 速度映射到 xy 方向
    from sklearn.decomposition import PCA
    X = np.asarray(adata.layers["Ms"] if "Ms" in adata.layers else adata.X, dtype=float)
    p = PCA(n_components=2, random_state=SEED).fit(X)
    # 用最小二乘把 PCA 空间对齐到物理 xy,使箭头方向可读
    Z = p.transform(X)
    Z1 = np.c_[Z, np.ones(len(Z))]
    beta, *_ = np.linalg.lstsq(Z1, xy, rcond=None)
    A = beta[:2]                                    # 2x2 线性映射 PCA->xy

    finite = np.concatenate([cons_raw[~np.isnan(cons_raw)], cons_sp[~np.isnan(cons_sp)]])
    vmin, vmax = (float(np.min(finite)), float(np.max(finite))) if finite.size else (-1.0, 1.0)

    fig, axes = plt.subplots(1, 2, figsize=(NATURE_W2, 3.4), constrained_layout=True)
    for ax, layer, cons, ttl in (
            (axes[0], "velocity", cons_raw, "scVelo (non-spatial)"),
            (axes[1], "velocity_spatial", cons_sp, "+ spatial kNN smoothing")):
        V = np.nan_to_num(np.asarray(adata.layers[layer], dtype=float))
        d = (p.transform(X + V) - Z) @ A
        n = np.linalg.norm(d, axis=1, keepdims=True)
        d = d / np.where(n == 0, 1, n)
        sc_ = ax.scatter(xy[:, 0], xy[:, 1], c=cons, cmap=CMAP_CONT, s=11,
                         vmin=vmin, vmax=vmax, linewidths=0, zorder=1)
        ax.quiver(xy[:, 0], xy[:, 1], d[:, 0], d[:, 1], color="#2B2B2B",
                  scale=52, width=0.0022, headwidth=3.5, headlength=4.5,
                  alpha=0.8, zorder=2)
        ax.set_title(ttl)
        ax.set_xlabel("spatial x"); ax.set_ylabel("spatial y")
        ax.set_aspect("equal")
    cb = fig.colorbar(sc_, ax=axes, shrink=0.85, pad=0.02)
    cb.set_label("cosine vs ground-truth direction")
    save_fig(fig, outstem); plt.close(fig)


def fig_consistency_raincloud(cons_raw, cons_sp, outstem):
    """每细胞一致性 raincloud(violin + box + jitter),对照两条基线。"""
    data = [cons_raw[~np.isnan(cons_raw)], cons_sp[~np.isnan(cons_sp)]]
    labels = ["scVelo\n(non-spatial)", "scVelo + spatial\nkNN smoothing"]
    cols = pal(2, "npg")
    rng = np.random.default_rng(SEED)

    fig, ax = plt.subplots(figsize=(NATURE_W1 + 0.8, 3.4), constrained_layout=True)
    vp = ax.violinplot(data, positions=[0, 1], widths=0.7, showextrema=False)
    for b, c in zip(vp["bodies"], cols):
        b.set_facecolor(c); b.set_alpha(0.35); b.set_edgecolor(c); b.set_linewidth(1.0)
    bp = ax.boxplot(data, positions=[0, 1], widths=0.12, showfliers=False,
                    patch_artist=True)
    for patch in bp["boxes"]:
        patch.set_facecolor("white"); patch.set_edgecolor("black")
    for i, (d, c) in enumerate(zip(data, cols)):
        ax.scatter(i + 0.26 + rng.normal(0, 0.035, len(d)), d, s=5, color=c,
                   alpha=0.45, linewidths=0)
    ax.axhline(0, color="#999999", lw=0.8, ls="--")
    ax.set_xticks([0, 1]); ax.set_xticklabels(labels)
    ax.set_ylabel("per-cell cosine vs ground truth")
    ax.set_title("Velocity direction accuracy")
    save_fig(fig, outstem); plt.close(fig)


def fig_knockout_lollipop(df, outstem, top: int = 18):
    """in-silico 敲除打分 lollipop(明确不用条形图)。"""
    d = df.head(top).iloc[::-1]
    y = np.arange(len(d))
    fig, ax = plt.subplots(figsize=(NATURE_W1 + 0.6, 0.19 * top + 1.3),
                           constrained_layout=True)
    ax.hlines(y, 0, d["score"], color="#C9C9C9", lw=1.1, zorder=1)
    sc_ = ax.scatter(d["score"], y, c=d["mean_abs_velocity"], cmap=CMAP_CONT,
                     s=42, zorder=2, edgecolors="black", linewidths=0.4)
    ax.set_yticks(y); ax.set_yticklabels(d["gene"])
    ax.set_xlabel("perturbation score  (1 - mean cosine after knockout)")
    ax.set_title("Baseline in-silico knockout ranking")
    cb = fig.colorbar(sc_, ax=ax, shrink=0.8, pad=0.02)
    cb.set_label("mean |velocity|")
    save_fig(fig, outstem); plt.close(fig)


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description="581 veloAgent baseline + guarded wrapper")
    ap.add_argument("--datadir", default=str(HERE / "example_data"),
                    help="含 spliced.csv/unspliced.csv/spatial_meta.csv 的目录")
    ap.add_argument("--outdir", default=str(HERE / "results"))
    ap.add_argument("--assets", default=str(HERE / "assets"))
    ap.add_argument("--knn", type=int, default=12, help="空间平滑近邻数")
    ap.add_argument("--tau", type=float, default=2.0, help="空间平滑距离衰减尺度")
    ap.add_argument("--run-veloagent", action="store_true",
                    help="尝试调用上游 veloAgent(需自行安装,见 README)")
    a = ap.parse_args()

    np.random.seed(SEED)
    set_pub_style()
    outdir, assets = Path(a.outdir), Path(a.assets)
    outdir.mkdir(parents=True, exist_ok=True); assets.mkdir(parents=True, exist_ok=True)

    print("[581] Step 1  读入空间 spliced/unspliced 数据")
    adata = load_data(Path(a.datadir))
    print(f"       cells={adata.n_obs}  genes={adata.n_vars}  spatial={adata.obsm['spatial'].shape}")

    print("[581] Step 2  基线 A:非空间 scVelo 速度场")
    mode = velocity_baseline(adata)
    print(f"       scVelo mode = {mode}")

    print("[581] Step 3  基线 B:空间 kNN 平滑(ABM 的朴素线性对照)")
    sm = spatial_smooth(adata, k=a.knn, tau=a.tau)
    print(f"       k={sm['k']}  tau={sm['tau']}")

    print("[581] Step 4  评估:逐细胞方向一致性")
    cons_raw = consistency(adata, "velocity")
    cons_sp = consistency(adata, "velocity_spatial")
    summary = {
        "n_cells": int(adata.n_obs), "n_genes": int(adata.n_vars),
        "scvelo_mode": mode, "spatial_knn": sm,
        "mean_cosine_nonspatial": round(float(np.nanmean(cons_raw)), 4),
        "mean_cosine_spatial": round(float(np.nanmean(cons_sp)), 4),
    }
    from scipy.stats import wilcoxon, spearmanr
    ok = ~(np.isnan(cons_raw) | np.isnan(cons_sp))
    if ok.sum() > 10 and np.any(cons_raw[ok] != cons_sp[ok]):
        summary["wilcoxon_p_spatial_vs_nonspatial"] = float(
            wilcoxon(cons_raw[ok], cons_sp[ok]).pvalue)
    try:
        import scvelo as scv
        scv.tl.velocity_graph(adata)
        scv.tl.velocity_pseudotime(adata)
        summary["pseudotime_vs_true_time_rho"] = round(float(
            spearmanr(adata.obs["velocity_pseudotime"], adata.obs["true_time"]).correlation), 3)
    except Exception as e:
        summary["pseudotime"] = f"skipped: {type(e).__name__}"
    for k, v in summary.items():
        print(f"       {k}: {v}")

    print("[581] Step 5  基线 C:逐基因 in-silico 敲除打分")
    ko = knockout_scores(adata, "velocity_spatial")
    ko.to_csv(outdir / "581_knockout_scores.csv", index=False)
    print(f"       top gene: {ko.iloc[0]['gene']}  score={ko.iloc[0]['score']:.4f}")

    print("[581] Step 6  出图")
    fig_spatial_quiver(adata, cons_raw, cons_sp, str(assets / "581_spatial_velocity_field"))
    fig_consistency_raincloud(cons_raw, cons_sp, str(assets / "581_consistency_raincloud"))
    fig_knockout_lollipop(ko, str(assets / "581_knockout_lollipop"))
    print(f"       PNG/PDF -> {assets}")

    pd.DataFrame({"cell": adata.obs_names,
                  "cosine_nonspatial": cons_raw,
                  "cosine_spatial": cons_sp,
                  "true_time": adata.obs["true_time"].values}).to_csv(
        outdir / "581_percell_consistency.csv", index=False)

    print("[581] Step 7  veloAgent 路径")
    va = run_veloagent(adata) if a.run_veloagent else {
        "status": "not requested", "hint": "加 --run-veloagent 以尝试调用上游"}
    print(f"       status: {va.get('status')}")
    if va.get("reason"):
        print(f"       reason: {va['reason']}")
        print(f"       install: {va['install']}")

    with open(outdir / "581_summary.json", "w", encoding="utf-8") as fh:
        json.dump({"baseline": summary, "veloagent": va,
                   "top10_knockout": ko.head(10).to_dict("records")},
                  fh, indent=1, ensure_ascii=False, default=str)
    print(f"[581] done -> {outdir / '581_summary.json'}")


if __name__ == "__main__":
    main()
