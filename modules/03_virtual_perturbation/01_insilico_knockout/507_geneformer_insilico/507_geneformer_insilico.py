# -*- coding: utf-8 -*-
# =============================================================================
# 507 · Geneformer 零样本嵌入 + in-silico 基因敲除 (foundation-model perturbation)
# -----------------------------------------------------------------------------
# 用基础模型 Geneformer(Theodoris 2023 Nature)做【零样本细胞嵌入】与【in-silico
# 基因删除】(删一个基因→量化细胞状态嵌入的偏移→排序候选调控子)。
#
# ★本机无 GPU / 未下载预训练权重 → Geneformer 主路径标记 "needs GPU"、默认不跑;
#   但内置一条【可本地跑通的诚实基线】(HVG+PCA 嵌入 + 朴素 zero-&-reproject 扰动),
#   生成展示图。两者的对照正是本模块的科学价值(见下"两把屠刀")。
#
# ★两把屠刀(投稿必防):
#   ① Genome Biol 2025 (10.1186/s13059-025-03574-x):零样本 FM 嵌入未必胜过 PCA/scVI
#      → 必须把 Geneformer 嵌入与 PCA 基线(本模块 fig1)同任务对照,别只报 FM。
#   ② Nat Methods 2025 (10.1038/s41592-025-02772-6):DL 扰动预测常输给简单线性基线
#      → Geneformer 的 in-silico 删除排序必须证明优于本模块的朴素基线(fig2)再下结论。
#
# Turnkey(基线): python 507_geneformer_insilico.py        # 本地出图
# Geneformer 全量: python 507_geneformer_insilico.py --run-geneformer  # 需 GPU+权重
# 复用 _framework/pubstyle.py;图中文字英文,注释中文。
# =============================================================================
from __future__ import annotations
import argparse
import sys
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")

# ---- 载入顶刊绘图框架(向上搜 _framework) ----------------------------------
HERE = Path(__file__).resolve().parent
for up in [HERE, *HERE.parents]:
    if (up / "_framework" / "pubstyle.py").exists():
        sys.path.insert(0, str(up / "_framework"))
        break
try:
    from pubstyle import set_pub_style, save_fig, pal, panel_labels, NATURE_W2
except Exception:
    def set_pub_style(*a, **k): pass
    def save_fig(fig, f, dpi=300):
        Path(f).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(str(f) + ".pdf"); fig.savefig(str(f) + ".png", dpi=dpi)
    def pal(n=None, name="npg"):
        import matplotlib; return list(matplotlib.cm.tab10.colors)
    def panel_labels(*a, **k): pass
    NATURE_W2 = 7.0

SEED = 42
np.random.seed(SEED)
N_HUB = 5            # 真驱动基因(hub)数
HUB_GENES = [f"G{i:03d}" for i in range(N_HUB)]


# ---- 合成数据:潜在"激活"程序 + 少数 hub 驱动 + 下游程序 + 背景 -----------
def make_synthetic(n_genes=200, n_cells=800):
    """一个连续"激活"潜变量驱动:N_HUB 个 hub 基因(与激活强相关)+ 30 个下游程序
    基因(由激活驱动)+ 背景。NB 过离散 + dropout。设计意图:hub 是真正的调控子,
    in-silico 删除它应造成最大状态偏移 —— 用于检验扰动方法能否把 hub 排到前面。"""
    rng = np.random.default_rng(SEED)
    act = rng.beta(1.5, 1.5, size=n_cells)                 # 每细胞激活程度 0..1
    log_mu = rng.normal(0.2, 0.25, size=(n_cells, n_genes))
    # hub:强随激活
    for j in range(N_HUB):
        log_mu[:, j] += 2.4 * act + rng.normal(0, 0.1, n_cells)
    # 下游程序 30 基因(由激活驱动,但本身不是"可删的因")
    for j in range(N_HUB, N_HUB + 30):
        log_mu[:, j] += 1.8 * act + rng.normal(0, 0.15, n_cells)
    mu = np.exp(log_mu)
    lam = rng.gamma(shape=2.0, scale=mu / 2.0)             # 负二项过离散
    counts = rng.poisson(lam)
    counts[rng.random(counts.shape) < 0.2] = 0            # 20% dropout
    genes = [f"G{i:03d}" for i in range(n_genes)]
    state = np.where(act < 0.33, "low", np.where(act < 0.66, "mid", "high"))
    return counts.astype(np.float32), genes, act, state


# ---- 诚实基线 ① 嵌入:HVG + PCA + UMAP -------------------------------------
def baseline_embedding(counts, genes):
    import scanpy as sc, anndata as ad
    adata = ad.AnnData(counts.copy()); adata.var_names = genes
    sc.pp.normalize_total(adata, target_sum=1e4); sc.pp.log1p(adata)
    adata.raw = adata
    sc.pp.highly_variable_genes(adata, n_top_genes=80)
    adata = adata[:, adata.var.highly_variable].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=20, random_state=SEED)
    sc.pp.neighbors(adata, use_rep="X_pca", random_state=SEED)
    sc.tl.umap(adata, random_state=SEED)
    return adata


# ---- 诚实基线 ② 朴素 in-silico 扰动:zero-&-reproject ----------------------
def baseline_perturb(counts, genes, candidates):
    """把候选基因在所有细胞中置零,投影回【原始 PCA 载荷】,量化平均位移。
    这是最朴素的扰动基线:只能捕捉基因自身方差的"直接"贡献,看不到下游调控传播
    —— 这一缺口正是 GRN(CellOracle 模块 069)/ FM(Geneformer)宣称的价值所在,
    但必须先证明它们的排序优于本基线。"""
    import numpy as np
    from sklearn.decomposition import PCA
    Xl = np.log1p(counts / counts.sum(1, keepdims=True).clip(min=1) * 1e4)
    mu, sd = Xl.mean(0), Xl.std(0).clip(min=1e-8)
    Z = (Xl - mu) / sd
    pca = PCA(n_components=20, random_state=SEED).fit(Z)
    emb0 = pca.transform(Z)
    gidx = {g: i for i, g in enumerate(genes)}
    res = []
    for g in candidates:
        Xp = counts.copy(); Xp[:, gidx[g]] = 0
        Xlp = np.log1p(Xp / Xp.sum(1, keepdims=True).clip(min=1) * 1e4)
        Zp = (Xlp - mu) / sd
        shift = np.linalg.norm(pca.transform(Zp) - emb0, axis=1).mean()
        res.append((g, float(shift)))
    res.sort(key=lambda x: -x[1])
    return res


# ---- Geneformer 全量路径(需 GPU + 预训练权重;本机默认不跑)----------------
def run_geneformer(counts, genes, candidates, outdir):
    """代表性实现:tokenize → EmbExtractor 零样本嵌入 → InSilicoPerturber 删除扰动。
    需 GPU、`geneformer` 包与 Hugging Face 预训练权重(ctheodoris/Geneformer,数 GB)。
    本机无 GPU/权重时给出清晰提示并退回基线;在 AutoDL RTX5090 等环境可启用。"""
    try:
        import torch
        from geneformer import TranscriptomeTokenizer, EmbExtractor, InSilicoPerturber, InSilicoPerturberStats
    except Exception as e:
        print(f"[geneformer] 依赖缺失({e});需 `pip install geneformer` + 预训练权重。已退回基线。")
        return None
    if not torch.cuda.is_available():
        print("[geneformer] 未检测到 GPU。Geneformer 嵌入/扰动需 GPU(本机无)→ 标记 needs-GPU,退回基线。")
        return None
    # —— 以下为代表性流程骨架(GPU 环境填好路径即可运行):
    # 1) 写出带 ensembl_id / n_counts 的 .h5ad/.loom,TranscriptomeTokenizer().tokenize_data(...)
    # 2) EmbExtractor(model_type="CellClassifier"/"Pretrained", emb_layer=-1).extract_embs(model_dir, token_dir, out, prefix)
    # 3) InSilicoPerturber(perturb_type="delete", genes_to_perturb=candidates_ensembl, ...).perturb_data(model_dir, token_dir, out, prefix)
    # 4) InSilicoPerturberStats(...).get_stats(...) → 每基因 cosine_shift,排序对照基线
    print("[geneformer] GPU 与权重就绪:请按 README 填入 model_dir / token 路径后启用四步流程。")
    return "ready"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-geneformer", action="store_true", help="启用 Geneformer 全量路径(需 GPU+权重)")
    ap.add_argument("--outdir", default=str(HERE / "results"))
    args = ap.parse_args()

    set_pub_style(base_size=11)
    DRES = Path(args.outdir); DAST = HERE / "assets"; DDAT = HERE / "example_data"
    for d in (DRES, DAST, DDAT):
        d.mkdir(parents=True, exist_ok=True)

    # ---- 数据 ---------------------------------------------------------------
    counts, genes, act, state = make_synthetic()
    np.savez(DDAT / "synthetic_counts.npz", X=counts, genes=genes, act=act, state=state)
    print(f"[gen] synthetic: {counts.shape[0]} cells x {counts.shape[1]} genes; "
          f"{N_HUB} hub drivers {HUB_GENES} (for demo only)")

    # ---- 基线 ① 嵌入 -------------------------------------------------------
    adata = baseline_embedding(counts, genes)

    # ---- 基线 ② 朴素 in-silico 扰动(hub + 随机非 hub 对照)----------------
    rng = np.random.default_rng(SEED)
    nonhub = list(rng.choice([g for g in genes if g not in HUB_GENES], 10, replace=False))
    candidates = HUB_GENES + nonhub
    pert = baseline_perturb(counts, genes, candidates)

    import pandas as pd
    dfp = pd.DataFrame(pert, columns=["gene", "shift"])
    dfp["is_hub"] = dfp["gene"].isin(HUB_GENES)
    dfp.to_csv(DRES / "baseline_insilico_ranking.csv", index=False)
    # hub 平均排名 vs 非 hub(基线是否把真驱动排前)
    ranks = {g: i for i, (g, _) in enumerate(pert)}
    hub_rank = np.mean([ranks[g] for g in HUB_GENES])
    print(f"[baseline] hub 平均排名 {hub_rank:.1f}/{len(candidates)}(越小越靠前=基线越能识别真驱动)")

    # ---- (可选)Geneformer 全量路径 --------------------------------------
    if args.run_geneformer:
        run_geneformer(counts, genes, candidates, DRES)
    else:
        print("[info] 未加 --run-geneformer;仅出诚实基线图。Geneformer 嵌入/扰动需 GPU+权重(见 README)。")

    # ---- 出图(Nature 主题)------------------------------------------------
    import matplotlib.pyplot as plt
    # Fig1: 基线 PCA-UMAP 嵌入(着色激活状态)—— Geneformer 零样本嵌入须与此对照
    fig, a = plt.subplots(figsize=(4.4, 3.6))
    U = adata.obsm["X_umap"]; states = ["low", "mid", "high"]
    cmap = {s: pal(3, "okabe_ito")[i] for i, s in enumerate(states)}
    for s in states:
        m = state == s
        a.scatter(U[m, 0], U[m, 1], s=10, c=cmap[s], label=s, alpha=0.8, linewidths=0)
    a.set_title("Baseline embedding (HVG + PCA)\ncompare Geneformer zero-shot against this")
    a.set_xlabel("UMAP1"); a.set_ylabel("UMAP2"); a.set_xticks([]); a.set_yticks([])
    a.legend(title="Activation", markerscale=1.5)
    fig.tight_layout(); save_fig(fig, DAST / "baseline_embedding"); plt.close(fig)

    # Fig2: 朴素 in-silico 扰动排序 lollipop(hub 高亮)—— Geneformer 删除须证明优于此。
    #       棒棒糖图替代条形图(顶刊更常用、低墨水、排序可读性更好)。
    fig, a = plt.subplots(figsize=(4.7, 4.4))
    order = dfp.sort_values("shift").reset_index(drop=True)
    y = np.arange(len(order))
    cols = [pal(2, "npg")[0] if h else "#9AA0A6" for h in order["is_hub"]]
    a.hlines(y, 0, order["shift"], color=cols, lw=1.6, alpha=0.65)
    a.scatter(order["shift"], y, c=cols, s=52, zorder=3, edgecolor="black", linewidth=0.5)
    a.set_yticks(y); a.set_yticklabels(order["gene"], fontsize=8)
    a.set_ylim(-0.6, len(order) - 0.4); a.set_xlim(0, order["shift"].max() * 1.10)
    a.set_xlabel("Mean embedding shift after in-silico KO")
    a.set_title("Naive baseline in-silico perturbation\n(red = true hub driver)")
    from matplotlib.lines import Line2D
    a.legend(handles=[Line2D([0], [0], marker="o", color="w", markerfacecolor=pal(2, "npg")[0],
                             markersize=8, label="hub driver"),
                      Line2D([0], [0], marker="o", color="w", markerfacecolor="#9AA0A6",
                             markersize=8, label="non-hub")], fontsize=8, loc="lower right")
    fig.tight_layout(); save_fig(fig, DAST / "baseline_insilico_ranking"); plt.close(fig)

    print("[fig] assets/: baseline_embedding, baseline_insilico_ranking (.pdf+.png)")
    # 依赖快照(铁律6)
    try:
        import session_info
        session_info.show(write_req_file=False)
    except Exception:
        import scanpy as sc
        print(f"[env] numpy={np.__version__}, scanpy={sc.__version__}; "
              f"Geneformer 路径见 README(needs GPU + pretrained weights)")


if __name__ == "__main__":
    main()
