# -*- coding: utf-8 -*-
# =============================================================================
# 568 · scPRINT —— 单细胞基础模型:细胞嵌入 / 去噪 / 注意力推断基因调控网络
# -----------------------------------------------------------------------------
# 上游论文 : Kalfon J, Samaran J, Peyré G, Cantini L.
#            "scPRINT: pre-training on 50 million cells allows robust gene
#            network predictions." Nat Commun. 2025 Apr 16;16(1):3607.
#            PMID 40240364 · doi:10.1038/s41467-025-58699-1  (已用 NCBI E-utilities 核实)
# 上游仓库 : https://github.com/jkobject/scPRINT  (README 指向新址 cantinilab/scPRINT)
# 本地源码 : C:\Users\fsy\Desktop\upstream-sources\568_scPRINT\
#
# 本模块定位(诚实边界):
#   scPRINT 官方三件套 Embedder / Denoiser / GNInfer 依赖 scdataloader + lamindb +
#   bionty + lightning + HuggingFace 权重(见 pyproject.toml),上游只在
#   MacOS / Ubuntu20.04 + Python3.10 上测过(无 triton 可跑 CPU,但载权重需
#   transformer="normal";实用规模仍应上 GPU)。
#   本机不装包,故官方路径写成 **守卫式封装**(--run-scprint),
#   模块主体是一条 **本机零改动可跑的朴素基线**,覆盖 scPRINT 的三个任务:
#     A) GRN 推断  —— 共表达/偏相关 打分 TF→target,对照真值网络算 AUPRC/AUROC
#     B) 去噪      —— 分子交叉验证(binomial split)+ kNN 池化,对照"不去噪"
#     C) 细胞嵌入  —— PCA 嵌入 + kNN 分类的交叉验证 macro-F1
#   任何"基础模型更好"的说法,必须先打赢这三条地板线
#   (Genome Biol 2025:单细胞 FM zero-shot 常打不过 PCA;
#    Nat Methods 2025:DL 扰动预测常打不过线性基线)。
#
# 依赖   : numpy pandas scipy scikit-learn matplotlib(本机已装)
# 输入   : example_data/ 的 3 个 CSV(synthetic, for demo only)
# 输出   : results/ 指标表 + JSON;assets/ 4 张顶刊风图(无条形图)
# =============================================================================
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.decomposition import PCA
from sklearn.metrics import (average_precision_score, f1_score,
                             mean_squared_error, roc_auc_score, roc_curve,
                             precision_recall_curve)
from sklearn.model_selection import StratifiedKFold
from sklearn.neighbors import KNeighborsClassifier, NearestNeighbors
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parent
FRAMEWORK = ROOT.parents[2] / "_framework"
sys.path.insert(0, str(FRAMEWORK))
from pubstyle import CMAP_CONT, NATURE_W2, pal, save_fig, set_pub_style  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402


# =============================================================================
# 0. 合成示例数据(synthetic, for demo only)
# =============================================================================
def make_synthetic(outdir: Path, seed: int = 42,
                   n_cells: int = 900, n_tf: int = 12, n_target: int = 108) -> None:
    """生成带 **已知真值 GRN** 的合成 scRNA 计数矩阵。

    设计意图(否则 GRN 评估会退化成"随便一个相关系数就满分",没有信息量):
      · TF 的"活性"= 细胞类型均值 + 3 个共享隐因子 + 细胞噪声
        → TF 之间彼此相关,共表达难以分清"哪个 TF"是真上游;
      · 每个 TF 调控 4–8 个 target,target 对数表达率 = 基线 + Σ w·TF活性;
      · TF 的 mRNA 只是活性的**含噪代理**(系数 0.45,噪声 0.8)
        → 复刻真实生物学里"转录本 ≠ 蛋白活性"的根本困难;
      · 未被任何 TF 选中的 target 成为真阴性边(候选边中约 94% 为阴性)。
    """
    rng = np.random.default_rng(seed)
    cts = np.repeat(["Ductal", "Immune", "Stromal"], n_cells // 3)
    ct_codes = pd.Categorical(cts).codes

    tf_names = [f"TF{i:02d}" for i in range(n_tf)]
    tg_names = [f"GENE{i:03d}" for i in range(n_target)]

    # --- TF 活性:细胞类型均值 + 共享隐因子(TF 之间相关)+ 细胞噪声 ---
    # 共享隐因子是关键:真实数据里 TF 彼此高度共表达,共表达打分很难分清
    # "哪一个 TF" 才是某个 target 的真正上游 —— 这正是 GRN 推断的核心困难。
    ct_mean = rng.normal(0.0, 0.90, size=(3, n_tf))
    n_factor = 3
    F = rng.normal(0.0, 1.0, size=(n_cells, n_factor))
    L = rng.normal(0.0, 0.85, size=(n_factor, n_tf))
    tf_act = ct_mean[ct_codes] + F @ L + rng.normal(0.0, 0.45, size=(n_cells, n_tf))

    # --- 真值 GRN:TF × target 稀疏权重 ---
    W = np.zeros((n_tf, n_target))
    for t in range(n_tf):
        n_edge = rng.integers(4, 9)
        tgts = rng.choice(n_target, size=n_edge, replace=False)
        W[t, tgts] = rng.choice([-1.0, 1.0], size=n_edge) * rng.uniform(0.55, 1.15, size=n_edge)

    # --- 表达率 ---
    tg_base = rng.uniform(0.1, 1.5, size=n_target)
    log_rate_tg = tg_base[None, :] + tf_act @ W + rng.normal(0, 0.45, size=(n_cells, n_target))
    # TF 的 mRNA 只是其"活性"的**含噪代理**(真实生物学:翻译后修饰、核转位等
    # 都不在转录本里)。这一条让共表达基线注定打不满,是模块的对照价值所在。
    tf_base = rng.uniform(0.5, 1.5, size=n_tf)
    log_rate_tf = tf_base[None, :] + 0.45 * tf_act + rng.normal(0, 0.8, size=(n_cells, n_tf))

    log_rate = np.hstack([log_rate_tf, log_rate_tg])
    depth = rng.lognormal(mean=0.0, sigma=0.28, size=(n_cells, 1))   # 文库深度差异
    counts = rng.poisson(np.exp(np.clip(log_rate, -6, 6)) * depth)

    genes = tf_names + tg_names
    cell_ids = [f"CELL{i:04d}" for i in range(n_cells)]
    outdir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(counts, index=cell_ids, columns=genes).to_csv(outdir / "counts.csv")
    pd.DataFrame({"cell_id": cell_ids, "cell_type": cts,
                  "organism_ontology_term_id": "NCBITaxon:9606"}).to_csv(
        outdir / "cell_meta.csv", index=False)

    src, dst, wt = [], [], []
    for t in range(n_tf):
        for g in np.nonzero(W[t])[0]:
            src.append(tf_names[t]); dst.append(tg_names[g]); wt.append(round(float(W[t, g]), 4))
    pd.DataFrame({"tf": src, "target": dst, "weight": wt}).to_csv(
        outdir / "true_grn_edges.csv", index=False)
    (outdir / "README.txt").write_text(
        "synthetic, for demo only —— 由 568_scprint_foundation_grn.py --regen-data 生成。\n"
        "counts.csv          : 细胞 × 基因 原始计数(行索引 cell_id)\n"
        "cell_meta.csv       : cell_id, cell_type, organism_ontology_term_id\n"
        "true_grn_edges.csv  : 真值调控边 tf,target,weight(仅合成数据才有)\n",
        encoding="utf-8")


# =============================================================================
# 1. 预处理工具
# =============================================================================
def cpm_log1p(counts: np.ndarray, target_sum: float = 1e4) -> np.ndarray:
    """等价 scanpy.pp.normalize_total(target_sum) + log1p 的最小实现。"""
    tot = counts.sum(axis=1, keepdims=True)
    tot[tot == 0] = 1.0
    return np.log1p(counts / tot * target_sum)


# --- 下面两个函数逐行对照上游源码实现,不是我自创的评估口径 ---
# 源码:scprint/tasks/denoise.py L359-L383  def split_molecules(...)
def split_molecules(umis: np.ndarray, data_split: float,
                    overlap_factor: float = 0.0,
                    random_state: np.random.RandomState | None = None):
    """分子交叉验证:把 UMI 计数二项拆成两份(上游 denoise.py::split_molecules 同式)。"""
    if random_state is None:
        random_state = np.random.RandomState()
    umis_X_disjoint = random_state.binomial(umis, data_split - overlap_factor)
    umis_Y_disjoint = random_state.binomial(
        umis - umis_X_disjoint, (1 - data_split) / (1 - data_split + overlap_factor))
    ov = umis - umis_X_disjoint - umis_Y_disjoint
    return umis_X_disjoint + ov, umis_Y_disjoint + ov


# 源码:scprint/tasks/denoise.py L355-L356  def poisson_nll_loss(...)
def poisson_nll_loss(y_pred: np.ndarray, y_true: np.ndarray) -> float:
    return float((y_pred - y_true * np.log(y_pred + 1e-6)).mean())


def knn_pool_denoise(counts: np.ndarray, k: int, n_pcs: int = 10, seed: int = 0) -> np.ndarray:
    """kNN 池化去噪:PCA→找 k 近邻(含自身)→原始计数求和。

    这是 Wagner & Yan kNN-smoothing 的**简化一步版**。上游把该算法整份 vendored 在
    scprint/tasks/knn_smooth.py(knn_smoothing(X, k, d=10, dither=0.03, seed=0)),
    并在 denoise.py L346-L350 的 withknn(adata, k=10) 里作为对照方法调用。
    本函数只用 numpy+sklearn 重写最小版本,便于本机零依赖跑通。
    """
    logx = cpm_log1p(counts)
    pcs = PCA(n_components=min(n_pcs, logx.shape[1] - 1), random_state=seed).fit_transform(logx)
    nn = NearestNeighbors(n_neighbors=k).fit(pcs)
    _, idx = nn.kneighbors(pcs)
    return counts[idx].sum(axis=1).astype(float)


# =============================================================================
# 2. 任务 A —— GRN 推断基线
# =============================================================================
def score_grn(logx: np.ndarray, genes: list[str], tfs: list[str]) -> dict[str, np.ndarray]:
    """三种朴素 TF→target 打分(全部 |·|,只比排序不比符号)。"""
    tf_idx = [genes.index(t) for t in tfs]
    Z = StandardScaler().fit_transform(logx)
    n = Z.shape[0]

    pear = (Z.T @ Z) / (n - 1)
    R = np.apply_along_axis(rankdata, 0, logx)
    Rz = StandardScaler().fit_transform(R)
    spear = (Rz.T @ Rz) / (n - 1)

    # 偏相关:收缩协方差取伪逆 → 精度矩阵 → 标准化为 partial correlation
    cov = pear + np.eye(pear.shape[0]) * 0.10          # Ledoit-Wolf 式对角收缩
    prec = np.linalg.pinv(cov)
    d = np.sqrt(np.abs(np.diag(prec)))
    part = -prec / np.outer(d, d)

    return {"pearson": np.abs(pear[tf_idx]),
            "spearman": np.abs(spear[tf_idx]),
            "partial_corr": np.abs(part[tf_idx])}


def eval_grn(scores: np.ndarray, truth: np.ndarray) -> dict:
    """在 TF × target 候选边集合上评估(排除 TF 自身列)。"""
    y = truth.ravel()
    s = np.nan_to_num(scores.ravel())
    return {"auroc": float(roc_auc_score(y, s)),
            "auprc": float(average_precision_score(y, s)),
            "random_auprc": float(y.mean()),
            "n_edges": int(y.sum()), "n_candidates": int(y.size)}


def filter_thresh(adj: np.ndarray) -> np.ndarray:
    """上游 grn.py::GNInfer.filter 的 'thresh' 分支:adj < 1/adj.shape[-1] 置 0。"""
    out = adj.copy()
    out[out < (1.0 / out.shape[-1])] = 0
    return out


# =============================================================================
# 3. 守卫式官方 scPRINT 路径
# =============================================================================
SCPRINT_HINT = """
[scPRINT 官方路径未启用]
下列调用签名逐一读自本地克隆源码 C:\\Users\\fsy\\Desktop\\upstream-sources\\568_scPRINT\\:

  scprint/__init__.py L3                from .model.model import scPrint
  scprint/model/model.py L37            class scPrint(L.LightningModule, PyTorchModelHubMixin)
  docs/index.md L97-L102                model = scPrint.load_from_checkpoint(ckpt, precpt_gene_emb=None,
                                                                            transformer="normal")  # 无 triton/CPU 时必须
  scprint/tasks/__init__.py             from .grn/.cell_emb/.denoise import *
  scprint/tasks/cell_emb.py L29-L85     Embedder(batch_size=64, num_workers=8, how="random expr",
                                                 max_len=2000, doclass=True, add_zero_genes=0,
                                                 precision="16-mixed", pred_embedding=[...],
                                                 devices=[0], dtype=torch.float16,
                                                 output_expression="none")
                                        embedder(model, adata, cache=False) -> adata(.obsm["scprint"]), metrics...
  scprint/tasks/denoise.py L27-L71      Denoiser(batch_size=10, num_workers=1, max_len=5000,
                                                 how="most var", predict_depth_mult=4,
                                                 downsample=None, devices=[0])
                                        denoiser(model, adata)
  scprint/tasks/grn.py L44-L129         GNInfer(layer=None, num_genes=3000, cell_type_col="cell_type",
                                                how="random expr", preprocess="softmax",
                                                head_agg="mean", filtration="thresh", k=10,
                                                symmetrize=False, max_cells=0, devices=[0])
                                        gninfer(model, adata, cell_type=None) -> GRNAnnData
                                        注:返回的邻接阵切了 [8:, 8:],前 8 个 token 是特殊 token。
  scprint/cli.py L12                    TASKS = [("embed", Embedder), ("gninfer", GNInfer), ("denoise", Denoiser)]
  docs/index.md L86                     → 命令行 `scprint embed|gninfer|denoise --config config/[medium|large|vlarge]`
                                        (仓库 config/ 实际文件名是 pretrain_medium.yml / pretrain_large.yml /
                                         pretrain_vlarge.yml / pretrain_xlarge.yml)

依赖(见 pyproject.toml):python 3.10 · torch 2.2 · lightning · lamindb · bionty ·
scDataLoader · GRnnData · BenGRN · scib-metrics · (可选 triton 走 GPU)。
权重在 HuggingFace https://huggingface.co/jkobject/scPRINT (docs/index.md L176)。
本模块 **不代为安装**,亦不伪造这条路径的运行结果。
"""


def try_scprint(args) -> dict:
    try:
        import scprint  # noqa: F401
        from scprint import scPrint  # noqa: F401
        from scprint.tasks import Denoiser, Embedder, GNInfer  # noqa: F401
    except Exception as e:      # noqa: BLE001
        print(f"[scPRINT] import 失败:{type(e).__name__}: {e}")
        print(SCPRINT_HINT)
        return {"status": "unavailable", "reason": f"{type(e).__name__}: {e}"}
    if not args.checkpoint or not Path(args.checkpoint).exists():
        print("[scPRINT] 包已装但未提供 --checkpoint(HuggingFace jkobject/scPRINT)。")
        print(SCPRINT_HINT)
        return {"status": "no_checkpoint"}
    print("[scPRINT] 包与权重均就绪 —— 本模块只做守卫,不代跑官方推理;"
          "请按上面签名在装好的环境里调用。")
    return {"status": "ready_but_not_run", "checkpoint": str(args.checkpoint)}


# =============================================================================
# 4. 出图
# =============================================================================
def fig_cells(pcs, cell_types, counts, cv_rows, pc_pair, out):
    fig, axes = plt.subplots(1, 3, figsize=(NATURE_W2 * 1.12, 2.7))
    cats = sorted(set(cell_types))
    cols = pal(len(cats), "npg")
    a, b = pc_pair
    for c, col in zip(cats, cols):
        m = np.array(cell_types) == c
        axes[0].scatter(pcs[m, a], pcs[m, b], s=7, c=col, label=c,
                        alpha=0.75, linewidths=0)
    axes[0].set(xlabel=f"PC{a + 1}", ylabel=f"PC{b + 1}",
                title="PCA embedding\n(most cell-type-informative PCs)")
    axes[0].legend(markerscale=1.6, loc="best")

    lib = np.log10(counts.sum(axis=1) + 1)
    data = [lib[np.array(cell_types) == c] for c in cats]
    vp = axes[1].violinplot(data, showextrema=False, widths=0.85)
    for b, col in zip(vp["bodies"], cols):
        b.set_facecolor(col); b.set_alpha(0.35); b.set_edgecolor(col); b.set_linewidth(1.0)
    rng = np.random.default_rng(0)
    for i, (d, col) in enumerate(zip(data, cols), start=1):
        axes[1].scatter(i + rng.normal(0, 0.055, d.size), d, s=3.5, c=col,
                        alpha=0.5, linewidths=0)
        axes[1].hlines(np.median(d), i - 0.28, i + 0.28, color="black", lw=1.4)
    axes[1].set_xticks(range(1, len(cats) + 1)); axes[1].set_xticklabels(cats, rotation=20)
    axes[1].set(ylabel="log10 total counts", title="Library size")

    axes[2].scatter(cv_rows["fold"], cv_rows["macro_f1"], s=42,
                    c=pal(3, "npg")[1], zorder=3, linewidths=0)
    axes[2].vlines(cv_rows["fold"], 0, cv_rows["macro_f1"],
                   color=pal(3, "npg")[1], lw=1.2, alpha=0.6)
    axes[2].axhline(np.mean(cv_rows["macro_f1"]), ls="--", lw=1.0, color="grey")
    axes[2].set(xlabel="CV fold", ylabel="macro-F1", ylim=(0, 1.05),
                title="Cell-type prediction\n(PCA + kNN floor)")
    fig.tight_layout(); save_fig(fig, out); plt.close(fig)


def fig_grn_curves(curves, metrics, prevalence, out):
    fig, axes = plt.subplots(1, 3, figsize=(NATURE_W2 * 1.12, 2.8))
    cols = pal(len(curves), "npg")
    for (name, cv), col in zip(curves.items(), cols):
        axes[0].plot(cv["fpr"], cv["tpr"], lw=1.5, color=col,
                     label=f"{name} ({metrics[name]['auroc']:.3f})")
        axes[1].plot(cv["recall"], cv["precision"], lw=1.5, color=col,
                     label=f"{name} ({metrics[name]['auprc']:.3f})")
    axes[0].plot([0, 1], [0, 1], ls=":", lw=1.0, color="grey")
    axes[0].set(xlabel="False positive rate", ylabel="True positive rate", title="ROC · TF→target")
    axes[0].legend(title="AUROC", loc="lower right")
    axes[1].axhline(prevalence, ls=":", lw=1.0, color="grey")
    axes[1].set(xlabel="Recall", ylabel="Precision", title="Precision–recall · TF→target")
    axes[1].legend(title="AUPRC", loc="upper right")

    names = list(metrics)
    y = np.arange(len(names))
    vals = [metrics[n]["auprc"] for n in names]
    axes[2].hlines(y, prevalence, vals, color="#BBBBBB", lw=1.6, zorder=1)
    axes[2].scatter([prevalence] * len(y), y, s=34, c="grey", zorder=3,
                    linewidths=0, label="random")
    axes[2].scatter(vals, y, s=54, c=cols[: len(y)], zorder=3, linewidths=0)
    axes[2].set_yticks(y); axes[2].set_yticklabels(names)
    axes[2].set(xlabel="AUPRC", title="Lift over random",
                xlim=(0, max(vals) * 1.25), ylim=(-0.6, len(y) - 0.4))
    axes[2].legend(loc="upper right", markerscale=1.2)
    fig.tight_layout(); save_fig(fig, out); plt.close(fig)


def fig_grn_heatmap(truth, score, best_name, tfs, out):
    """左:推断得分热图 + 真值边就地叠加(红圈);右:真/假边得分分布对比。"""
    fig, axes = plt.subplots(1, 2, figsize=(NATURE_W2 * 1.12, 3.1),
                             gridspec_kw={"width_ratios": [2.6, 1]})
    sc = filter_thresh(score / (score.max() + 1e-12))
    im = axes[0].imshow(sc, aspect="auto", cmap=CMAP_CONT, interpolation="nearest")
    ti, tj = np.nonzero(truth)
    axes[0].scatter(tj, ti, s=22, facecolors="none", edgecolors="#E64B35",
                    linewidths=0.8, label="true edge")
    axes[0].set(title=f"Inferred TF→target · {best_name} (thresh-filtered)\n"
                      "red circles = ground-truth edges",
                xlabel="target gene", ylabel="TF")
    axes[0].set_yticks(range(len(tfs))); axes[0].set_yticklabels(tfs, fontsize=5)
    axes[0].legend(loc="upper right", markerscale=1.2, labelcolor="white")
    fig.colorbar(im, ax=axes[0], fraction=0.025, pad=0.015, label="scaled score")

    raw = score / (score.max() + 1e-12)
    grp = [raw[truth == 0], raw[truth == 1]]
    cols = [pal(4, "npg")[3], pal(4, "npg")[0]]
    vp = axes[1].violinplot(grp, showextrema=False, widths=0.85)
    for b, c in zip(vp["bodies"], cols):
        b.set_facecolor(c); b.set_alpha(0.35); b.set_edgecolor(c); b.set_linewidth(1.0)
    rng = np.random.default_rng(1)
    for i, (g, c) in enumerate(zip(grp, cols), start=1):
        s = g if g.size <= 800 else rng.choice(g, 800, replace=False)
        axes[1].scatter(i + rng.normal(0, 0.055, s.size), s, s=3, c=c, alpha=0.4, linewidths=0)
        axes[1].hlines(np.median(g), i - 0.28, i + 0.28, color="black", lw=1.4)
    axes[1].set_xticks([1, 2]); axes[1].set_xticklabels(["non-edge", "true edge"])
    axes[1].set(ylabel="scaled score", title="Score separation")
    fig.tight_layout(); save_fig(fig, out); plt.close(fig)


def fig_denoise(den, out):
    fig, axes = plt.subplots(1, 2, figsize=(NATURE_W2 * 0.82, 2.8))
    raw = [d for d in den if d["method"] == "raw (no denoising)"][0]
    knn = [d for d in den if d["method"] != "raw (no denoising)"]
    ks = [d["k"] for d in knn]
    col = pal(3, "npg")[2]
    axes[0].plot(ks, [d["poisson_nll"] for d in knn], "-o", color=col, ms=6, lw=1.4)
    axes[0].axhline(raw["poisson_nll"], ls="--", lw=1.2, color="grey")
    axes[0].text(ks[-1], raw["poisson_nll"], " raw", va="bottom", ha="right",
                 color="grey", fontsize=8)
    axes[0].set(xlabel="k (neighbours pooled)", ylabel="Poisson NLL (held-out half)",
                title="Denoising · molecular CV")

    y = np.arange(len(knn))
    delta = [raw["mse_log"] - d["mse_log"] for d in knn]
    axes[1].hlines(y, 0, delta, color="#BBBBBB", lw=1.6, zorder=1)
    axes[1].scatter(delta, y, s=54, c=col, zorder=3, linewidths=0)
    axes[1].axvline(0, ls=":", lw=1.0, color="grey")
    axes[1].set_yticks(y); axes[1].set_yticklabels([f"kNN k={k}" for k in ks])
    axes[1].set(xlabel="MSE reduction vs raw  (log-CPM)", title="Gain over no-denoising")
    fig.tight_layout(); save_fig(fig, out); plt.close(fig)


# =============================================================================
# 5. 主流程
# =============================================================================
def main() -> None:
    p = argparse.ArgumentParser(description="568 · scPRINT 基础模型三任务(可跑基线 + 守卫式官方封装)")
    p.add_argument("--counts", default=str(ROOT / "example_data" / "counts.csv"))
    p.add_argument("--meta", default=str(ROOT / "example_data" / "cell_meta.csv"))
    p.add_argument("--true_grn", default=str(ROOT / "example_data" / "true_grn_edges.csv"))
    p.add_argument("--label_key", default="cell_type")
    p.add_argument("--tf_prefix", default="TF", help="示例数据里 TF 以此前缀命名")
    p.add_argument("--n_pcs", type=int, default=20)
    p.add_argument("--folds", type=int, default=5)
    p.add_argument("--knn_k", default="3,5,10,20,40", help="去噪 kNN 的 k 网格")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--outdir", default=str(ROOT / "results"))
    p.add_argument("--assets", default=str(ROOT / "assets"))
    p.add_argument("--regen-data", action="store_true", help="重新生成 example_data/")
    p.add_argument("--run-scprint", action="store_true", help="尝试官方 scPRINT 路径(需装包+权重+GPU)")
    p.add_argument("--checkpoint", default="", help="官方 scPRINT .ckpt 路径")
    args = p.parse_args()

    np.random.seed(args.seed)
    set_pub_style(base_size=9)
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    assets = Path(args.assets); assets.mkdir(parents=True, exist_ok=True)

    ex = ROOT / "example_data"
    if args.regen_data or not (ex / "counts.csv").exists():
        print("Step 0  生成合成示例数据 (synthetic, for demo only)")
        make_synthetic(ex, seed=args.seed)

    print("Step 1  读入计数矩阵与元数据")
    cnt_df = pd.read_csv(args.counts, index_col=0)
    meta = pd.read_csv(args.meta)
    meta = meta.set_index("cell_id").loc[cnt_df.index]
    counts = cnt_df.to_numpy().astype(float)
    genes = list(cnt_df.columns)
    labels = meta[args.label_key].astype(str).to_numpy()
    tfs = [g for g in genes if g.startswith(args.tf_prefix)]
    targets = [g for g in genes if not g.startswith(args.tf_prefix)]
    print(f"        {counts.shape[0]} cells × {counts.shape[1]} genes · "
          f"{len(tfs)} TF · {len(targets)} targets · {len(set(labels))} cell types")

    logx = cpm_log1p(counts)

    # ---------------- 任务 C:细胞嵌入 + 标签预测 ----------------
    print("Step 2  任务 C · PCA 嵌入 + kNN 标签预测(交叉验证 macro-F1)")
    pcs = PCA(n_components=min(args.n_pcs, counts.shape[1] - 1),
              random_state=args.seed).fit_transform(StandardScaler().fit_transform(logx))
    skf = StratifiedKFold(n_splits=args.folds, shuffle=True, random_state=args.seed)
    cv = {"fold": [], "macro_f1": []}
    for i, (tr, te) in enumerate(skf.split(pcs, labels), start=1):
        clf = KNeighborsClassifier(n_neighbors=15).fit(pcs[tr], labels[tr])
        cv["fold"].append(i)
        cv["macro_f1"].append(float(f1_score(labels[te], clf.predict(pcs[te]), average="macro")))
    print(f"        macro-F1 = {np.mean(cv['macro_f1']):.3f} ± {np.std(cv['macro_f1']):.3f}")
    pd.DataFrame(cv).to_csv(outdir / "568_celltype_cv_macro_f1.csv", index=False)

    # 仅供【作图】挑选最能分开细胞类型的两个 PC(方差最大的 PC 被 TF 共享隐因子占据)。
    # 注意:这一步用到标签,因此只影响可视化,**绝不参与上面的 CV 指标**。
    f_stat = []
    for j in range(pcs.shape[1]):
        grp = [pcs[labels == c, j] for c in sorted(set(labels))]
        gm = pcs[:, j].mean()
        between = sum(len(g) * (g.mean() - gm) ** 2 for g in grp)
        within = sum(((g - g.mean()) ** 2).sum() for g in grp)
        f_stat.append(between / (within + 1e-12))
    pc_pair = tuple(int(v) for v in np.argsort(f_stat)[::-1][:2])
    print(f"        (作图用最有判别力的 PC: PC{pc_pair[0] + 1} / PC{pc_pair[1] + 1})")

    # ---------------- 任务 A:GRN 推断 ----------------
    print("Step 3  任务 A · 共表达 / 偏相关 打分 TF→target,对照真值网络")
    edges = pd.read_csv(args.true_grn)
    tgt_pos = {g: j for j, g in enumerate(genes)}
    truth_full = np.zeros((len(tfs), len(genes)))
    tf_pos = {t: i for i, t in enumerate(tfs)}
    for _, r in edges.iterrows():
        truth_full[tf_pos[r["tf"]], tgt_pos[r["target"]]] = 1
    tcols = [tgt_pos[g] for g in targets]
    truth = truth_full[:, tcols]                     # 候选边只取 TF × non-TF,排除 TF-TF 与自身

    all_scores = score_grn(logx, genes, tfs)
    grn_metrics, curves = {}, {}
    for name, S in all_scores.items():
        Ssub = S[:, tcols]
        grn_metrics[name] = eval_grn(Ssub, truth)
        fpr, tpr, _ = roc_curve(truth.ravel(), np.nan_to_num(Ssub.ravel()))
        pr, rc, _ = precision_recall_curve(truth.ravel(), np.nan_to_num(Ssub.ravel()))
        curves[name] = {"fpr": fpr, "tpr": tpr, "precision": pr, "recall": rc}
        m = grn_metrics[name]
        print(f"        {name:14s} AUROC={m['auroc']:.3f}  AUPRC={m['auprc']:.3f} "
              f"(random {m['random_auprc']:.3f})")
    best = max(grn_metrics, key=lambda n: grn_metrics[n]["auprc"])
    pd.DataFrame(grn_metrics).T.to_csv(outdir / "568_grn_benchmark.csv")

    # 导出 top 边(把打分变成可用产物,不只是指标)
    Sb = all_scores[best][:, tcols]
    ii, jj = np.unravel_index(np.argsort(Sb, axis=None)[::-1][:200], Sb.shape)
    hit = truth[ii, jj].astype(int)
    pd.DataFrame({"tf": [tfs[i] for i in ii], "target": [targets[j] for j in jj],
                  "score": Sb[ii, jj], "in_true_grn": hit}).to_csv(
        outdir / "568_grn_top_edges.csv", index=False)

    # precision@k:排在最前面的边有多少是真的 —— 比 AUPRC 更贴近"我只做前 N 条去验证"的实际用法
    prec_at_k = {f"precision@{k}": float(hit[:k].mean()) for k in (20, 50, 100, 200)}
    print("        " + "  ".join(f"{k}={v:.2f}" for k, v in prec_at_k.items()))

    # ---------------- 任务 B:去噪 ----------------
    print("Step 4  任务 B · 分子交叉验证(binomial split)+ kNN 池化去噪")
    rs = np.random.RandomState(args.seed)
    Xh, Yh = split_molecules(counts.astype(int), 0.5, 0.0, rs)
    y_tot = Yh.sum(axis=1, keepdims=True); y_tot[y_tot == 0] = 1.0

    def score_pred(pred: np.ndarray) -> dict:
        p_tot = pred.sum(axis=1, keepdims=True); p_tot[p_tot == 0] = 1.0
        scaled = pred / p_tot * y_tot                     # 上游 open_benchmark 的深度对齐
        return {"poisson_nll": poisson_nll_loss(scaled + 1e-8, Yh),
                "mse_log": float(mean_squared_error(cpm_log1p(Yh), cpm_log1p(pred)))}

    den = [{"method": "raw (no denoising)", "k": 1, **score_pred(Xh.astype(float))}]
    for k in [int(v) for v in args.knn_k.split(",")]:
        den.append({"method": f"kNN pooling k={k}", "k": k,
                    **score_pred(knn_pool_denoise(Xh.astype(float), k=k, seed=args.seed))})
    for d in den:
        print(f"        {d['method']:22s} PoissonNLL={d['poisson_nll']:.4f}  MSE(log)={d['mse_log']:.4f}")
    pd.DataFrame(den).to_csv(outdir / "568_denoise_benchmark.csv", index=False)
    best_den = min(den, key=lambda d: d["poisson_nll"])

    # ---------------- 官方路径守卫 ----------------
    print("Step 5  官方 scPRINT 路径")
    sc_status = try_scprint(args) if args.run_scprint else {"status": "not_requested"}
    if not args.run_scprint:
        print("        未请求(加 --run-scprint 查看真实 API 签名与依赖清单)")

    # ---------------- 出图 ----------------
    print("Step 6  出图(lollipop / violin / 曲线 / heatmap,无条形图)")
    fig_cells(pcs, labels, counts, cv, pc_pair, assets / "fig1_cells_embedding_and_labels")
    fig_grn_curves(curves, grn_metrics, float(truth.mean()), assets / "fig2_grn_benchmark")
    fig_grn_heatmap(truth, all_scores[best][:, tcols], best, tfs, assets / "fig3_grn_heatmap")
    fig_denoise(den, assets / "fig4_denoise_molecular_cv")
    for f in assets.glob("fig*.p*"):
        (outdir / f.name).write_bytes(f.read_bytes())

    summary = {
        "module": "568_scprint_foundation_grn",
        "upstream": {"repo": "https://github.com/jkobject/scPRINT",
                     "paper": "Kalfon J et al. Nat Commun 2025;16:3607",
                     "pmid": "40240364", "doi": "10.1038/s41467-025-58699-1"},
        "params": {k: v for k, v in vars(args).items()},
        "data": {"n_cells": int(counts.shape[0]), "n_genes": int(counts.shape[1]),
                 "n_tf": len(tfs), "n_targets": len(targets),
                 "n_true_edges": int(truth.sum())},
        "task_C_cell_embedding": {"macro_f1_mean": float(np.mean(cv["macro_f1"])),
                                  "macro_f1_sd": float(np.std(cv["macro_f1"]))},
        "task_A_grn": {"metrics": grn_metrics, "best_scorer": best,
                       "precision_at_k": prec_at_k},
        "task_B_denoise": {"rows": den, "best": best_den["method"]},
        "scprint_official": sc_status,
        # 依赖版本快照(铁律6 可复现):关键数字全部由本次运行生成,不手填
        "session_info": {
            "python": sys.version.split()[0],
            "numpy": np.__version__, "pandas": pd.__version__,
            "scipy": __import__("scipy").__version__,
            "scikit-learn": __import__("sklearn").__version__,
            "matplotlib": __import__("matplotlib").__version__,
        },
    }
    (outdir / "568_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n完成 · 指标与图见 {outdir}  展示图见 {assets}")


if __name__ == "__main__":
    main()
