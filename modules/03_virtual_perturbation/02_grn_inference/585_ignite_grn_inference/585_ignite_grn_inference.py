# -*- coding: utf-8 -*-
"""585 · IGNITE — 动力学 Ising 反问题推断基因调控网络(GRN)+ 敲除预测.

IGNITE = Inference of Gene Networks using Inverse kinetic Theory and Experiments。
它把伪时序排序后的单细胞表达二值化成自旋 (±1),用**非对称动力学 Ising 模型的
反问题**(asynchronous Glauber dynamics + 最大似然 + L1/L2 正则)学出耦合矩阵 J
和外场 h,再用学到的网络模拟单/多基因敲除。

本模块结构(与库规矩一致:任何"更好"的方法都必须带一个本机能跑的朴素对照):
  · 基线(永远可跑,只用本机已有依赖):三个朴素 GRN 推断对照
      B1 Pearson 相关网络         —— 最经典也最朴素的共表达网络
      B2 GraphicalLassoCV 偏相关  —— 稀疏精度矩阵,去掉间接相关
      B3 滞后岭回归 s(t+1)~s(t)   —— 非对称线性耦合,IGNITE 动力学的线性对应物
    三者都在同一份"已知真值网络"的合成数据上用 AUROC/AUPRC 评边恢复。
  · IGNITE 路径(守卫式引用封装,--ignite-repo 指向本地 clone):
    IGNITE 没有 pip/conda 包,上游是仓库里的 notebook + lib/ 源码,只能克隆后用。

上游(已核实,见 README):
  论文 Corridori et al., PLoS Comput Biol 2026, doi:10.1371/journal.pcbi.1014067, PMID 41984780
  仓库 https://github.com/CleliaCorridori/IGNITE
真实 API 逐条核对自上游源码本体(2026-07-21 复核,非二手文档):
  IGNITE_and_SCODE_notebooks/lib/ml_wrapper.py   —— asynch_reconstruction 类及其方法
  IGNITE_and_SCODE_notebooks/lib/fun_asynch.py   —— 输入形状/±1 约定/J 的方向
  IGNITE_and_SCODE_notebooks/lib/funcs_ko.py     —— 敲除模拟函数
  IGNITE_and_SCODE_notebooks/lib/funcs_IsingPars.py —— grid_search 超参搜索
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

# 批处理脚本统一走非交互后端。这不只是习惯问题:上游 asynch_reconstruction.reconstruct()
# 内部会调 plot_fields_and_couplings() → plt.show()(见 lib/ml_wrapper.py:66 与 :122),
# 在交互后端下会弹窗并阻塞整个脚本。Agg 下退化成一条 UserWarning,不阻塞。
import matplotlib
matplotlib.use("Agg")

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")
ASSETS = os.path.join(HERE, "assets")
EXAMPLE = os.path.join(HERE, "example_data")

# 出图统一走框架样式
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..", "..", "_framework")))
from pubstyle import (  # noqa: E402
    set_pub_style, pal, save_fig, CMAP_CONT, CMAP_DIVERGE, NATURE_W1, NATURE_W2,
)


# =============================================================================
# 合成示例数据:已知真值的非对称调控网络 + 异步 Glauber 轨迹
# -----------------------------------------------------------------------------
# 注意:这里的模拟器**只用于造带 ground truth 的 demo 数据**,不是 IGNITE 的推断,
# 也不是对 IGNITE 的复现。IGNITE 做的是反向问题(从数据学 J),这里是正向出数据。
# =============================================================================
def simulate_kinetic_ising(n_genes: int = 24, n_steps: int = 4000, density: float = 0.15,
                           delta_t: float = 0.1, seed: int = 0):
    """正向模拟:给定稀疏非对称 J 与外场 h,跑异步 Glauber 动力学出 ±1 自旋轨迹。

    返回 spins (n_genes × n_steps, 取值 ±1) 与真值 J (n_genes × n_genes)。
    """
    rng = np.random.default_rng(seed)
    J = np.zeros((n_genes, n_genes))
    n_edges = int(density * n_genes * (n_genes - 1))
    # 非对称:调控是有向的,J[i, j] = j 对 i 的作用,J[i,j] != J[j,i]
    idx = [(i, j) for i in range(n_genes) for j in range(n_genes) if i != j]
    for k in rng.choice(len(idx), size=n_edges, replace=False):
        i, j = idx[k]
        J[i, j] = rng.choice([-1.0, 1.0]) * rng.uniform(0.6, 1.4)
    h = rng.normal(0, 0.25, n_genes)

    spins = np.empty((n_genes, n_steps))
    s = rng.choice([-1.0, 1.0], n_genes)
    for t in range(n_steps):
        # 异步更新:每步只翻一个基因,符合 Glauber 异步约定
        i = rng.integers(n_genes)
        field = h[i] + J[i] @ s
        p_up = 1.0 / (1.0 + np.exp(-2.0 * field))
        s = s.copy()
        s[i] = 1.0 if rng.random() < p_up else -1.0
        spins[:, t] = s
    return spins, J, h


def write_example_data(path: str, seed: int = 0):
    """把合成自旋轨迹与真值网络写进 example_data/(文件头标注 synthetic)。"""
    spins, J, h = simulate_kinetic_ising(seed=seed)
    genes = [f"Gene{i:02d}" for i in range(spins.shape[0])]
    os.makedirs(path, exist_ok=True)

    expr = os.path.join(path, "spins_pseudotime_ordered.csv")
    with open(expr, "w", encoding="utf-8") as fh:
        fh.write("# synthetic, for demo only -- kinetic Ising simulation, not real scRNA-seq\n")
        fh.write("# rows = genes, columns = pseudotime-ordered cells, values = binarized spin (+1/-1)\n")
        pd.DataFrame(spins, index=genes,
                     columns=[f"Cell{t:05d}" for t in range(spins.shape[1])]).to_csv(fh)

    truth = os.path.join(path, "true_network.csv")
    with open(truth, "w", encoding="utf-8") as fh:
        fh.write("# synthetic, for demo only -- ground-truth asymmetric coupling matrix J[target, regulator]\n")
        pd.DataFrame(J, index=genes, columns=genes).to_csv(fh)
    return expr, truth


# =============================================================================
# 基线:三个朴素 GRN 推断对照(全部只用本机已有依赖)
# =============================================================================
def baseline_correlation(spins: np.ndarray) -> np.ndarray:
    """B1 Pearson 共表达网络:最朴素的对照,对称、无方向、含大量间接边。"""
    C = np.corrcoef(spins)
    return np.nan_to_num(C)


def baseline_partial_correlation(spins: np.ndarray) -> np.ndarray:
    """B2 GraphicalLassoCV 稀疏偏相关:从精度矩阵读条件独立,压掉间接相关。"""
    from sklearn.covariance import GraphicalLassoCV
    X = spins.T  # sklearn 要 样本 × 特征
    X = X + np.random.default_rng(0).normal(0, 1e-3, X.shape)  # 破坏二值退化,保证协方差可逆
    try:
        model = GraphicalLassoCV(alphas=4, max_iter=200).fit(X)
        P = -model.precision_
    except Exception as e:  # 数值不收敛时退回经验精度矩阵,不让基线整体失败
        print(f"       [B2] GraphicalLassoCV failed ({type(e).__name__}), falling back to pinv(cov)")
        P = -np.linalg.pinv(np.cov(spins))
    d = np.sqrt(np.abs(np.diag(P)))
    P = P / np.outer(d, d)
    np.fill_diagonal(P, 0.0)
    return np.nan_to_num(P)


def baseline_lagged_ridge(spins: np.ndarray, alpha: float = 1.0) -> np.ndarray:
    """B3 滞后岭回归 s(t+1) ~ s(t):非对称线性耦合。

    这是本模块最关键的对照——它和 IGNITE 用的是同一份时序信息、同一个"上一时刻
    决定下一时刻"的因果方向,只是把 Glauber 的非线性最大似然换成线性最小二乘。
    IGNITE 若要声称有增益,必须先赢过它,而不是只赢过 Pearson。
    """
    from sklearn.linear_model import Ridge
    X = spins[:, :-1].T   # (T-1) × n_genes,时刻 t
    Y = spins[:, 1:].T    # (T-1) × n_genes,时刻 t+1
    W = Ridge(alpha=alpha, fit_intercept=True).fit(X, Y).coef_  # (n_genes × n_genes)
    np.fill_diagonal(W, 0.0)
    return np.nan_to_num(W)


def score_edges(J_true: np.ndarray, W: np.ndarray) -> dict:
    """按 |权重| 排序做边恢复评估;只看非对角元素,真值边定义为 J_true != 0。"""
    from sklearn.metrics import roc_auc_score, average_precision_score, roc_curve
    off = ~np.eye(J_true.shape[0], dtype=bool)
    y = (J_true[off] != 0).astype(int)
    s = np.abs(W[off])
    fpr, tpr, _ = roc_curve(y, s)
    return {
        "auroc": float(roc_auc_score(y, s)),
        "auprc": float(average_precision_score(y, s)),
        "prevalence": float(y.mean()),
        "_fpr": fpr, "_tpr": tpr,
    }


# =============================================================================
# IGNITE 路径:守卫式引用封装
# =============================================================================
def run_ignite(spins: np.ndarray, genes: list[str], repo: str | None,
               delta_t: float = 0.1, lam: float = 0.01, n_epochs: int = 200):
    """尝试调用真实 IGNITE。

    IGNITE **没有 pip / conda 包**:上游是一个仓库,推断实现放在
    `IGNITE_and_SCODE_notebooks/lib/`,要用必须先 clone 再把该 lib 目录的父目录
    加进 sys.path。传 --ignite-repo <本地clone路径> 启用。

    已核实的真实 API(逐条对照上游源码,复核日期 2026-07-21;行号为上游 main 分支):
        from lib.ml_wrapper import asynch_reconstruction        # ml_wrapper.py:7 (class)
        model = asynch_reconstruction(x, delta_t, LAMBDA, MOM=None, gamma=1,
                                      opt='NADAM', reg='L1', ax_names=[])  # :14 __init__
        model.reconstruct(x, Nepochs, start_lr=1, drop=0.99, edrop=20)     # :50
        model.h            # 外场向量,__init__ :38 初始化,reconstruct 后被赋值
        model.J            # 耦合矩阵 Nvar × Nvar,__init__ :33
        model.generate_samples(t_size=None, seed=1)                        # :124
        model.find_likelihood(x)                                           # :69
        model.load_parameters(h, J, plot=False)                            # :74
    J 的方向:上游 fun_asynch.find_theta() 算的是 `h + J @ x`(fun_asynch.py:77),
    即 theta_i = h_i + Σ_j J[i,j] x_j,所以 **J[i, j] = 调控者 j 对靶基因 i 的作用**,
    与本模块 example_data/true_network.csv 的 J[target, regulator] 约定一致。
    输入 x 的形状/取值(fun_asynch.py:80-104 asynch_glauber_dynamics 的 docstring 与
    generate_samples_asynch 的 ±1 初始化):**Nvar × Nsteps 的 ±1 自旋矩阵**,
    列按伪时序排列。敲除模拟见上游 lib/funcs_ko.py 的 KO_wrap(:24)/ info_KO(:58)/
    KO_avg_weighted(:86)/ KO_diff_sim(:117)。

    ⚠️ 上游副作用:reconstruct() 会打印逐 epoch 表格,并在结束时调
    plot_fields_and_couplings() 画一张 h/J 图再 plt.show()(ml_wrapper.py:66, :122)。
    本模块已在文件头强制 Agg 后端,并在调用后 plt.close("all") 清掉这张图。

    上游仓库无 requirements 文件、也无 LICENSE 文件(2026-07-21 查 GitHub 仓库树确认),
    依赖(numba)与超参数标定以官方 notebook `01_IGNITE_InferenceMethod.ipynb` 为准,
    此处**未固定**。
    """
    install_hint = (
        "IGNITE 无 pip 包;用法:\n"
        "  git clone https://github.com/CleliaCorridori/IGNITE.git\n"
        "  python 585_ignite_grn_inference.py --ignite-repo /path/to/IGNITE\n"
        "  (推断实现在 IGNITE_and_SCODE_notebooks/lib/,需要 numba)"
    )
    if not repo:
        return {"status": "skipped", "reason": "未提供 --ignite-repo", "how_to_run": install_hint}

    libdir = os.path.join(repo, "IGNITE_and_SCODE_notebooks")
    if not os.path.isdir(os.path.join(libdir, "lib")):
        return {"status": "skipped",
                "reason": f"{libdir}/lib 不存在,--ignite-repo 似乎不是 IGNITE 的 clone",
                "how_to_run": install_hint}

    sys.path.insert(0, libdir)
    try:
        from lib.ml_wrapper import asynch_reconstruction
    except ImportError as e:
        return {"status": "skipped",
                "reason": f"导入 lib.ml_wrapper 失败 ({e});上游依赖 numba,本机可能缺",
                "how_to_run": install_hint}

    # 真实调用。签名已核实,但超参数(LAMBDA/Nepochs/delta_t)的标定以官方 notebook 为准。
    # 上游用 numba njit,输入必须是 float64;ax_names 长度必须等于基因数(会当 tick label)。
    import matplotlib.pyplot as plt
    x = np.ascontiguousarray(spins, dtype=np.float64)
    try:
        model = asynch_reconstruction(x, delta_t, lam, opt="NADAM", reg="L1", ax_names=list(genes))
        model.reconstruct(x, n_epochs)
    except Exception as e:  # 不让上游报错拖垮整个模块——基线结果照常写出
        plt.close("all")
        return {"status": "failed", "reason": f"{type(e).__name__}: {e}",
                "how_to_run": install_hint}
    plt.close("all")  # 清掉 reconstruct() 内部 plot_fields_and_couplings() 留下的图
    return {"status": "ok", "J": np.asarray(model.J), "h": np.asarray(model.h)}


# =============================================================================
# 出图(全部非条形图:heatmap / ROC 曲线 / dumbbell / 散点 / violin)
# =============================================================================
def fig_networks(J_true, mats: dict, outstem: str):
    """Fig 1 · 真值网络与各基线推断矩阵并排 heatmap。"""
    import matplotlib.pyplot as plt
    items = [("Ground truth J", J_true)] + list(mats.items())
    fig, axes = plt.subplots(1, len(items), figsize=(NATURE_W2 * 1.05, 2.5))
    for ax, (name, M) in zip(np.ravel(axes), items):
        v = np.percentile(np.abs(M[~np.eye(M.shape[0], dtype=bool)]), 98) or 1.0
        im = ax.imshow(M, cmap=CMAP_DIVERGE, vmin=-v, vmax=v, interpolation="nearest")
        ax.set_title(name, fontsize=8)
        ax.set_xlabel("Regulator", fontsize=7)
        if ax is np.ravel(axes)[0]:
            ax.set_ylabel("Target", fontsize=7)
        ax.set_xticks([]); ax.set_yticks([])
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03).ax.tick_params(labelsize=6)
    fig.suptitle("Inferred gene-regulatory couplings vs ground truth", fontsize=9, y=1.06)
    save_fig(fig, outstem)
    plt.close(fig)


def fig_roc_and_auc(scores: dict, outstem: str):
    """Fig 2 · 左:边恢复 ROC 曲线;右:AUROC/AUPRC dumbbell(替代条形图)。"""
    import matplotlib.pyplot as plt
    # wspace 拉开:右图的方法名较长,否则会压到左图 ROC 区域
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(NATURE_W2 * 1.1, 3.0),
                                   gridspec_kw={"wspace": 0.75})
    colors = pal(len(scores), "npg")

    for (name, sc), c in zip(scores.items(), colors):
        ax1.plot(sc["_fpr"], sc["_tpr"], lw=1.6, color=c,
                 label=f"{name} (AUROC {sc['auroc']:.3f})")
    ax1.plot([0, 1], [0, 1], ls="--", lw=1.0, color="#999999")
    ax1.set_xlabel("False positive rate"); ax1.set_ylabel("True positive rate")
    ax1.set_title("Edge recovery", fontsize=9)
    ax1.legend(fontsize=6, loc="lower right")

    names = list(scores)
    ypos = np.arange(len(names))
    for y, name, c in zip(ypos, names, colors):
        a, p = scores[name]["auroc"], scores[name]["auprc"]
        ax2.plot([min(a, p), max(a, p)], [y, y], lw=1.4, color="#CCCCCC", zorder=1)
        ax2.scatter(a, y, s=52, color=c, zorder=3, label="AUROC" if y == 0 else None)
        ax2.scatter(p, y, s=52, color=c, zorder=3, marker="D",
                    facecolors="white", edgecolors=c, linewidths=1.4,
                    label="AUPRC" if y == 0 else None)
    prev = list(scores.values())[0]["prevalence"]
    ax2.axvline(0.5, ls="--", lw=1.0, color="#999999")
    ax2.axvline(prev, ls=":", lw=1.0, color="#666666")
    ax2.text(prev, len(names) - 0.35, "edge prevalence", fontsize=6, rotation=90,
             ha="right", va="top", color="#666666")
    ax2.set_yticks(ypos); ax2.set_yticklabels(names, fontsize=7)
    ax2.set_ylim(-0.6, len(names) - 0.4)   # 留白,避免首尾标记贴边
    ax2.set_xlabel("Score"); ax2.set_xlim(0, 1.02)
    ax2.set_title("AUROC vs AUPRC", fontsize=9)
    # 图例放左下(分数都偏右侧,左下是空区),避免压住 B1 的标记
    ax2.legend(fontsize=6, loc="lower left")
    save_fig(fig, outstem)
    plt.close(fig)


def fig_weight_agreement(J_true, W, name: str, outstem: str):
    """Fig 3 · 左:真值耦合 vs 推断权重散点;右:真边/非边的权重 violin。"""
    import matplotlib.pyplot as plt
    off = ~np.eye(J_true.shape[0], dtype=bool)
    t, w = J_true[off], W[off]
    is_edge = t != 0
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(NATURE_W2, 3.0))
    c = pal(3, "npg")

    ax1.scatter(t[~is_edge], w[~is_edge], s=9, alpha=0.35, color="#BBBBBB",
                label="non-edge", edgecolors="none")
    ax1.scatter(t[is_edge], w[is_edge], s=14, alpha=0.85, color=c[0],
                label="true edge", edgecolors="none")
    ax1.axhline(0, lw=0.8, color="#999999"); ax1.axvline(0, lw=0.8, color="#999999")
    if is_edge.sum() > 2:
        from scipy.stats import spearmanr
        rho = spearmanr(t[is_edge], w[is_edge]).correlation
        ax1.set_title(f"{name}: sign & magnitude recovery (rho = {rho:.2f})", fontsize=8)
    ax1.set_xlabel("Ground-truth coupling"); ax1.set_ylabel("Inferred weight")
    ax1.legend(fontsize=6)

    parts = ax2.violinplot([np.abs(w[~is_edge]), np.abs(w[is_edge])],
                           showextrema=False, showmedians=True, widths=0.8)
    for body, col in zip(parts["bodies"], ["#BBBBBB", c[0]]):
        body.set_facecolor(col); body.set_alpha(0.65); body.set_edgecolor("black")
    parts["cmedians"].set_color("black")
    for i, arr in enumerate([np.abs(w[~is_edge]), np.abs(w[is_edge])], start=1):
        jit = np.random.default_rng(0).normal(i, 0.045, min(len(arr), 400))
        ax2.scatter(jit, np.random.default_rng(1).choice(arr, len(jit), replace=False),
                    s=4, alpha=0.3, color="black", edgecolors="none", zorder=3)
    ax2.set_xticks([1, 2]); ax2.set_xticklabels(["Non-edge", "True edge"])
    ax2.set_ylabel("|Inferred weight|")
    ax2.set_title("Weight separation", fontsize=8)
    save_fig(fig, outstem)
    plt.close(fig)


# =============================================================================
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--spins", default=os.path.join(EXAMPLE, "spins_pseudotime_ordered.csv"),
                    help="基因 × 伪时序细胞 的 ±1 自旋矩阵 CSV(行=基因,列=细胞)")
    ap.add_argument("--truth", default=os.path.join(EXAMPLE, "true_network.csv"),
                    help="可选:真值耦合矩阵 CSV,用于评边恢复;真实数据通常没有")
    ap.add_argument("--ignite-repo", default=None,
                    help="本地 IGNITE clone 路径;不给则只跑基线")
    ap.add_argument("--delta-t", type=float, default=0.1, help="IGNITE 时间步长")
    ap.add_argument("--lam", type=float, default=0.01, help="IGNITE 正则强度 LAMBDA")
    ap.add_argument("--epochs", type=int, default=200, help="IGNITE 重构轮数 Nepochs")
    ap.add_argument("--outdir", default=RESULTS)
    ap.add_argument("--figdir", default=None,
                    help="出图目录;默认:跑默认 outdir 时写 assets/,自定义 outdir 时写 outdir")
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()

    np.random.seed(a.seed)
    set_pub_style(base_size=9)
    os.makedirs(a.outdir, exist_ok=True)
    # 自定义 outdir 的运行不应覆盖仓库里提交的展示图
    figdir = a.figdir or (ASSETS if os.path.abspath(a.outdir) == os.path.abspath(RESULTS)
                          else a.outdir)
    os.makedirs(figdir, exist_ok=True)

    print("[585] Step 1 · 载入输入")
    if not os.path.exists(a.spins):
        print("       example_data 缺失,先生成合成数据")
        write_example_data(EXAMPLE, seed=a.seed)
    df = pd.read_csv(a.spins, index_col=0, comment="#")
    spins = df.to_numpy(dtype=float)
    genes = list(df.index)
    print(f"       {spins.shape[0]} genes × {spins.shape[1]} pseudotime-ordered cells")
    uniq = np.unique(spins)
    if not np.all(np.isin(uniq, [-1.0, 1.0])):
        print("       输入非 ±1,按中位数二值化(IGNITE 要求 ±1 自旋)")
        spins = np.where(spins > np.median(spins, axis=1, keepdims=True), 1.0, -1.0)

    J_true = None
    if a.truth and os.path.exists(a.truth):
        J_true = pd.read_csv(a.truth, index_col=0, comment="#").to_numpy(dtype=float)

    print("[585] Step 2 · 基线 GRN 推断(3 个朴素对照)")
    mats = {
        "B1 Pearson": baseline_correlation(spins),
        "B2 Partial corr (glasso)": baseline_partial_correlation(spins),
        "B3 Lagged ridge": baseline_lagged_ridge(spins),
    }
    for k in mats:
        print(f"       {k}: done")

    print("[585] Step 3 · IGNITE 路径")
    ig = run_ignite(spins, genes, a.ignite_repo, a.delta_t, a.lam, a.epochs)
    if ig["status"] == "ok":
        mats["IGNITE (kinetic Ising)"] = ig["J"]
        pd.DataFrame(ig["J"], index=genes, columns=genes).to_csv(
            os.path.join(a.outdir, "ignite_J.csv"))
        pd.Series(ig["h"], index=genes).to_csv(os.path.join(a.outdir, "ignite_h.csv"))
        print("       IGNITE 推断完成,已写出 J / h")
    else:
        print(f"       skipped: {ig['reason']}")
        print("       " + ig["how_to_run"].replace("\n", "\n       "))

    print("[585] Step 4 · 写出网络矩阵")
    for name, M in mats.items():
        slug = name.split(" ")[0].lower().strip("()")
        pd.DataFrame(M, index=genes, columns=genes).to_csv(
            os.path.join(a.outdir, f"network_{slug}.csv"))

    summary = {"n_genes": len(genes), "n_cells": int(spins.shape[1]),
               "ignite": {k: v for k, v in ig.items() if k not in ("J", "h")},
               "methods": list(mats)}

    print("[585] Step 5 · 评估与出图")
    if J_true is not None:
        scores = {name: score_edges(J_true, M) for name, M in mats.items()}
        for name, sc in scores.items():
            print(f"       {name}: AUROC {sc['auroc']:.3f} · AUPRC {sc['auprc']:.3f}")
        summary["edge_recovery"] = {
            k: {kk: vv for kk, vv in v.items() if not kk.startswith("_")}
            for k, v in scores.items()}
        pd.DataFrame(summary["edge_recovery"]).T.to_csv(
            os.path.join(a.outdir, "edge_recovery_scores.csv"))
        fig_roc_and_auc(scores, os.path.join(figdir, "585_fig2_edge_recovery"))
        best = max(scores, key=lambda k: scores[k]["auroc"])
        fig_weight_agreement(J_true, mats[best], best,
                             os.path.join(figdir, "585_fig3_weight_agreement"))
        summary["best_method_by_auroc"] = best
        fig_networks(J_true, mats, os.path.join(figdir, "585_fig1_networks"))
    else:
        print("       未提供真值网络,跳过边恢复评估,只出网络 heatmap")
        fig_networks(np.zeros_like(list(mats.values())[0]), mats,
                     os.path.join(figdir, "585_fig1_networks"))

    with open(os.path.join(a.outdir, "585_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=1, ensure_ascii=False, default=str)
    print(f"[585] done · results -> {a.outdir} · figures -> {figdir}")


if __name__ == "__main__":
    main()
