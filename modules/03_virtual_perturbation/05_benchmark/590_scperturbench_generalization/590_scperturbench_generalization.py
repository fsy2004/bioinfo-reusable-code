# -*- coding: utf-8 -*-
"""590 · scPerturBench — 单细胞扰动响应预测的评测尺 (evaluation metrics wrapper).

本模块**不安装、不复现** scPerturBench 里的 27 个预测方法(那需要 Podman 镜像 + GPU)。
它复用的是该 benchmark 的**评测层**:把「你自己的扰动预测」按 scPerturBench 的同一套
指标打分,并**强制和朴素基线对照**(controlBaseline / trainMean),回答唯一重要的问题:
    你的深度扰动模型,真的赢过"什么都不预测"和"预测训练集均值"吗?

上游: Wei Z, Wang Y, Gao Y, ..., Liu Q. Benchmarking algorithms for generalizable
      single-cell perturbation response prediction. Nature Methods 2026 Feb.
      PMID 41381899 · doi:10.1038/s41592-025-02980-0
      repo https://github.com/bm2-lab/scPerturBench

指标口径来源(实际读过的源码 URL,均已核对):
  - scPerturBench 打分脚本
    repos/bm2-lab/scPerturBench/contents/Cellular_context_generalization/o.o.d./calPerformance.py
    repos/bm2-lab/scPerturBench/contents/Perturbation_generalization/calPerformance_genetic.py
    → 真实调用: pt.tools.Distance(metric=<m>, layer_key='X')
                 .onesided_distances(adata, groupby=<col>, selected_group='imputed',
                                     groups=['stimulated'])
      metrics = ['mse','pearson_distance','edistance','sym_kldiv','wasserstein']
      pearson_distance 前先 calculateDelta(减 control 均值);sym_kldiv 后取 log2(x+1);
      **仅** edistance / wasserstein / sym_kldiv(+mean_var_distribution)这几个分布型指标
      会先把各组 subsample 至 2000 细胞(np.random.seed(42)),mse / pearson_distance 用全量;
      先按预计算的 top-N DEG(100 / 5000)子集基因。
  - 朴素基线
    Cellular_context_generalization/o.o.d./baseControl.py  (方法名 controlMean)
    Cellular_context_generalization/o.o.d./trainMean.py    (generateExp:高斯采样,非广播均值)
  - pertpy 指标定义 (公式逐行对齐)
    repos/scverse/pertpy/contents/pertpy/tools/_distances/_distances.py

本机没有 pertpy(见 --use-pertpy 守卫路径),故默认走**本地实现**:公式逐行照抄上述
pertpy 源码,不是臆造。唯一无法本地等价的是 wasserstein —— pertpy 用 ott-jax 的
Sinkhorn 熵正则 OT,本地改用 sliced_wasserstein(随机投影 1D W2),**指标名刻意不同**,
数值不可与上游论文表格直接比较。装了 pertpy 时加 --use-pertpy 可拿到上游同名指标。

用法:
    python 590_scperturbench_generalization.py
    python 590_scperturbench_generalization.py --observed my_obs.csv --predicted my_pred.csv \
        --n-deg 100 --outdir results/run1
"""
from __future__ import annotations

import argparse
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist, pdist
from scipy.stats import pearsonr

warnings.filterwarnings("ignore")

ROOT = Path(__file__).parent
sys.path.insert(0, str(ROOT.parents[2] / "_framework"))
from pubstyle import (CMAP_CONT, NATURE_W1, NATURE_W2, pal, save_fig,  # noqa: E402
                      set_pub_style)

SEED = 42

# 距离型指标(越小越好) vs 相似型指标(越大越好)
LOWER_IS_BETTER = {
    "mse": True,
    "pearson_distance": True,
    "edistance": True,
    "sym_kldiv": True,
    "sliced_wasserstein": True,
    "wasserstein": True,
    "common_deg": False,
}
DEFAULT_METRICS = ["mse", "pearson_distance", "edistance", "sym_kldiv",
                   "sliced_wasserstein", "common_deg"]

# 上游 calPerformance.py:109 / calPerformance_genetic.py:110 的分支条件:
#   if doSubSample and metric in ['edistance','wasserstein','mean_var_distribution','sym_kldiv']
# 其余指标(mse / pearson_distance)在上游用全量细胞,不 subsample。
SUBSAMPLED_METRICS = {"edistance", "wasserstein", "sliced_wasserstein",
                      "mean_var_distribution", "sym_kldiv"}


# =============================================================================
# 指标实现 —— 公式逐行对齐 pertpy/tools/_distances/_distances.py
# =============================================================================
def _pairwise_mean_within(X: np.ndarray) -> float:
    """组内平均两两欧氏距离(不含对角线)。对齐 pertpy _euclidean_pairwise_mean_within。"""
    if X.shape[0] < 2:
        return 0.0
    return float(pdist(X, metric="euclidean").mean())


def _pairwise_mean_between(X: np.ndarray, Y: np.ndarray) -> float:
    """组间平均两两欧氏距离(全矩阵均值)。对齐 pertpy _euclidean_pairwise_mean_between。"""
    return float(cdist(X, Y, metric="euclidean").mean())


def m_edistance(X: np.ndarray, Y: np.ndarray) -> float:
    """Energy distance: 2*between - within_X - within_Y (pertpy Edistance.__call__)。"""
    return 2 * _pairwise_mean_between(X, Y) - _pairwise_mean_within(X) - _pairwise_mean_within(Y)


def m_mse(X: np.ndarray, Y: np.ndarray) -> float:
    """pseudobulk 均值向量的 L2 范数平方 / 基因数 (pertpy MeanSquaredDistance)。"""
    return float(np.linalg.norm(X.mean(axis=0) - Y.mean(axis=0), ord=2) ** 2 / X.shape[1])


def m_pearson_distance(X: np.ndarray, Y: np.ndarray) -> float:
    """1 - pearsonr(pseudobulk_X, pseudobulk_Y) (pertpy PearsonDistance)。

    上游在调用它之前会先做 calculateDelta(减去 control 均值),所以这里传入的应当
    已经是 delta 矩阵 —— 即论文里的 PCC-delta。
    """
    a, b = X.mean(axis=0), Y.mean(axis=0)
    if np.std(a) == 0 or np.std(b) == 0:      # 常数预测 → 相关系数无定义
        return float("nan")
    return float(1 - pearsonr(a, b)[0])


def m_sym_kldiv(X: np.ndarray, Y: np.ndarray, epsilon: float = 1e-8) -> float:
    """逐基因高斯假设下的对称 KL,再对基因取平均 (pertpy SymmetricKLDivergence)。

    上游 calPerformance 在拿到该值后再取 log2(x+1),本函数返回**原始值**,
    log2 变换在 score_one() 里显式做,和上游一致。
    """
    xm, xs = X.mean(axis=0), X.std(axis=0) + epsilon
    ym, ys = Y.mean(axis=0), Y.std(axis=0) + epsilon
    kl = np.log(ys / xs) + (xs ** 2 + (xm - ym) ** 2) / (2 * ys ** 2) - 0.5
    klr = np.log(xs / ys) + (ys ** 2 + (ym - xm) ** 2) / (2 * xs ** 2) - 0.5
    return float(np.mean(kl + klr))


def m_sliced_wasserstein(X: np.ndarray, Y: np.ndarray, n_proj: int = 200,
                         seed: int = SEED) -> float:
    """Sliced Wasserstein-2:随机方向投影后的 1D W2 距离取平均。

    ⚠️ 这**不是**上游/pertpy 的 `wasserstein`(那是 ott-jax Sinkhorn 熵正则 OT 的
    reg_ot_cost)。本机无 ott-jax,故用刻意改名的可本地计算的替代量。
    只用于"同一批方法之间横向排序",不可与论文数值直接对比。
    """
    rng = np.random.default_rng(seed)
    d = X.shape[1]
    V = rng.normal(size=(d, n_proj))
    V /= np.linalg.norm(V, axis=0, keepdims=True)
    px, py = X @ V, Y @ V
    n = min(px.shape[0], py.shape[0])
    q = (np.arange(n) + 0.5) / n                      # 用分位数对齐不等长样本
    out = 0.0
    for j in range(n_proj):
        a = np.quantile(px[:, j], q)
        b = np.quantile(py[:, j], q)
        out += np.mean((a - b) ** 2)
    return float(np.sqrt(out / n_proj))


def _top_deg(mat_a: np.ndarray, mat_b: np.ndarray, genes: list[str], n: int) -> list[str]:
    """按 Welch t 统计量绝对值排序取 top-n 基因(a vs b)。"""
    ma, mb = mat_a.mean(axis=0), mat_b.mean(axis=0)
    va = mat_a.var(axis=0, ddof=1) / max(mat_a.shape[0], 1)
    vb = mat_b.var(axis=0, ddof=1) / max(mat_b.shape[0], 1)
    t = (mb - ma) / np.sqrt(va + vb + 1e-12)
    order = np.argsort(-np.abs(t))[:n]
    return [genes[i] for i in order]


# =============================================================================
# 主流程
# =============================================================================
def load_inputs(observed: Path, predicted: Path):
    obs = pd.read_csv(observed)
    pred = pd.read_csv(predicted)
    for col in ("cell_id", "group"):
        if col not in obs.columns:
            raise SystemExit(f"observed.csv 缺列 `{col}`")
    if "method" not in pred.columns:
        raise SystemExit("predicted.csv 缺列 `method`")
    genes = [c for c in obs.columns if c not in ("cell_id", "group")]
    missing = [g for g in genes if g not in pred.columns]
    if missing:
        raise SystemExit(f"predicted.csv 缺 {len(missing)} 个基因列,例: {missing[:5]}")
    return obs, pred, genes


def split_and_baselines(obs: pd.DataFrame, genes: list[str], seed: int = SEED):
    """把 stimulated 切成 train/eval 两半,并用 train 半构造朴素基线。

    ★ 防数据泄漏:评估用的 ground truth 只用 eval 半;trainMean 与 DEG 选择只看 train 半。
    """
    rng = np.random.default_rng(seed)
    ctrl = obs[obs["group"] == "control"][genes].to_numpy(float)
    stim = obs[obs["group"] == "stimulated"][genes].to_numpy(float)
    idx = rng.permutation(stim.shape[0])
    half = stim.shape[0] // 2
    stim_train, stim_eval = stim[idx[:half]], stim[idx[half:]]

    # ── controlBaseline ──────────────────────────────────────────────────────
    # 对齐上游 baseControl.py(方法名 'controlMean'):
    #   imputed_adata = control_adata.copy()  → 直接把 held-out 的 control 细胞
    #   原样当作预测,不重采样、不改细胞数。即"什么都没发生"。
    control_baseline = ctrl.copy()

    # ── trainMean ────────────────────────────────────────────────────────────
    # 对齐上游 trainMean.py 的 generateExp():
    #   train_mean   = 训练集(其他 cellular context)同一扰动下的 per-gene 均值
    #   control_std  = held-out control 细胞的 per-gene 标准差(NaN→0)
    #   expression   = np.random.normal(loc=train_mean[i], scale=control_std[i],
    #                                   size=cellNum)   , cellNum = control 细胞数
    # ★注意:上游**不是**把均值广播成常数矩阵 —— 它按 control 的方差加了高斯噪声,
    #   所以分布型指标(edistance / sym_kldiv)在上游是有定义的,不会因零方差被特判。
    train_mean = stim_train.mean(axis=0)
    control_std = np.nan_to_num(ctrl.std(axis=0), nan=0.0)
    n_cells = ctrl.shape[0]
    train_mean_pred = rng.normal(loc=train_mean, scale=control_std,
                                 size=(n_cells, train_mean.shape[0]))

    baselines = {"controlBaseline": control_baseline, "trainMean": train_mean_pred}
    return ctrl, stim_train, stim_eval, baselines


def _subsample(X: np.ndarray, n: int, rng: np.random.Generator) -> np.ndarray:
    """对齐上游 f_subSample:每组最多 n 个细胞。"""
    if X.shape[0] <= n:
        return X
    return X[rng.choice(X.shape[0], n, replace=False)]


def score_one(pred_X: np.ndarray, truth_X: np.ndarray, ctrl_X: np.ndarray,
              genes: list[str], metric: str, subsample: int,
              n_deg: int, seed: int = SEED,
              full: tuple | None = None) -> float:
    """给单个 (方法, 指标) 打分。

    pred_X/truth_X/ctrl_X 已经子集到 top-N DEG 基因;`full` 传入未子集的
    (pred, truth, ctrl, genes),只有 common_deg 用它 —— 在已经是 DEG 的基因里
    再选 DEG 会恒等于 1,必须回到全基因空间算。
    """
    rng = np.random.default_rng(seed)
    # 上游 calPerformance.py:108-112 —— 只有分布型指标走 subsample 后的 adata,
    # mse / pearson_distance 用全量细胞。这里照抄该分支条件。
    if metric in SUBSAMPLED_METRICS:
        P = _subsample(pred_X, subsample, rng)
        T = _subsample(truth_X, subsample, rng)
        C = _subsample(ctrl_X, subsample, rng)
    else:
        P, T, C = pred_X, truth_X, ctrl_X

    if metric == "common_deg":
        # 本地实现:在**全基因空间**里,预测 DEG 与真实 DEG 的 top-N 重叠比例(越大越好)。
        # ⚠️ scPerturBench README 列出 "Common-DEGs",但已发布的 calPerformance*.py
        #    中未见其实现代码,故此处口径为本库自定义,不保证与论文数值一致。
        if full is None:
            raise ValueError("common_deg 需要全基因矩阵")
        fP, fT, fC, fg = full
        k = min(n_deg, len(fg))
        true_deg = set(_top_deg(fC, fT, fg, k))
        pred_deg = set(_top_deg(fC, fP, fg, k))
        return len(true_deg & pred_deg) / max(len(true_deg), 1)

    if metric == "pearson_distance":
        # 上游 calculateDelta:imputed / stimulated 各自减去 control 的均值
        cm = C.mean(axis=0)
        return m_pearson_distance(P - cm, T - cm)
    if metric == "mse":
        return m_mse(P, T)
    if metric == "edistance":
        return m_edistance(P, T)
    if metric == "sym_kldiv":
        return float(np.log2(m_sym_kldiv(P, T) + 1))     # 上游同款 log2(x+1)
    if metric == "sliced_wasserstein":
        return m_sliced_wasserstein(P, T, seed=seed)
    raise ValueError(f"未知指标 {metric}")


def score_with_pertpy(preds: dict, truth_X, ctrl_X, genes, metrics, subsample):
    """守卫式引用封装:装了 pertpy 才走,拿到与论文完全同名同实现的指标。

    真实 API(已核对上游源码,见文件头):
        pt.tools.Distance(metric=<m>, layer_key='X')
          .onesided_distances(adata, groupby='perturbation',
                              selected_group='imputed', groups=['stimulated'])
    """
    try:
        import anndata as ad
        import pertpy as pt
    except ImportError as e:
        return {"status": "skipped",
                "reason": f"未安装 {e.name};如需上游同名指标(含 ott-jax Sinkhorn "
                          f"wasserstein)请自行 `pip install pertpy`",
                "note": "本模块默认不依赖 pertpy,本地实现已逐行对齐其公式"}
    out = {"status": "ran", "pertpy_version": getattr(pt, "__version__", "?"), "values": {}}
    for name, P in preds.items():
        for m in metrics:
            # 上游 calPerformance.py:101-102 —— pearson_distance 之前先 calculateDelta:
            # imputed / stimulated 各自减去 control 的均值,control 组保持原样。
            if m == "pearson_distance":
                cm = ctrl_X.mean(axis=0)
                Pm, Tm = P - cm, truth_X - cm
            else:
                Pm, Tm = P, truth_X
            X = np.vstack([ctrl_X, Pm, Tm])
            obs = pd.DataFrame({"perturbation":
                                ["control"] * ctrl_X.shape[0] + ["imputed"] * Pm.shape[0] +
                                ["stimulated"] * Tm.shape[0]})
            A = ad.AnnData(X, obs=obs)
            A.var_names = genes
            A.layers["X"] = A.X
            try:
                D = pt.tools.Distance(metric=m, layer_key="X")
                df = D.onesided_distances(A, groupby="perturbation",
                                          selected_group="imputed", groups=["stimulated"])
                # onesided_distances 返回 DataFrame(index=selected_group);上游写的是
                # pairwise_df['stimulated'],在新版 pandas 里 float(Series) 会告警/报错,
                # 故这里显式取标量。
                val = float(np.asarray(df["stimulated"]).ravel()[0])
                if m == "sym_kldiv":                     # 上游同款 log2(x+1)
                    val = float(np.log2(val + 1))
                out["values"].setdefault(name, {})[m] = val
            except Exception as e:                       # 某些指标需额外依赖(ott-jax 等)
                out["values"].setdefault(name, {})[m] = f"failed: {type(e).__name__}"
    return out


# =============================================================================
# 出图(顶刊风格,禁用条形图)
# =============================================================================
def _rank_frame(perf: pd.DataFrame) -> pd.DataFrame:
    """把每个指标转成排名(1 = 该指标下最好),方向按 LOWER_IS_BETTER 校正。"""
    wide = perf.pivot(index="method", columns="metric", values="performance")
    ranks = pd.DataFrame(index=wide.index, columns=wide.columns, dtype=float)
    for m in wide.columns:
        asc = LOWER_IS_BETTER.get(m, True)
        ranks[m] = wide[m].rank(ascending=asc, na_option="bottom")
    return wide, ranks


def fig_rank_heatmap(ranks: pd.DataFrame, outdir: Path):
    import matplotlib.pyplot as plt
    order = ranks.mean(axis=1).sort_values().index
    R = ranks.loc[order]
    fig, ax = plt.subplots(figsize=(NATURE_W1 + 1.4, 0.42 * len(R) + 1.7))
    im = ax.imshow(R.to_numpy(), cmap=CMAP_CONT + "_r", aspect="auto")
    ax.set_xticks(range(R.shape[1]))
    ax.set_xticklabels(R.columns, rotation=40, ha="right")
    ax.set_yticks(range(R.shape[0]))
    ax.set_yticklabels(R.index)
    for i in range(R.shape[0]):
        for j in range(R.shape[1]):
            v = R.iat[i, j]
            ax.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=8,
                    color="white" if v > R.to_numpy().mean() else "black")
    ax.set_title("Rank per metric (1 = best)")
    fig.colorbar(im, ax=ax, shrink=0.75, label="rank")
    save_fig(fig, outdir / "fig1_rank_heatmap")
    plt.close(fig)


def fig_vs_baseline(wide: pd.DataFrame, baseline: str, outdir: Path):
    """点图:每个指标下各方法相对 controlBaseline 的比值。虚线 1.0 = 基线地板。"""
    import matplotlib.pyplot as plt
    metrics = [m for m in wide.columns if np.isfinite(wide.loc[baseline, m])
               and wide.loc[baseline, m] != 0]
    methods = [m for m in wide.index if m != baseline]
    colors = pal(len(methods), "npg")
    fig, ax = plt.subplots(figsize=(NATURE_W2 * 0.8, 0.55 * len(metrics) + 1.6))
    for yi, met in enumerate(metrics):
        ref = wide.loc[baseline, met]
        ax.axhline(yi, color="#DDDDDD", lw=0.8, zorder=0)
        for k, meth in enumerate(methods):
            v = wide.loc[meth, met]
            if not np.isfinite(v):
                continue
            ax.scatter(v / ref, yi + (k - (len(methods) - 1) / 2) * 0.16,
                       s=46, color=colors[k], zorder=3,
                       label=meth if yi == 0 else None, edgecolor="white", lw=0.6)
    ax.axvline(1.0, color="black", ls="--", lw=1.2, zorder=2)
    ax.set_xscale("log")
    ax.set_yticks(range(len(metrics)))
    ax.set_yticklabels(metrics)
    ax.set_xlabel(f"value relative to {baseline}, log scale\n"
                  "(left of the dashed line = beats the naive floor, for distance metrics)")
    ax.set_title("Does the model beat the naive floor?")
    ax.set_ylim(-0.6, len(metrics) - 0.4)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0))
    save_fig(fig, outdir / "fig2_vs_baseline_dotplot")
    plt.close(fig)


def fig_delta_scatter(preds: dict, truth_X, ctrl_X, outdir: Path):
    """每基因 delta 的预测 vs 真实散点,每方法一个 panel,标注 PCC-delta。"""
    import matplotlib.pyplot as plt
    cm = ctrl_X.mean(axis=0)
    true_d = truth_X.mean(axis=0) - cm
    names = list(preds)
    colors = pal(len(names), "npg")
    fig, axes = plt.subplots(1, len(names), figsize=(2.35 * len(names), 2.6), sharey=True)
    axes = np.atleast_1d(axes)
    lim = None
    for ax, name, c in zip(axes, names, colors):
        pd_ = preds[name].mean(axis=0) - cm
        ax.axhline(0, color="#CCCCCC", lw=0.8)
        ax.axvline(0, color="#CCCCCC", lw=0.8)
        ax.scatter(true_d, pd_, s=16, color=c, alpha=0.8, edgecolor="white", lw=0.4)
        if np.std(pd_) > 0:
            r = pearsonr(true_d, pd_)[0]
            ax.set_title(f"{name}\nPCC-delta = {r:.2f}", fontsize=9)
        else:
            ax.set_title(f"{name}\nPCC-delta = n/a (constant)", fontsize=9)
        if lim is None:
            lim = max(np.abs(true_d).max(), 1e-6) * 1.25
        ax.plot([-lim, lim], [-lim, lim], ls=":", color="black", lw=0.9)
        ax.set_xlim(-lim, lim)
        ax.set_xlabel("observed delta")
    axes[0].set_ylabel("predicted delta")
    save_fig(fig, outdir / "fig3_delta_scatter")
    plt.close(fig)


# =============================================================================
def main():
    ap = argparse.ArgumentParser(description="scPerturBench 风格的扰动预测评测")
    ap.add_argument("--observed", default=str(ROOT / "example_data" / "observed.csv"))
    ap.add_argument("--predicted", default=str(ROOT / "example_data" / "predicted.csv"))
    ap.add_argument("--outdir", default=str(ROOT / "results"))
    ap.add_argument("--n-deg", type=int, default=30,
                    help="先把基因子集到 control-vs-stimulated 的 top-N DEG(上游用 100/5000)")
    ap.add_argument("--subsample", type=int, default=2000, help="每组最多细胞数(上游 2000)")
    ap.add_argument("--metrics", default=",".join(DEFAULT_METRICS))
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--use-pertpy", action="store_true",
                    help="额外用 pertpy 官方 Distance 复算(需自行安装 pertpy)")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    set_pub_style()
    np.random.seed(args.seed)
    metrics = [m.strip() for m in args.metrics.split(",") if m.strip()]

    print("Step 1  读取输入")
    obs, pred, genes = load_inputs(Path(args.observed), Path(args.predicted))
    print(f"        {obs.shape[0]} observed cells · {len(genes)} genes · "
          f"{pred['method'].nunique()} candidate method(s)")

    print("Step 2  切分 train/eval + 构造朴素基线 (controlBaseline / trainMean)")
    ctrl, stim_train, stim_eval, baselines = split_and_baselines(obs, genes, args.seed)
    preds = {m: g[genes].to_numpy(float) for m, g in pred.groupby("method")}
    for k, v in preds.items():
        if v.shape[0] != stim_eval.shape[0]:
            print(f"        note: {k} 有 {v.shape[0]} 个预测细胞,ground truth 有 "
                  f"{stim_eval.shape[0]} 个 —— 分布型指标允许不等样本量")
    preds.update(baselines)
    print(f"        参评方法: {', '.join(preds)}")

    print(f"Step 3  按 top-{args.n_deg} DEG 子集基因(只用 control + stim_train,防泄漏)")
    deg = _top_deg(ctrl, stim_train, genes, min(args.n_deg, len(genes)))
    gi = [genes.index(g) for g in deg]
    ctrl_d, stim_eval_d = ctrl[:, gi], stim_eval[:, gi]
    preds_d = {k: v[:, gi] for k, v in preds.items()}

    print("Step 4  打分")
    rows = []
    for name, P in preds_d.items():
        for met in metrics:
            val = score_one(P, stim_eval_d, ctrl_d, deg, met,
                            args.subsample, args.n_deg, args.seed,
                            full=(preds[name], stim_eval, ctrl, genes))
            rows.append({"method": name, "metric": met, "performance": round(float(val), 4),
                         "n_deg": args.n_deg, "n_pred": P.shape[0],
                         "n_truth": stim_eval_d.shape[0], "seed": args.seed})
            print(f"        {name:<16} {met:<20} {val:.4f}")
    perf = pd.DataFrame(rows)
    perf.to_csv(outdir / "performance.tsv", sep="\t", index=False)

    print("Step 5  排名与出图")
    wide, ranks = _rank_frame(perf)
    wide.to_csv(outdir / "performance_wide.csv")
    ranks.to_csv(outdir / "rank_matrix.csv")
    fig_rank_heatmap(ranks, outdir)
    fig_vs_baseline(wide, "controlBaseline", outdir)
    fig_delta_scatter(preds_d, stim_eval_d, ctrl_d, outdir)

    pertpy_report = {"status": "not requested"}
    if args.use_pertpy:
        print("Step 6  pertpy 官方 Distance 复算")
        pertpy_report = score_with_pertpy(preds_d, stim_eval_d, ctrl_d, deg,
                                          [m for m in metrics if m != "common_deg"],
                                          args.subsample)
        print("        ", pertpy_report.get("status"), pertpy_report.get("reason", ""))

    mean_rank = ranks.mean(axis=1).sort_values()
    beats = {}
    for m in wide.index:
        if m in ("controlBaseline", "trainMean"):
            continue
        wins = []
        for met in wide.columns:
            a, b = wide.loc[m, met], wide.loc["controlBaseline", met]
            if not (np.isfinite(a) and np.isfinite(b)):
                continue
            wins.append(a < b if LOWER_IS_BETTER.get(met, True) else a > b)
        beats[m] = f"{sum(wins)}/{len(wins)} metrics beat controlBaseline"
    summary = {
        "module": "590_scperturbench_generalization",
        "upstream": {"paper": "Nat Methods 2026, PMID 41381899",
                     "doi": "10.1038/s41592-025-02980-0",
                     "repo": "https://github.com/bm2-lab/scPerturBench"},
        "metrics": metrics,
        "n_deg": args.n_deg, "seed": args.seed,
        "mean_rank": {k: round(float(v), 2) for k, v in mean_rank.items()},
        "verdict_vs_naive_floor": beats,
        "pertpy_crosscheck": pertpy_report,
        "session": {                                     # 依赖版本快照(铁律6:可复现)
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": __import__("scipy").__version__,
            "matplotlib": __import__("matplotlib").__version__,
        },
        "caveats": [
            "sliced_wasserstein 不是 pertpy/上游的 Sinkhorn wasserstein,不可与论文数值直接比较",
            "common_deg 为本库自定义口径;上游已发布脚本中未见其实现",
        ],
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False),
                                         encoding="utf-8")
    print("Step 7  完成 ->", outdir)
    for k, v in beats.items():
        print(f"        {k}: {v}")


if __name__ == "__main__":
    main()
