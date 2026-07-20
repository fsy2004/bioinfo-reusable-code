"""591 · scArchon — 单细胞扰动响应预测的基准评测(leave-one-batch-out)。

scArchon 是**基准平台**,不是预测方法:它用 Snakemake + Singularity 把 CellOT / CPA / scGen /
scVIDR / scPRAM / scPreGAN / scDisInFACT / trVAE / SCREEN 以及一个 **linear** 基线包起来,
在同一份数据、同一套指标上做对照(Radig et al., Genome Biology 2026)。

本模块做两件事:
  ① **本机可跑的评测骨架**(默认路径,零依赖外部包):按 scArchon 的实验设计
     (留出一个 batch → 用其余 batch 学扰动响应 → 预测留出 batch 的 control 细胞)
     跑三个朴素预测器,并用与 scArchon 同族的指标打分、出图。
     其中 `linear` 就是 scArchon 仓库里 scripts/linear/snkmk_linear.py 的加性 delta 基线,
     `control` 就是 scripts/metrics/metrics_control.py 的"不做任何预测"地板。
  ② **守卫式引用封装**(--check-scarchon):检查 snakemake/singularity 是否可用,
     打印真实的调用方式与 config/datasets.tsv 列名。scArchon 只有 Snakemake 工作流入口,
     没有可 import 的 Python API,因此这里不伪造任何函数签名。

真实接口核对来源(逐行读过 hdsu-bioquant/scArchon @ main 的本地克隆):
  README.md                           安装/运行命令、batch 不得含空格、60 GB 镜像表
  config/datasets.tsv                 9 个列名(见 check_scarchon)
  Snakefile:93-343                    rule run_{scgen,scdisinfact,scpregan,scpram,screen,
                                      scvidr,cpa,trvae,cellot,linear} → 被封装的 10 个工具
  scripts/linear/snkmk_linear.py:104-118   linear 加性 delta 基线
  scripts/metrics/metrics_control.py:56-83 control 地板(留出 batch 的 ctrl vs stim)
  scripts/metrics/metrics.py:41-190        Metrics 类:pertpy.tl.Distance / common_degs / r2_scores
  论文          PMID 42121287 · doi:10.1186/s13059-026-04104-z · Genome Biol 27(1):162, 2026-05-12
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
import warnings

warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")
ASSETS = os.path.join(HERE, "assets")
EXAMPLE = os.path.join(HERE, "example_data", "synthetic_perturb.h5ad")
FRAMEWORK = os.path.abspath(os.path.join(HERE, "..", "..", "..", "_framework"))
sys.path.insert(0, FRAMEWORK)

SEED = 2026


# ---------------------------------------------------------------- 合成示例数据
def make_example(path: str = EXAMPLE, n_cells: int = 1800, n_genes: int = 300,
                 n_batches: int = 4, seed: int = SEED):
    """合成一个 control/stimulated × 多 batch 的小数据集(synthetic, for demo only)。

    设计要点(为了让基准评测本身有意义):
      · 扰动响应 = 一个 **共享的**基础 delta + 每个 batch 自己的**偏移**
        → 于是"跨 batch 迁移扰动响应"这件事既非平凡也非不可能,基线之间能分出高下;
      · control 状态本身有 batch 效应 → 排除掉"直接背下 stimulated 均值"就能赢的退化情形。
    """
    import anndata as ad
    import numpy as np
    import pandas as pd

    rng = np.random.default_rng(seed)
    genes = [f"G{i:03d}" for i in range(n_genes)]

    base = rng.gamma(2.0, 1.2, n_genes)                       # 基础表达
    shared_delta = np.zeros(n_genes)                          # 共享扰动响应(可迁移部分)
    resp_idx = rng.choice(n_genes, size=40, replace=False)    # 真正响应的基因
    shared_delta[resp_idx] = rng.normal(0, 1.4, 40)

    X, obs = [], []
    per = n_cells // (n_batches * 2)
    for b in range(n_batches):
        batch_shift = rng.normal(0, 0.35, n_genes)            # control 侧 batch 效应
        batch_delta = shared_delta + rng.normal(0, 0.45, n_genes)  # batch 特异的响应偏移
        for cond, extra in (("control", 0.0), ("stimulated", 1.0)):
            mu = base + batch_shift + extra * batch_delta
            mat = rng.normal(mu, 0.55, size=(per, n_genes))
            X.append(mat)
            obs.append(pd.DataFrame({
                "condition": cond,
                "batch": f"donor{b+1}",
            }, index=[f"{cond[:4]}_{b}_{i}" for i in range(per)]))

    A = ad.AnnData(np.vstack(X).astype("float32"))
    A.obs = pd.concat(obs)
    A.var_names = genes
    A.uns["synthetic"] = "synthetic, for demo only — 591 scArchon module"
    A.uns["true_response_genes"] = [genes[i] for i in sorted(resp_idx)]
    os.makedirs(os.path.dirname(path), exist_ok=True)
    A.write_h5ad(path)
    return path


# ---------------------------------------------------------------- 预测器(基线)
def predict_control(ctrl_test, ctrl_train, stim_train):
    """地板 0:不做任何预测,直接把 control 当成预测值。

    对应 scripts/metrics/metrics_control.py:56-59 —— 上游取留出 batch(target)的
    control 与 stimulated 两组,直接算二者的距离,即"完全不预测"时的分数。
    任何方法若打不过它,说明它没有学到任何扰动信息。
    """
    return ctrl_test.copy()


def predict_linear(ctrl_test, ctrl_train, stim_train):
    """scArchon 自带的 linear 基线:加性 delta。

    逐行对应 scripts/linear/snkmk_linear.py:99-118 ——
      train = adata[~(adata.obs[batch] == target)]                       # L99
      gene_diff_train = stimulated_avg_train - control_avg_train         # L112
      predicted = ctrl_adata.copy(); predicted.X += gene_diff_train      # L117-118
    即 delta = mean(stim_train) - mean(ctrl_train);pred = control_test + delta。
    """
    delta = stim_train.mean(axis=0) - ctrl_train.mean(axis=0)
    return ctrl_test + delta[None, :]


def predict_global_mean(ctrl_test, ctrl_train, stim_train):
    """地板 2:忽略细胞自身状态,一律预测训练集的 stimulated 平均谱。

    用来检验指标是否被"只要均值对就得分"的退化解骗到(细胞级别的异质性完全丢失)。
    """
    import numpy as np
    return np.repeat(stim_train.mean(axis=0)[None, :], ctrl_test.shape[0], axis=0)


PREDICTORS = {
    "control": predict_control,
    "linear": predict_linear,
    "global_mean": predict_global_mean,
}


# ---------------------------------------------------------------- 指标
def _rbf_mmd(A, B, n_sub=250, seed=SEED):
    """MMD²(RBF 核,median heuristic 定带宽)。子采样以控制 O(n²) 开销。"""
    import numpy as np
    from sklearn.metrics.pairwise import rbf_kernel
    from scipy.spatial.distance import pdist

    rng = np.random.default_rng(seed)
    a = A[rng.choice(A.shape[0], min(n_sub, A.shape[0]), replace=False)]
    b = B[rng.choice(B.shape[0], min(n_sub, B.shape[0]), replace=False)]
    med = np.median(pdist(np.vstack([a, b])[: 2 * n_sub]))
    gamma = 1.0 / (2 * med ** 2 + 1e-12)
    kaa, kbb, kab = rbf_kernel(a, a, gamma), rbf_kernel(b, b, gamma), rbf_kernel(a, b, gamma)
    return float(kaa.mean() + kbb.mean() - 2 * kab.mean())


def score(pred, true_stim, ctrl_test, n_deg=(20, 100)):
    """与 scArchon scripts/metrics/metrics.py 同族、但**并非同一公式**的指标。

    上游实际做法(metrics.py):
      · 距离项 pertpy.tl.Distance(m).bootstrap(...)  (L73-79),m ∈ {mse, wasserstein,
        pearson_distance, mmd, t_test, cosine_distance};
      · common_degs (L83-103):sc.tl.rank_genes_groups(method="t-test", reference=control)
        后取 stimulated / predicted 两组排名的 top-20 / top-100 名单求交;
      · r2_scores (L107-190):对 top-N DEG 的**均值表达向量**算 stim 与 pred 的
        np.corrcoef(...)**2,bootstrap 重采样 80% 取均值。

    本模块的差异(本机未装 pertpy / 不跑 rank_genes_groups,故为自建实现):
      · 距离项用 numpy/scipy/sklearn 重写同名同义的量,**数值不与 pertpy 逐位一致**,
        且 pertpy 各指标的内部定义未逐条核对;只用于同一批方法之间的相对排序;
      · DEG 名单按 |delta| = |mean(stim) − mean(ctrl)| 排序取 top-N,不做统计检验;
      · r2_delta_* 是 **delta 空间**的决定系数 1 − SS_res/SS_tot,比上游"均值向量相关
        系数平方"更严格 —— 后者会因 control 与 stimulated 的基线表达高度相关而恒接近 1,
        分不出方法高下,故此处刻意换成 delta 空间。**因此不可与论文中的 r2_* 数值对比。**
    """
    import numpy as np
    from scipy.stats import pearsonr, wasserstein_distance

    m_pred, m_true, m_ctrl = pred.mean(0), true_stim.mean(0), ctrl_test.mean(0)
    out = {}

    # —— 分布层面的距离(越小越好)
    out["mse"] = float(np.mean((m_pred - m_true) ** 2))
    out["wasserstein_1d_mean"] = float(np.mean([
        wasserstein_distance(pred[:, g], true_stim[:, g]) for g in range(pred.shape[1])
    ]))
    out["mmd_rbf"] = _rbf_mmd(pred, true_stim)
    cs = float(np.dot(m_pred, m_true) / (np.linalg.norm(m_pred) * np.linalg.norm(m_true) + 1e-12))
    out["cosine_distance"] = 1.0 - cs
    out["pearson_distance"] = 1.0 - float(pearsonr(m_pred, m_true)[0])

    # —— DEG 层面:真扰动方向上有没有学对(delta 空间才是难的部分)
    d_true, d_pred = m_true - m_ctrl, m_pred - m_ctrl
    rank_true = np.argsort(-np.abs(d_true))
    for k in n_deg:
        idx = rank_true[:k]
        ss_res = np.sum((d_pred[idx] - d_true[idx]) ** 2)
        ss_tot = np.sum((d_true[idx] - d_true[idx].mean()) ** 2)
        out[f"r2_delta_top{k}"] = float(1 - ss_res / (ss_tot + 1e-12))
        rank_pred = np.argsort(-np.abs(d_pred))[:k]
        out[f"common_degs_top{k}"] = int(len(set(idx.tolist()) & set(rank_pred.tolist())))
    out["r2_delta_all"] = float(1 - np.sum((d_pred - d_true) ** 2) /
                                (np.sum((d_true - d_true.mean()) ** 2) + 1e-12))
    return out


# 指标方向:True = 越大越好
HIGHER_BETTER = {
    "mse": False, "wasserstein_1d_mean": False, "mmd_rbf": False,
    "cosine_distance": False, "pearson_distance": False,
    "r2_delta_top20": True, "r2_delta_top100": True, "r2_delta_all": True,
    "common_degs_top20": True, "common_degs_top100": True,
}


# ---------------------------------------------------------------- 评测主循环
def leave_one_batch_out(adata, cond_key="condition", ctrl="control",
                        stim="stimulated", batch_key="batch"):
    """scArchon 的实验设计:每次留出一个 batch(datasets.tsv 里的 `target`)。"""
    import numpy as np
    import pandas as pd

    Xd = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    obs = adata.obs
    rows, per_gene = [], []
    for held in sorted(obs[batch_key].unique()):
        te = obs[batch_key] == held
        ctrl_test = Xd[(te & (obs[cond_key] == ctrl)).values]
        stim_test = Xd[(te & (obs[cond_key] == stim)).values]
        ctrl_train = Xd[(~te & (obs[cond_key] == ctrl)).values]
        stim_train = Xd[(~te & (obs[cond_key] == stim)).values]
        if min(map(len, (ctrl_test, stim_test, ctrl_train, stim_train))) < 5:
            print(f"       [skip] batch {held}: 某一格细胞数 < 5")
            continue
        for name, fn in PREDICTORS.items():
            pred = fn(ctrl_test, ctrl_train, stim_train)
            s = score(pred, stim_test, ctrl_test)
            s.update(method=name, held_out_batch=held)
            rows.append(s)
            d_true = stim_test.mean(0) - ctrl_test.mean(0)
            d_pred = pred.mean(0) - ctrl_test.mean(0)
            per_gene.append(pd.DataFrame({
                "method": name, "held_out_batch": held,
                "gene": adata.var_names, "delta_true": d_true, "delta_pred": d_pred,
                "abs_err": np.abs(d_pred - d_true),
            }))
    return pd.DataFrame(rows), pd.concat(per_gene, ignore_index=True)


# ---------------------------------------------------------------- 出图(禁用条形图)
def make_figures(scores, per_gene, outdir):
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    from pubstyle import CMAP_DIVERGE, pal, save_fig, set_pub_style

    set_pub_style()
    cols = pal(len(PREDICTORS))
    cmap_m = dict(zip(PREDICTORS, cols))
    metric_cols = [c for c in scores.columns if c in HIGHER_BETTER]
    figs = {}

    # ① dumbbell:每个留出 batch,control 地板 → linear 的 R²(top100) 位移
    fig, ax = plt.subplots(figsize=(5.6, 3.4))
    piv = scores.pivot(index="held_out_batch", columns="method", values="r2_delta_top100")
    y = np.arange(len(piv))
    for i, (b, r) in enumerate(piv.iterrows()):
        ax.plot([r["control"], r["linear"]], [i, i], color="0.65", lw=2, zorder=1)
    for m in ("control", "linear", "global_mean"):
        ax.scatter(piv[m], y, s=70, color=cmap_m[m], label=m, zorder=3,
                   edgecolor="white", linewidth=0.8)
    ax.axvline(0, ls="--", lw=1, color="0.4")
    ax.set_yticks(y); ax.set_yticklabels(piv.index)
    ax.set_xlabel(r"$R^2$ of predicted perturbation delta (top-100 DEGs)")
    ax.set_ylabel("Held-out batch")
    ax.set_title("Leave-one-batch-out: does the method beat the control floor?")
    ax.set_ylim(-0.6, len(piv) - 0.4)
    ax.legend(frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    save_fig(fig, os.path.join(outdir, "591_dumbbell_r2_by_batch")); plt.close(fig)
    figs["591_dumbbell_r2_by_batch.png"] = "dumbbell · 各留出 batch 上 delta-R² 相对 control 地板的位移"

    # ② heatmap:方法 × 指标,按指标方向统一成 0-1(1 = 该指标下最好)
    agg = scores.groupby("method")[metric_cols].mean()
    norm = agg.copy()
    for c in metric_cols:
        v = agg[c].values.astype(float)
        rng_ = v.max() - v.min()
        z = (v - v.min()) / rng_ if rng_ > 1e-12 else np.zeros_like(v)
        norm[c] = z if HIGHER_BETTER[c] else 1 - z
    order = [m for m in ["control", "linear", "global_mean"] if m in norm.index]
    # agg 与 norm 必须同序,否则格内标注的原始值会错配到别的方法行上
    norm, agg = norm.loc[order], agg.loc[order]
    fig, ax = plt.subplots(figsize=(7.6, 2.6))
    im = ax.imshow(norm.values, cmap="viridis", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(metric_cols)))
    ax.set_xticklabels(metric_cols, rotation=40, ha="right", fontsize=8)
    ax.set_yticks(range(len(norm))); ax.set_yticklabels(norm.index)
    for i in range(norm.shape[0]):
        for j in range(norm.shape[1]):
            ax.text(j, i, f"{agg.iloc[i, j]:.2f}", ha="center", va="center",
                    fontsize=7, color="white" if norm.values[i, j] < 0.55 else "black")
    fig.colorbar(im, ax=ax, shrink=0.8, label="normalised (1 = best)")
    ax.set_title("Method × metric (cell values = raw metric, colour = normalised)")
    save_fig(fig, os.path.join(outdir, "591_heatmap_method_metric")); plt.close(fig)
    figs["591_heatmap_method_metric.png"] = "heatmap · 方法 × 指标(格内为原始值,配色为归一化后的好坏)"

    # ③ 散点:真 delta vs 预测 delta(逐基因),看方法是把响应学对了还是压平了
    fig, axes = plt.subplots(1, len(PREDICTORS), figsize=(10.5, 3.4), sharex=True, sharey=True)
    for ax, m in zip(np.ravel(axes), PREDICTORS):
        d = per_gene[per_gene.method == m]
        ax.scatter(d.delta_true, d.delta_pred, s=6, alpha=0.35, color=cmap_m[m],
                   edgecolor="none")
        lim = [per_gene.delta_true.min() * 1.1, per_gene.delta_true.max() * 1.1]
        ax.plot(lim, lim, ls="--", lw=1, color="0.4")
        # control 的预测 delta 恒为 0(零方差)→ 相关系数无定义,如实标注而不是打印 nan
        if d.delta_pred.std() < 1e-12:
            ax.set_title(f"{m}  (r undefined: constant prediction)", fontsize=9)
        else:
            r = np.corrcoef(d.delta_true, d.delta_pred)[0, 1]
            ax.set_title(f"{m}  (r = {r:.2f})", fontsize=10)
        ax.set_xlabel("True delta (stim − ctrl)")
    np.ravel(axes)[0].set_ylabel("Predicted delta")
    save_fig(fig, os.path.join(outdir, "591_scatter_delta_recovery")); plt.close(fig)
    figs["591_scatter_delta_recovery.png"] = "散点 · 逐基因真实 delta vs 预测 delta(对角线=完美)"

    # ④ raincloud:逐基因绝对误差分布(violin + 抖动散点 + 箱)
    fig, ax = plt.subplots(figsize=(5.6, 3.6))
    rng = np.random.default_rng(SEED)
    data = [per_gene.loc[per_gene.method == m, "abs_err"].values for m in PREDICTORS]
    vp = ax.violinplot(data, positions=np.arange(len(data)), showextrema=False, widths=0.75)
    for body, c in zip(vp["bodies"], cols):
        body.set_facecolor(c); body.set_alpha(0.35); body.set_edgecolor(c)
    for i, (d, c) in enumerate(zip(data, cols)):
        sub = d[rng.choice(len(d), min(400, len(d)), replace=False)]
        ax.scatter(np.full(len(sub), i) + rng.normal(0, 0.045, len(sub)) - 0.25, sub,
                   s=4, alpha=0.3, color=c, edgecolor="none")
    bp = ax.boxplot(data, positions=np.arange(len(data)) + 0.22, widths=0.12,
                    showfliers=False, patch_artist=True)
    for patch in bp["boxes"]:
        patch.set_facecolor("white")
    ax.set_xticks(range(len(data))); ax.set_xticklabels(list(PREDICTORS))
    ax.set_ylabel("Per-gene |predicted − true| delta")
    ax.set_title("Error distribution, pooled over held-out batches")
    save_fig(fig, os.path.join(outdir, "591_raincloud_gene_error")); plt.close(fig)
    figs["591_raincloud_gene_error.png"] = "raincloud · 逐基因绝对误差分布(violin+散点+箱)"

    for f in figs:
        shutil.copy(os.path.join(outdir, f), os.path.join(ASSETS, f))
    return figs


# ---------------------------------------------------------------- 守卫式:真 scArchon
def check_scarchon():
    """检查能否跑真正的 scArchon 工作流。不伪造 API:它只有 Snakemake 入口。"""
    info = {
        "snakemake": shutil.which("snakemake"),
        "singularity": shutil.which("singularity") or shutil.which("apptainer"),
    }
    info["ready"] = bool(info["snakemake"] and info["singularity"])
    info["repo"] = "https://github.com/hdsu-bioquant/scArchon"
    info["install"] = ("conda create -c conda-forge -c bioconda -n snakemake_env snakemake && "
                       "conda activate snakemake_env   # 另需 CUDA>=12.4, Singularity>=3.6, ~60GB 磁盘")
    info["run"] = ("snakemake --use-singularity --singularity-args '--nv -B .:/dum' "
                   "--cores all --jobs 1 --keep-going")
    info["config_tsv_columns"] = ["file_path", "condition", "condition_control",
                                  "condition_stimulated", "batch", "target",
                                  "experiment_name", "output_dir", "Tools"]
    info["wrapped_tools"] = ["cellot", "cpa", "scgen", "scvidr", "scpram", "scpregan",
                             "scdisinfact", "trvae", "screen", "linear"]
    info["note"] = ("scArchon 无可 import 的 Python API,只能通过 Snakemake 驱动;"
                    "工具名小写、batch 取值不得含空格(README 明确要求)。")
    return info


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--h5ad", default=EXAMPLE, help="AnnData:含 condition/batch 两列 obs")
    ap.add_argument("--condition-key", default="condition")
    ap.add_argument("--control", default="control")
    ap.add_argument("--stimulated", default="stimulated")
    ap.add_argument("--batch-key", default="batch")
    ap.add_argument("--outdir", default=RESULTS)
    ap.add_argument("--make-example", action="store_true", help="重新生成合成示例数据后退出")
    ap.add_argument("--check-scarchon", action="store_true", help="只检查真 scArchon 工作流是否可跑")
    a = ap.parse_args()

    if a.make_example:
        print("[591] 生成合成示例数据 ->", make_example()); return

    if a.check_scarchon:
        print("[591] scArchon 工作流可用性检查")
        for k, v in check_scarchon().items():
            print(f"       {k}: {v}")
        return

    os.makedirs(a.outdir, exist_ok=True); os.makedirs(ASSETS, exist_ok=True)
    import anndata as ad

    if not os.path.exists(a.h5ad):
        if a.h5ad == EXAMPLE:
            print("[591] 示例数据缺失,先生成"); make_example()
        else:
            sys.exit(f"找不到输入文件:{a.h5ad}")
    adata = ad.read_h5ad(a.h5ad)
    for k in (a.condition_key, a.batch_key):
        if k not in adata.obs:
            sys.exit(f"obs 缺少列 '{k}' —— scArchon 的实验设计需要 condition 与 batch 两列")
    print(f"[591] 输入 {adata.n_obs} cells × {adata.n_vars} genes · "
          f"batches={sorted(adata.obs[a.batch_key].unique())}")

    print("[591] leave-one-batch-out 评测(control / linear / global_mean)")
    scores, per_gene = leave_one_batch_out(adata, a.condition_key, a.control,
                                           a.stimulated, a.batch_key)
    if scores.empty:
        sys.exit("没有任何 batch 通过最小细胞数检查")
    scores.to_csv(os.path.join(a.outdir, "591_scores_per_batch.csv"), index=False)
    per_gene.to_csv(os.path.join(a.outdir, "591_per_gene_delta.csv"), index=False)

    metric_cols = [c for c in scores.columns if c in HIGHER_BETTER]
    agg = scores.groupby("method")[metric_cols].mean()
    agg.to_csv(os.path.join(a.outdir, "591_scores_mean.csv"))
    print(agg[["mse", "pearson_distance", "r2_delta_top100", "common_degs_top100"]]
          .round(3).to_string())

    print("[591] 出图")
    figs = make_figures(scores, per_gene, a.outdir)

    sc = check_scarchon()
    # 依赖版本快照(铁律6:可复现)
    import importlib
    versions = {"python": sys.version.split()[0]}
    for p in ("numpy", "scipy", "pandas", "sklearn", "anndata", "matplotlib"):
        try:
            versions[p] = getattr(importlib.import_module(p), "__version__", "?")
        except ImportError:
            versions[p] = "not installed"
    with open(os.path.join(a.outdir, "591_summary.json"), "w", encoding="utf-8") as fh:
        json.dump({"input": os.path.basename(a.h5ad),
                   "n_cells": int(adata.n_obs), "n_genes": int(adata.n_vars),
                   "batches": sorted(map(str, adata.obs[a.batch_key].unique())),
                   "mean_scores": agg.round(4).to_dict(),
                   "figures": list(figs),
                   "seed": SEED,
                   "versions": versions,
                   "scarchon_workflow_ready": sc["ready"]},
                  fh, indent=1, ensure_ascii=False, default=str)
    print(f"[591] 结果写入 {a.outdir};展示图复制到 assets/")
    if not sc["ready"]:
        print("[591] 真 scArchon 工作流不可用(缺 snakemake/singularity);"
              "本次仅为本机基线评测。安装见 --check-scarchon")


if __name__ == "__main__":
    main()
