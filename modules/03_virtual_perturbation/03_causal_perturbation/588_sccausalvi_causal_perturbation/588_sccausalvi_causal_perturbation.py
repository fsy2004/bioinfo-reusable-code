"""588 · scCausalVI — 因果解耦的扰动响应建模(背景状态 vs 处理效应)+ 可跑线性基线.

scCausalVI(An et al., *Cell Systems* 2025;16(11):101443, PMID 41197632,
doi:10.1016/j.cels.2025.101443)用「深层结构因果网络」把单细胞表达拆成两组隐变量:
  · background latent —— 细胞固有状态(去掉处理影响)
  · treatment-effect latent —— 处理诱导的、细胞状态特异的转录改变
并支持 cross-condition 反事实预测(把 stim 细胞映射回 ctrl,或反之)与响应细胞识别。

本模块两条路:
  A. **线性基线(默认,本机依赖即可跑通)** —— 用条件均值中心化 + PCA 做"背景表示",
     用全局/细胞类型特异 Δ 做反事实预测,用共享空间 kNN 距离做响应细胞打分。
     这是本库的规矩:任何深度扰动模型都必须先跟朴素线性对照比(DL 扰动模型输给线性
     基线是已发表的常见结果),没有基线就不许单报深度模型。
  B. **scCausalVI 真实路径(--run-sccausalvi,需 pip install scCausalVI)** ——
     签名全部逐行核对自上游源码(见 README「API 来源」的 文件:行号 表),未臆造;
     但本机未安装该包,**这条路径未在本机实跑验证**。

图中文字英文,代码注释中文。固定随机种子。
"""
from __future__ import annotations

import argparse
import json
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
    from pubstyle import set_pub_style, save_fig, pal, NATURE_W2, CMAP_CONT
except Exception:                                     # 框架缺失时的最小降级(不影响分析)
    def set_pub_style(*a, **k): pass

    def save_fig(fig, f, dpi=300):
        Path(f).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(str(f) + ".pdf"); fig.savefig(str(f) + ".png", dpi=dpi)

    def pal(n=None, name="npg"):
        base = ["#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4"]
        return base[:n] if n else base
    NATURE_W2 = 7.0
    CMAP_CONT = "viridis"

SEED = 42
EXAMPLE = HERE / "example_data" / "synthetic_counts.npz"
RESULTS = HERE / "results"
ASSETS = HERE / "assets"


# =============================================================================
# 0. 合成示例数据(synthetic, for demo only)
# =============================================================================
def make_example(path: Path = EXAMPLE, n_per_group: int = 160, n_genes: int = 300,
                 seed: int = SEED) -> Path:
    """造一个"处理效应是细胞类型特异的"case-control 数据集。

    这正是 scCausalVI 声称要解决的场景:全局差异表达会把细胞类型特异的响应抹平。
    真值:CT_A 强响应、CT_B 中等响应、CT_C 不响应(阴性对照细胞类型)。
    另加一个与 condition 平衡(不混杂)的 batch 效应,用来检验背景表示是否被技术变异污染。
    """
    rng = np.random.default_rng(seed)
    cts = ["CT_A", "CT_B", "CT_C"]
    conds = ["ctrl", "stim"]

    # 细胞类型固有表达程序(背景因子的真值)
    ct_mean = {c: rng.gamma(2.0, 1.2, n_genes) + 0.2 for c in cts}

    # 细胞类型特异的处理效应:各自作用在不同基因块上,强度不同
    te = {c: np.zeros(n_genes) for c in cts}
    te["CT_A"][0:40] = rng.uniform(0.9, 1.6, 40)          # 强响应
    te["CT_B"][25:60] = rng.uniform(0.35, 0.7, 35)        # 中等响应,基因块与 A 部分重叠
    # CT_C 全 0 —— 阴性对照

    # batch 效应:乘性,基因维度;两个 batch 在 condition 上均衡分配(不与处理混杂)
    batch_f = {0: np.ones(n_genes), 1: np.exp(rng.normal(0, 0.30, n_genes))}

    X, ct_lab, cond_lab, bat_lab = [], [], [], []
    for ct in cts:
        for cd in conds:
            eff = te[ct] if cd == "stim" else np.zeros(n_genes)
            for i in range(n_per_group):
                b = i % 2                                  # batch 与 condition 正交
                lam = ct_mean[ct] * np.exp(eff) * batch_f[b]
                lam = lam * rng.uniform(0.7, 1.4)          # 细胞文库大小波动
                X.append(rng.poisson(lam))
                ct_lab.append(ct); cond_lab.append(cd); bat_lab.append(b)

    X = np.asarray(X, dtype=np.float32)
    ct_lab = np.asarray(ct_lab); cond_lab = np.asarray(cond_lab)
    bat_lab = np.asarray(bat_lab, dtype=np.int8)
    # 真值:该细胞是否属于"会响应"的细胞类型(仅用于评估,不参与打分)
    responder = np.isin(ct_lab, ["CT_A", "CT_B"]) & (cond_lab == "stim")

    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path, X=X, cell_type=ct_lab, condition=cond_lab, batch=bat_lab,
        true_responder=responder,
        gene_names=np.asarray([f"G{i:03d}" for i in range(n_genes)]),
        note=np.asarray(["synthetic, for demo only -- module 588 scCausalVI"]),
    )
    return path


def load_data(path: Path):
    if not path.exists():
        print(f"[588] example data missing, generating -> {path}")
        make_example(path)
    d = np.load(path, allow_pickle=True)
    return (d["X"].astype(np.float64), d["cell_type"].astype(str),
            d["condition"].astype(str), d["batch"].astype(int),
            d["true_responder"].astype(bool), d["gene_names"].astype(str))


# =============================================================================
# 1. 预处理(CPM + log1p)—— 仅用于本模块的线性基线
# 注:上游 ifnb tutorial 用的是 target_sum=1e6,上游 responsive_cells() 的默认值是 1e4;
#     基线只在自身内部横向比较,故这里取 1e4,不声称与 tutorial 口径一致。
# =============================================================================
def lognorm(X: np.ndarray, target_sum: float = 1e4) -> np.ndarray:
    lib = X.sum(1, keepdims=True)
    lib[lib == 0] = 1.0
    return np.log1p(X / lib * target_sum)


# =============================================================================
# 2. 基线 A:背景表示 —— 原始 PCA vs 条件中心化 PCA
# =============================================================================
def background_representations(L: np.ndarray, cond: np.ndarray, n_pcs: int = 20,
                               seed: int = SEED):
    """两种"背景隐空间"的朴素做法:

    raw     —— 直接 PCA(不做任何解耦),处理效应必然泄漏进去
    centered—— 每个 condition 各自减去自己的基因均值再 PCA(线性去处理效应)
               这是 scCausalVI 的 background latent 的线性对照物
    """
    from sklearn.decomposition import PCA
    reps = {}
    reps["raw"] = PCA(n_components=n_pcs, random_state=seed).fit_transform(L)

    Lc = L.copy()
    for c in np.unique(cond):                                  # 逐条件去均值
        m = cond == c
        Lc[m] -= Lc[m].mean(0, keepdims=True)
    reps["condition_centered"] = PCA(n_components=n_pcs, random_state=seed).fit_transform(Lc)
    return reps


def disentanglement_metrics(rep: np.ndarray, cond: np.ndarray, ct: np.ndarray,
                            seed: int = SEED) -> dict:
    """解耦好不好 = ①条件信息泄漏低 ②细胞类型信息保留高。

    两个指标都用交叉验证外样本预测,避免"用同一批细胞既拟合又评估"的循环分析。
    condition_leakage_auc: 0.5 = 背景表示里读不出处理条件(理想)
    celltype_knn_acc     : 越高 = 固有细胞状态保留越好
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.neighbors import KNeighborsClassifier
    from sklearn.model_selection import cross_val_predict, StratifiedKFold
    from sklearn.preprocessing import StandardScaler
    from sklearn.pipeline import make_pipeline
    from sklearn.metrics import roc_auc_score, accuracy_score

    cv = StratifiedKFold(5, shuffle=True, random_state=seed)
    y = (cond == "stim").astype(int)
    clf = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000))
    p = cross_val_predict(clf, rep, y, cv=cv, method="predict_proba")[:, 1]
    # 取 max(auc, 1-auc):方向无所谓,能被解码出来就是泄漏(<0.5 同样是可解码信号)
    auc = max(roc_auc_score(y, p), 1 - roc_auc_score(y, p))

    knn = make_pipeline(StandardScaler(), KNeighborsClassifier(15))
    pred = cross_val_predict(knn, rep, ct, cv=cv)
    return {"condition_leakage_auc": float(auc),
            "celltype_knn_acc": float(accuracy_score(ct, pred))}


# =============================================================================
# 3. 基线 B:反事实预测(ctrl -> stim)的线性对照
# =============================================================================
def counterfactual_baselines(L: np.ndarray, cond: np.ndarray, ct: np.ndarray,
                             seed: int = SEED):
    """两个线性反事实模型,在**留出细胞**上评估,杜绝数据泄漏。

    global   : Δ 用全体细胞的 stim-ctrl 基因均值差(不看细胞类型)
    celltype : Δ 按细胞类型各算一份(更强的线性基线)
    评估:在留出集上比较"预测的每类每基因 Δ"与"观测到的 Δ"(Pearson r / R²)。
    """
    rng = np.random.default_rng(seed)
    n = L.shape[0]
    fit = np.zeros(n, bool)
    # 按 celltype × condition 分层各取一半做拟合,另一半做评估
    for c in np.unique(ct):
        for d in np.unique(cond):
            idx = np.where((ct == c) & (cond == d))[0]
            fit[rng.permutation(idx)[: len(idx) // 2]] = True
    ev = ~fit

    def gmean(mask):
        return L[mask].mean(0)

    d_global = gmean(fit & (cond == "stim")) - gmean(fit & (cond == "ctrl"))
    d_ct = {c: gmean(fit & (ct == c) & (cond == "stim")) - gmean(fit & (ct == c) & (cond == "ctrl"))
            for c in np.unique(ct)}

    rows, preds = [], {"global": [], "celltype": []}
    obs = []
    for c in np.unique(ct):
        o = gmean(ev & (ct == c) & (cond == "stim")) - gmean(ev & (ct == c) & (cond == "ctrl"))
        obs.append(o); preds["global"].append(d_global); preds["celltype"].append(d_ct[c])
        rows.append(c)
    obs = np.concatenate(obs)
    out = {}
    for k, v in preds.items():
        pv = np.concatenate(v)
        r = float(np.corrcoef(pv, obs)[0, 1])
        ss = 1 - float(((obs - pv) ** 2).sum() / ((obs - obs.mean()) ** 2).sum())
        out[k] = {"pearson_r": r, "r2": ss, "pred": pv}
    return out, obs, rows, d_ct


# =============================================================================
# 4. 基线 C:响应细胞打分(scCausalVI responsive_cells 的非参数对照)
# =============================================================================
def responsive_score_baseline(L: np.ndarray, cond: np.ndarray, n_pcs: int = 20,
                              k: int = 15, seed: int = SEED):
    """每个 stim 细胞到 control 流形的 kNN 距离,用 control-control 的同类距离作零分布标准化。

    对应 scCausalVI 的思路:处理引起的偏移要跟"生成模型自身的重构不确定性"比;
    这里用最朴素的可跑替代(不需要任何标签,因此评估不循环)。
    """
    from sklearn.decomposition import PCA
    from sklearn.neighbors import NearestNeighbors
    Z = PCA(n_components=n_pcs, random_state=seed).fit_transform(L)
    ctrl = Z[cond == "ctrl"]
    nn = NearestNeighbors(n_neighbors=k + 1).fit(ctrl)

    # 零分布:control 细胞到 control 流形的距离(去掉自身)
    d_null = nn.kneighbors(ctrl)[0][:, 1:].mean(1)
    mu, sd = d_null.mean(), d_null.std() + 1e-9
    d_all = nn.kneighbors(Z)[0]
    # 对 control 细胞要去掉"自己"这一近邻,对 stim 细胞不需要
    d_cell = np.where(cond == "ctrl", d_all[:, 1:].mean(1), d_all[:, :k].mean(1))
    return (d_cell - mu) / sd


# =============================================================================
# 5. scCausalVI 真实路径(守卫式;签名取自上游源码/tutorial,本机未安装故未实跑)
# =============================================================================
def run_sccausalvi(X, cond, ct, control: str = "ctrl", max_epochs: int = 200,
                   outdir: Path = RESULTS) -> dict:
    """按官方 tutorial 的顺序调用真实 API(签名逐行核对自上游源码,行号见 README)。

    API 来源(上游 commit 本地克隆,文件:行号):
      · scCausalVI/__init__.py:3,16                     —— 包导出 scCausalVIModel
      · model/scCausalVI.py:54                          —— scCausalVIModel.__init__
      · model/scCausalVI.py:116                         —— setup_anndata (classmethod)
      · model/scCausalVI.py:162                         —— get_latent_representation -> (bg, te)
      · model/scCausalVI.py:683                         —— responsive_cells(仅源码,tutorial 未演示)
      · model/base/training_mixin.py:9                  —— train(group_indices_list, ...)
      · docs/source/tutorial/scCausalVI-ifnb.ipynb      —— condition2int 构造与调用顺序
    包未安装时优雅退出并打印真实安装命令,绝不伪造结果。
    """
    try:
        import anndata as ad
        import pandas as pd
        import torch
        from scCausalVI import scCausalVIModel
    except ImportError as e:
        return {"status": "skipped",
                "reason": f"{getattr(e, 'name', e)} not installed",
                "install": "pip install scCausalVI   # 或 pip install git+https://github.com/ShaokunAn/scCausalVI.git"}

    A = ad.AnnData(np.asarray(X, dtype=np.float32))
    A.obs["condition"] = pd.Categorical(cond)
    A.obs["cell_type"] = pd.Categorical(ct)
    A.layers["counts"] = A.X.copy()

    # 条件顺序固定为「对照在前」,与 tutorial 的 conditions=['ctrl','stim'] 一致(可复现)
    conditions = [control] + sorted(c for c in set(map(str, cond)) if c != control)
    if control not in set(map(str, cond)):
        return {"status": "error",
                "reason": f"control condition '{control}' not found in condition labels "
                          f"{sorted(set(map(str, cond)))}"}
    group_indices_list = [np.where(cond == g)[0] for g in conditions]

    scCausalVIModel.setup_anndata(A, condition_key="condition", layer="counts")
    condition2int = (A.obs.groupby("condition", observed=False)["_scvi_condition"]
                     .first().to_dict())

    model = scCausalVIModel(A, condition2int=condition2int, control=control,
                            n_background_latent=10, n_te_latent=10,
                            n_layers=2, n_hidden=128,
                            use_mmd=True, mmd_weight=10, norm_weight=0.2)
    model.train(group_indices_list, use_gpu=bool(torch.cuda.is_available()),
                max_epochs=max_epochs)

    A.obsm["latent_bg"], A.obsm["latent_t"] = model.get_latent_representation()
    treat = [c for c in conditions if c != control][0]
    resp = model.responsive_cells(A, treatment_condition=treat,
                                  control_condition=control)

    outdir.mkdir(parents=True, exist_ok=True)
    np.save(outdir / "sccausalvi_latent_bg.npy", A.obsm["latent_bg"])
    np.save(outdir / "sccausalvi_latent_te.npy", A.obsm["latent_t"])
    resp.obs.to_csv(outdir / "sccausalvi_responsive_cells.csv")
    return {"status": "ok",
            "n_background_latent": int(A.obsm["latent_bg"].shape[1]),
            "n_te_latent": int(A.obsm["latent_t"].shape[1]),
            "n_responsive_returned": int(resp.n_obs)}


# =============================================================================
# 6. 出图(全部非条形图:散点 / slopegraph / violin / lollipop)
# =============================================================================
def figures(reps, metrics, L, cond, ct, score, true_resp, cf, obs, d_ct, genes, outdir: Path):
    import matplotlib.pyplot as plt
    set_pub_style(base_size=10)
    outdir.mkdir(parents=True, exist_ok=True)
    cts = sorted(np.unique(ct))
    cols_ct = dict(zip(cts, pal(len(cts), "npg")))
    cols_cd = {"ctrl": "#4DBBD5", "stim": "#E64B35"}
    made = []

    # --- Fig 1:两种背景表示的散点(按条件/按细胞类型上色) ---
    fig, axes = plt.subplots(2, 2, figsize=(NATURE_W2, 6.2))
    for j, (name, Z) in enumerate([("Raw PCA (no disentanglement)", reps["raw"]),
                                   ("Condition-centered PCA (linear baseline)",
                                    reps["condition_centered"])]):
        for i, (lab, colmap, arr) in enumerate([("Condition", cols_cd, cond),
                                                ("Cell type", cols_ct, ct)]):
            ax = axes[i, j]
            for g in sorted(np.unique(arr)):
                m = arr == g
                ax.scatter(Z[m, 0], Z[m, 1], s=5, alpha=.65, lw=0,
                           color=colmap[g], label=str(g))
            ax.set_xlabel("PC1"); ax.set_ylabel("PC2")
            ax.set_title(f"{name}\ncoloured by {lab.lower()}", fontsize=8)
            ax.legend(markerscale=2.5, fontsize=7, loc="best")
    fig.tight_layout(); save_fig(fig, outdir / "fig1_background_embedding"); plt.close(fig)
    made.append("fig1_background_embedding.png")

    # --- Fig 2:解耦指标 slopegraph(raw -> centered) ---
    fig, ax = plt.subplots(figsize=(3.6, 4.0))
    keys = [("condition_leakage_auc", "Condition leakage AUC\n(lower is better, 0.5 = ideal)"),
            ("celltype_knn_acc", "Cell-type kNN accuracy\n(higher is better)")]
    xs = [0, 1]
    for (k, lab), c in zip(keys, pal(2, "npg")):
        ys = [metrics["raw"][k], metrics["condition_centered"][k]]
        ax.plot(xs, ys, "-o", color=c, lw=2, ms=8, label=lab)
        for x, y in zip(xs, ys):
            ax.annotate(f"{y:.3f}", (x, y), textcoords="offset points",
                        xytext=(0, 9), ha="center", fontsize=8)
    ax.axhline(0.5, ls=":", c="grey", lw=1)
    ax.set_xticks(xs); ax.set_xticklabels(["Raw PCA", "Condition-\ncentered"])
    ax.set_xlim(-.35, 1.35); ax.set_ylim(0.3, 1.08); ax.set_ylabel("Metric value")
    ax.set_title("Linear disentanglement baseline")
    ax.legend(fontsize=7, loc="lower left")
    fig.tight_layout(); save_fig(fig, outdir / "fig2_disentanglement_slopegraph"); plt.close(fig)
    made.append("fig2_disentanglement_slopegraph.png")

    # --- Fig 3:响应打分的 violin + 抖动点(按细胞类型 × 条件) ---
    fig, ax = plt.subplots(figsize=(NATURE_W2 * .78, 3.6))
    rng = np.random.default_rng(SEED)
    pos, labs = [], []
    p = 0
    for c in cts:
        for d in ["ctrl", "stim"]:
            v = score[(ct == c) & (cond == d)]
            parts = ax.violinplot([v], positions=[p], widths=.8, showextrema=False)
            for b in parts["bodies"]:
                b.set_facecolor(cols_cd[d]); b.set_alpha(.45); b.set_edgecolor("black"); b.set_lw(.6)
            ax.scatter(p + rng.normal(0, .07, v.size), v, s=3, alpha=.35, lw=0, color="black")
            ax.hlines(np.median(v), p - .3, p + .3, color="black", lw=1.6)
            pos.append(p); labs.append(f"{c}\n{d}"); p += 1
        p += .5
    ax.axhline(0, ls=":", c="grey", lw=1)
    ax.set_xticks(pos); ax.set_xticklabels(labs, fontsize=7)
    ax.set_ylabel("Perturbation score (z vs control null)")
    ax.set_title("Cell-state-specific response (CT_C is the non-responding control)")
    fig.tight_layout(); save_fig(fig, outdir / "fig3_response_score_violin"); plt.close(fig)
    made.append("fig3_response_score_violin.png")

    # --- Fig 4:反事实预测 散点(预测 Δ vs 观测 Δ,留出细胞) ---
    fig, axes = plt.subplots(1, 2, figsize=(NATURE_W2, 3.3), sharey=True)
    for ax, (k, ttl) in zip(axes, [("global", "Global Δ (condition mean shift)"),
                                   ("celltype", "Cell-type-specific Δ")]):
        ax.scatter(cf[k]["pred"], obs, s=6, alpha=.45, lw=0, color="#3C5488")
        lo = min(cf[k]["pred"].min(), obs.min()); hi = max(cf[k]["pred"].max(), obs.max())
        ax.plot([lo, hi], [lo, hi], ls="--", c="#E64B35", lw=1.2)
        ax.set_xlabel("Predicted Δ (held-out)"); ax.set_title(ttl, fontsize=9)
        ax.annotate(f"r = {cf[k]['pearson_r']:.3f}\nR² = {cf[k]['r2']:.3f}",
                    (.04, .92), xycoords="axes fraction", va="top", fontsize=8)
    axes[0].set_ylabel("Observed Δ (held-out)")
    fig.tight_layout(); save_fig(fig, outdir / "fig4_counterfactual_scatter"); plt.close(fig)
    made.append("fig4_counterfactual_scatter.png")

    # --- Fig 5:最强响应细胞类型的 top 基因 lollipop(不用条形图) ---
    top_ct = max(d_ct, key=lambda c: np.abs(d_ct[c]).sum())
    d = d_ct[top_ct]
    idx = np.argsort(-np.abs(d))[:20][::-1]
    fig, ax = plt.subplots(figsize=(3.9, 5.0))
    y = np.arange(len(idx))
    ax.hlines(y, 0, d[idx], color="#B0B7C3", lw=1.2)
    sc_ = ax.scatter(d[idx], y, c=np.abs(d[idx]), cmap=CMAP_CONT, s=46, zorder=3,
                     edgecolor="black", lw=.4)
    ax.axvline(0, color="black", lw=.8)
    ax.set_yticks(y); ax.set_yticklabels(genes[idx], fontsize=7)
    ax.set_xlabel("Treatment Δ (log-normalised)")
    ax.set_title(f"Top response genes · {top_ct}")
    fig.colorbar(sc_, ax=ax, label="|Δ|", fraction=.045, pad=.03)
    fig.tight_layout(); save_fig(fig, outdir / "fig5_top_response_genes_lollipop"); plt.close(fig)
    made.append("fig5_top_response_genes_lollipop.png")
    return made


# =============================================================================
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", default=str(EXAMPLE),
                    help=".npz with X / cell_type / condition / batch (默认合成示例)")
    ap.add_argument("--outdir", default=str(RESULTS))
    ap.add_argument("--assets", default=str(ASSETS))
    ap.add_argument("--control", default="ctrl", help="对照条件名")
    ap.add_argument("--n-pcs", type=int, default=20)
    ap.add_argument("--run-sccausalvi", action="store_true",
                    help="尝试真实 scCausalVI 路径(需 pip install scCausalVI)")
    ap.add_argument("--max-epochs", type=int, default=200)
    a = ap.parse_args()

    np.random.seed(SEED)
    outdir, assets = Path(a.outdir), Path(a.assets)
    outdir.mkdir(parents=True, exist_ok=True); assets.mkdir(parents=True, exist_ok=True)

    print("[588] Step 1  载入数据")
    X, ct, cond, batch, true_resp, genes = load_data(Path(a.input))
    print(f"       {X.shape[0]} cells × {X.shape[1]} genes · "
          f"conditions={sorted(set(cond))} · cell types={sorted(set(ct))}")

    print("[588] Step 2  归一化 (CPM+log1p)")
    L = lognorm(X)

    print("[588] Step 3  基线 A:背景表示与解耦指标(5 折交叉验证,防循环分析)")
    reps = background_representations(L, cond, n_pcs=a.n_pcs)
    metrics = {k: disentanglement_metrics(v, cond, ct) for k, v in reps.items()}
    for k, v in metrics.items():
        print(f"       {k:<19} leakage AUC={v['condition_leakage_auc']:.3f}  "
              f"cell-type kNN acc={v['celltype_knn_acc']:.3f}")

    print("[588] Step 4  基线 B:线性反事实预测(留出细胞评估)")
    cf, obs, ct_order, d_ct = counterfactual_baselines(L, cond, ct)
    for k, v in cf.items():
        print(f"       delta model = {k:<9} r={v['pearson_r']:.3f}  R2={v['r2']:.3f}")

    print("[588] Step 5  基线 C:响应细胞打分 + 与真值比 AUROC")
    score = responsive_score_baseline(L, cond, n_pcs=a.n_pcs)
    from sklearn.metrics import roc_auc_score
    m_stim = cond == "stim"
    auroc = float(roc_auc_score(true_resp[m_stim], score[m_stim])) if true_resp[m_stim].std() else float("nan")
    print(f"       responder AUROC (stim cells, truth = responding cell types) = {auroc:.3f}")

    print("[588] Step 6  出图")
    figs = figures(reps, metrics, L, cond, ct, score, true_resp, cf, obs, d_ct, genes, assets)
    for f in figs:
        print(f"       {f}")

    # ---- 落盘 ----
    import pandas as pd
    pd.DataFrame([{"representation": k, **v} for k, v in metrics.items()]).to_csv(
        outdir / "disentanglement_metrics.csv", index=False)
    pd.DataFrame({"cell_type": ct, "condition": cond, "batch": batch,
                  "perturbation_score": score, "true_responder": true_resp}).to_csv(
        outdir / "per_cell_perturbation_score.csv", index=False)
    pd.DataFrame({"gene": genes, **{f"delta_{c}": d_ct[c] for c in sorted(d_ct)}}).to_csv(
        outdir / "celltype_treatment_deltas.csv", index=False)

    # 依赖版本快照(铁律:可复现)
    import platform, sklearn, matplotlib
    session = {"python": platform.python_version(), "numpy": np.__version__,
               "pandas": pd.__version__, "scikit-learn": sklearn.__version__,
               "matplotlib": matplotlib.__version__}

    summary = {
        "module": "588_sccausalvi_causal_perturbation",
        "seed": SEED,
        "session_info": session,
        "n_cells": int(X.shape[0]), "n_genes": int(X.shape[1]),
        "baseline_disentanglement": metrics,
        "baseline_counterfactual": {k: {"pearson_r": v["pearson_r"], "r2": v["r2"]}
                                    for k, v in cf.items()},
        "baseline_responder_auroc": auroc,
        "figures": figs,
    }

    if a.run_sccausalvi:
        print("[588] Step 7  scCausalVI 真实路径")
        r = run_sccausalvi(X, cond, ct, control=a.control,
                           max_epochs=a.max_epochs, outdir=outdir)
        for k, v in r.items():
            print(f"       {k}: {v}")
        summary["sccausalvi"] = r
    else:
        summary["sccausalvi"] = {"status": "not requested (--run-sccausalvi)"}
        print("[588] Step 7  未请求 scCausalVI 路径,只跑基线")

    with open(outdir / "588_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=1, ensure_ascii=False, default=str)
    print(f"[588] done -> {outdir / '588_summary.json'}")


if __name__ == "__main__":
    main()
