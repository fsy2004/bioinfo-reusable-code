"""582 · D-SPIN — 从多重扰动 scRNA-seq 反推自旋网络(Ising)调控模型。

D-SPIN 把细胞状态写成基因程序上的 Ising 自旋模型:一个跨样本共享的耦合矩阵 J
(程序 i 与程序 j 的相互作用)+ 每个扰动样本自己的场向量 h(扰动响应向量)。
论文:Jiang J, et al. D-SPIN constructs regulatory network models from scRNA-seq
that reveal organizing principles of perturbation response. Cell 2026.
doi:10.1016/j.cell.2026.04.028 · PMID 42127893
仓库:https://github.com/JialongJiang/DSPIN   PyPI: https://pypi.org/project/dspin

本模块的结构(与库内 561 一致):
  · 基线(必跑,只用本机已有依赖):
      B0  朴素相关网络 —— 把所有细胞混在一起算程序状态的 Pearson 相关;
      B1  朴素平均场逆 Ising —— 先在每个样本内部算协方差再跨样本平均,
          再做 J_mf = -C^-1(Kappen & Rodríguez 1998 的 naive mean-field 逆问题;
          Nguyen, Zecchina & Berg, Adv Phys 2017 综述里的 nMF)。
      两者都不是 D-SPIN 的伪似然求解器,只是"直接耦合 vs 间接共变"这件事的
      最朴素对照。示例数据带真值,所以两条基线都会被真的评一遍分。
  · D-SPIN 正式路径(--run-dspin,需要 pip install dspin):
      守卫式调用,签名逐字来自上游源码 dspin/dspin.py(见 dspin_call 里的注释)。
      没装包就打印真实安装命令后跳过,不做任何静默降级。

出图统一走 modules/_framework/pubstyle.py。
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

HERE = os.path.dirname(os.path.abspath(__file__))
FRAMEWORK = os.path.abspath(os.path.join(HERE, "..", "..", "..", "_framework"))
sys.path.insert(0, FRAMEWORK)

from pubstyle import (  # noqa: E402  框架统一顶刊样式
    CMAP_DIVERGE, NATURE_W1, NATURE_W2, pal, save_fig, set_pub_style,
)

import matplotlib.pyplot as plt  # noqa: E402

EXAMPLE = os.path.join(HERE, "example_data")
RESULTS = os.path.join(HERE, "results")
ASSETS = os.path.join(HERE, "assets")
SEED = 582


# --------------------------------------------------------------------------- #
# 1. 读数据
# --------------------------------------------------------------------------- #
def load_inputs(expr_csv: str, meta_csv: str):
    """读 log 归一化表达(行=细胞)与细胞元数据。"""
    X = pd.read_csv(expr_csv, index_col=0)
    meta = pd.read_csv(meta_csv)
    if "cell_id" in meta.columns:
        meta = meta.set_index("cell_id").loc[X.index]
    need = {"sample_id", "batch", "if_control"}
    missing = need - set(meta.columns)
    if missing:
        sys.exit(f"cell_meta 缺列: {sorted(missing)}(D-SPIN 需要 sample_id/batch/if_control)")
    meta["if_control"] = meta["if_control"].astype(str).str.lower().isin(
        {"true", "1", "yes", "t"})
    return X, meta


# --------------------------------------------------------------------------- #
# 2. 基因程序发现 + 三态离散化(oNMF 的朴素替身)
# --------------------------------------------------------------------------- #
def discover_programs(X: pd.DataFrame, num_spin: int, seed: int = SEED):
    """sklearn NMF 分解出 num_spin 个基因程序,再离散成 {-1,0,+1}。

    D-SPIN 正式流程用的是重复 oNMF(正交 NMF)+ KMeans 三态离散
    (ProgramDSPIN.gene_program_discovery / AbstractDSPIN.discretize)。
    这里用普通 NMF 作为只依赖 sklearn 的朴素替身(已明确标注,不冒充 oNMF),
    但离散化沿用上游的 KMeans 三态思路 —— 固定分位数切分会强行让三个状态各占
    1/3,把程序本身的活性不平衡抹掉;KMeans 是数据自适应的。
    """
    from sklearn.cluster import KMeans
    from sklearn.decomposition import NMF

    nmf = NMF(n_components=num_spin, init="nndsvda", random_state=seed,
              max_iter=800, tol=1e-5)
    A = nmf.fit_transform(X.values)          # 细胞 × 程序 活性
    W = nmf.components_                      # 程序 × 基因 载荷

    # 三态离散:每个程序独立跑 KMeans(k=3),按簇心高低映射为 -1 / 0 / +1
    S = np.zeros_like(A)
    for p in range(num_spin):
        v = A[:, p].reshape(-1, 1)
        km = KMeans(n_clusters=3, n_init=10, random_state=seed).fit(v)
        rank = np.argsort(np.argsort(km.cluster_centers_.ravel()))  # 簇 → 0/1/2
        S[:, p] = (rank[km.labels_] - 1).astype(float)
    prog_names = [f"P{p + 1}" for p in range(num_spin)]
    return (pd.DataFrame(A, index=X.index, columns=prog_names),
            pd.DataFrame(S, index=X.index, columns=prog_names),
            pd.DataFrame(W, index=prog_names, columns=X.columns))


def top_genes(W: pd.DataFrame, k: int = 5) -> dict:
    """每个程序的前 k 个高载荷基因,用于给程序取一个可读名字。"""
    return {p: list(W.loc[p].sort_values(ascending=False).index[:k]) for p in W.index}


# --------------------------------------------------------------------------- #
# 3. 基线 B0 / B1
# --------------------------------------------------------------------------- #
def baseline_correlation(S: pd.DataFrame) -> np.ndarray:
    """B0 朴素相关网络:所有细胞混池的 Pearson 相关,忽略样本结构。"""
    C = np.corrcoef(S.values.T)
    np.fill_diagonal(C, 0.0)
    return C


def baseline_mean_field(S: pd.DataFrame, sample_id: pd.Series):
    """B1 朴素平均场逆 Ising:样本内协方差跨样本平均后求逆取负。

    J_mf = -C^{-1}(off-diagonal);场 h_s 由每个样本的平均自旋经平均场方程反解。
    先在样本内算协方差,是为了不把"不同扰动把整体活性推到不同水平"这种
    跨样本的场驱动共变误当成耦合——这正是朴素相关网络的主要失真来源。
    """
    P = S.shape[1]
    covs, ws = [], []
    for sid, idx in S.groupby(sample_id.values).groups.items():
        sub = S.loc[idx].values
        if sub.shape[0] < P + 2:            # 样本太小,协方差不可靠
            continue
        covs.append(np.cov(sub, rowvar=False))
        ws.append(sub.shape[0])
    C = np.average(np.array(covs), axis=0, weights=np.array(ws, float))
    C_reg = C + np.eye(P) * (1e-3 * np.trace(C) / P)     # 轻正则,保证可逆
    J = -np.linalg.inv(C_reg)
    np.fill_diagonal(J, 0.0)

    # 每个样本的场:h_i = arctanh(<s_i>) - sum_j J_ij <s_j>(平均场自洽方程)
    h = {}
    for sid, idx in S.groupby(sample_id.values).groups.items():
        m = np.clip(S.loc[idx].values.mean(axis=0), -0.98, 0.98)
        h[sid] = np.arctanh(m) - J @ m
    H = pd.DataFrame(h, index=S.columns).T
    return J, H


def response_relative_to_control(H: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    """扰动响应向量 = 样本场 - 同 batch 对照样本场的均值。

    对应 D-SPIN 的 AbstractDSPIN.response_relative_to_control(sample_id_key,
    if_control_key, batch_key);这里用同样的"按 batch 减对照"逻辑作用在
    平均场反解出的 h 上。
    """
    smeta = meta.groupby("sample_id").agg(batch=("batch", "first"),
                                          if_control=("if_control", "first"))
    out = {}
    for sid in H.index:
        b = smeta.loc[sid, "batch"]
        ctrl = smeta[(smeta.batch == b) & (smeta.if_control)].index
        ref = H.loc[[c for c in ctrl if c in H.index]].mean(axis=0) if len(ctrl) else 0.0
        out[sid] = H.loc[sid] - ref
    R = pd.DataFrame(out).T
    return R.loc[[s for s in H.index if not smeta.loc[s, "if_control"]]]


# --------------------------------------------------------------------------- #
# 4. 把恢复出的程序对齐到真值程序(只有示例数据才有真值)
# --------------------------------------------------------------------------- #
def align_to_truth(W: pd.DataFrame, n_true: int, genes_per: int):
    """按基因载荷把 NMF 程序匹配到真值程序(匈牙利算法最大化载荷相关)。"""
    from scipy.optimize import linear_sum_assignment

    truth_load = np.zeros((n_true, W.shape[1]))
    for p in range(n_true):
        truth_load[p, p * genes_per:(p + 1) * genes_per] = 1.0
    cost = -np.corrcoef(np.vstack([W.values, truth_load]))[:W.shape[0], W.shape[0]:]
    r, c = linear_sum_assignment(cost)
    order = np.empty(n_true, dtype=int)
    order[c] = r                       # order[真值程序] = 恢复程序下标
    return order, float(-cost[r, c].mean())


def score_network(J_est: np.ndarray, J_true: np.ndarray) -> dict:
    """两项评分:与真值耦合的 Spearman;把非零真值边当正类的 AUROC。"""
    from scipy.stats import spearmanr
    from sklearn.metrics import roc_auc_score

    iu = np.triu_indices_from(J_true, k=1)
    est, tru = J_est[iu], J_true[iu]
    rho = float(spearmanr(est, tru).correlation)
    y = (np.abs(tru) > 1e-9).astype(int)
    auc = float(roc_auc_score(y, np.abs(est))) if 0 < y.sum() < len(y) else float("nan")
    return {"spearman_vs_truth": round(rho, 3), "edge_auroc": round(auc, 3)}


# --------------------------------------------------------------------------- #
# 5. D-SPIN 正式路径(守卫式)
# --------------------------------------------------------------------------- #
def dspin_call(X: pd.DataFrame, meta: pd.DataFrame, num_spin: int, outdir: str) -> dict:
    """真正的 D-SPIN。没装包就带上真实安装命令退出,绝不静默降级。

    以下调用签名逐字取自上游源码(2026-07-20 读取):
      https://raw.githubusercontent.com/JialongJiang/DSPIN/main/dspin/dspin.py
      https://raw.githubusercontent.com/JialongJiang/DSPIN/main/dspin/plot.py
    · DSPIN(adata, save_path, num_spin=None, filter_threshold=0.02, **kwargs)
      —— __new__ 按基因数自动派发 GeneDSPIN / ProgramDSPIN
    · ProgramDSPIN.gene_program_discovery(num_repeat=10, seed=0,
          cluster_key='leiden', mode='compute_summary', prior_programs=None,
          params={}, discretize_params={})
    · AbstractDSPIN.network_inference(sample_id_key='sample_id', method='auto',
          directed=False, params=None, sample_list_ordered=None,
          prior_network=None, perturb_matrix=None, if_control_key='if_control',
          batch_key='batch', run_with_matlab=False)
      method 可选 'maximum_likelihood' / 'mcmc_maximum_likelihood' /
      'pseudo_likelihood' / 'auto'(默认)
    · AbstractDSPIN.response_relative_to_control(sample_id_key='sample_id',
          if_control_key='if_control', batch_key='batch')
    · 属性:model.network, model.responses, model.relative_responses,
          model.program_representation, model.sample_list, model.name_list
    · dspin.plot:plot_network_heatmap, plot_response_heatmap,
          plot_network_diagram, create_undirected_network, compute_modules
    params 的默认键由 AbstractDSPIN.default_params(method) 生成
    (dspin/dspin.py:206-256):num_epoch / cur_j / cur_h / save_path / rec_gap /
    seed / save_log / lambda_l1_j / lambda_l1_h / lambda_l2_j / lambda_l2_h /
    backtrack_gap / backtrack_tol,再按 method 追加 stepsz
    (ML 0.2 / MCMC-ML 0.02 + mcmc_* 三项 / PL 0.05)。
    本模块不覆写这些超参,一律用上游默认值。
    """
    try:
        import dspin
    except ImportError:
        return {"status": "skipped",
                "reason": "dspin 未安装",
                "install": "pip install dspin  # 需要 anndata scanpy matplotlib tqdm igraph leidenalg"}
    try:
        import anndata as ad
        import scanpy as sc

        adata = ad.AnnData(X.values.astype("float32"))
        adata.var_names = list(X.columns)
        adata.obs_names = list(X.index)
        for k in ("sample_id", "batch", "if_control"):
            adata.obs[k] = meta[k].values
        sc.pp.pca(adata, n_comps=min(30, adata.n_vars - 1))
        sc.pp.neighbors(adata)
        sc.tl.leiden(adata, key_added="leiden")

        save_path = os.path.join(outdir, "dspin_run")
        os.makedirs(save_path, exist_ok=True)
        model = dspin.DSPIN(adata, save_path, num_spin=num_spin)
        if hasattr(model, "gene_program_discovery"):        # ProgramDSPIN 分支
            model.gene_program_discovery(num_repeat=10, seed=SEED, cluster_key="leiden")
        model.network_inference(sample_id_key="sample_id", method="pseudo_likelihood",
                                if_control_key="if_control", batch_key="batch")
        model.response_relative_to_control(sample_id_key="sample_id",
                                           if_control_key="if_control",
                                           batch_key="batch")
        np.savetxt(os.path.join(outdir, "dspin_network_J.csv"),
                   np.asarray(model.network), delimiter=",")
        np.savetxt(os.path.join(outdir, "dspin_relative_responses.csv"),
                   np.asarray(model.relative_responses), delimiter=",")
        return {"status": "ok",
                "dspin_version": getattr(dspin, "__version__", "?"),
                "network_shape": list(np.asarray(model.network).shape),
                "class": type(model).__name__}
    except Exception as e:                                   # 上游报错如实抛出,不吞
        return {"status": "error", "error": f"{type(e).__name__}: {e}",
                "hint": "上游官方 Colab demo 见仓库 README.md 的 Demos 一节:"
                        "https://github.com/JialongJiang/DSPIN"}


# --------------------------------------------------------------------------- #
# 6. 出图(无条形图:heatmap / 散点 / 网络图 / lollipop / dumbbell)
# --------------------------------------------------------------------------- #
def _hm(ax, M, names, vmax, title):
    im = ax.imshow(M, cmap=CMAP_DIVERGE, vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(len(names))); ax.set_xticklabels(names, rotation=90)
    ax.set_yticks(range(len(names))); ax.set_yticklabels(names)
    ax.set_title(title)
    for s in ax.spines.values():
        s.set_visible(True)
    return im


def fig_network_heatmaps(mats: dict, names, path):
    """并排热图:真值耦合 vs 朴素相关 vs 平均场耦合。"""
    n = len(mats)
    fig, axes = plt.subplots(1, n, figsize=(NATURE_W2, NATURE_W2 / n + 0.9))
    axes = np.atleast_1d(axes)
    for ax, (t, M) in zip(axes, mats.items()):
        Z = M / (np.abs(M).max() + 1e-12)                 # 各自标准化到 [-1,1] 便于比形状
        im = _hm(ax, Z, names, 1.0, t)
    fig.colorbar(im, ax=axes.tolist(), fraction=0.02, pad=0.02,
                 label="Scaled coupling")
    save_fig(fig, path); plt.close(fig)


def fig_edge_scatter(J_corr, J_mf, J_true, scores, path):
    """散点:每条程序对的真值耦合 vs 推断值,直接边 / 间接对分色。"""
    iu = np.triu_indices_from(J_true, k=1)
    tru = J_true[iu]
    direct = np.abs(tru) > 1e-9
    cols = pal(3, "npg")
    fig, axes = plt.subplots(1, 2, figsize=(NATURE_W2, NATURE_W1 * 0.95), sharex=True)
    for ax, (name, J) in zip(axes, [("Pooled correlation (B0)", J_corr),
                                    ("Mean-field Ising (B1)", J_mf)]):
        est = J[iu]
        est = est / (np.abs(est).max() + 1e-12)
        ax.axhline(0, color="#BBBBBB", lw=0.8, zorder=0)
        ax.axvline(0, color="#BBBBBB", lw=0.8, zorder=0)
        ax.scatter(tru[~direct], est[~direct], s=26, c=cols[1], alpha=.75,
                   edgecolor="none", label="Indirect pair (true J = 0)")
        ax.scatter(tru[direct], est[direct], s=42, c=cols[0], alpha=.9,
                   edgecolor="white", linewidth=.5, label="Direct edge")
        s = scores[name.split(" (")[0]]
        ax.set_title(f"{name}\nrho={s['spearman_vs_truth']}  AUROC={s['edge_auroc']}",
                     fontsize=9)
        ax.set_xlabel("True coupling")
        ax.set_ylabel("Inferred (scaled)")
    axes[0].legend(loc="upper left", fontsize=7)
    save_fig(fig, path); plt.close(fig)


def fig_auroc_dumbbell(scores, path):
    """Dumbbell:两条基线在同一指标上的位置对比(不用条形图)。"""
    metrics = ["spearman_vs_truth", "edge_auroc"]
    labels = ["Spearman vs true J", "Edge-detection AUROC"]
    cols = pal(3, "npg")
    fig, ax = plt.subplots(figsize=(NATURE_W1 * 1.35, NATURE_W1 * 0.7))
    for i, (m, lab) in enumerate(zip(metrics, labels)):
        a = scores["Pooled correlation"][m]
        b = scores["Mean-field Ising"][m]
        ax.plot([a, b], [i, i], color="#9A9A9A", lw=2, zorder=1, solid_capstyle="round")
        ax.scatter([a], [i], s=110, c=cols[1], zorder=2, edgecolor="white", lw=1)
        ax.scatter([b], [i], s=110, c=cols[0], zorder=2, edgecolor="white", lw=1)
        ax.text(b, i + .22, f"{b:.2f}", ha="center", fontsize=8, color=cols[0])
        ax.text(a, i - .32, f"{a:.2f}", ha="center", fontsize=8, color=cols[1])
    ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels)
    ax.set_ylim(-.6, len(labels) - .4)
    ax.set_xlabel("Score (higher = closer to true direct coupling)")
    ax.scatter([], [], s=90, c=cols[1], label="B0 pooled correlation")
    ax.scatter([], [], s=90, c=cols[0], label="B1 mean-field Ising")
    ax.legend(loc="lower right", fontsize=7)
    save_fig(fig, path); plt.close(fig)


def fig_response_heatmap(R: pd.DataFrame, path):
    """扰动 × 程序 的相对响应向量热图。"""
    M = R.values
    v = np.abs(M).max()
    fig, ax = plt.subplots(figsize=(NATURE_W1 * 1.5, NATURE_W1 * 1.05))
    im = ax.imshow(M, cmap=CMAP_DIVERGE, vmin=-v, vmax=v, aspect="auto")
    ax.set_xticks(range(R.shape[1])); ax.set_xticklabels(R.columns, rotation=90)
    ax.set_yticks(range(R.shape[0])); ax.set_yticklabels(R.index)
    ax.set_title("Perturbation response vectors (relative to batch control)")
    for s in ax.spines.values():
        s.set_visible(True)
    fig.colorbar(im, ax=ax, fraction=0.035, pad=0.03, label="Field shift h")
    save_fig(fig, path); plt.close(fig)


def fig_response_lollipop(R: pd.DataFrame, path, top: int = 12):
    """Lollipop:响应幅度最大的 程序×扰动 组合。"""
    L = R.stack().rename("h").reset_index()
    L.columns = ["sample", "program", "h"]
    L["lab"] = L["sample"] + " · " + L["program"]
    L = L.reindex(L.h.abs().sort_values(ascending=False).index)[:top].iloc[::-1]
    cols = pal(3, "npg")
    c = [cols[0] if v > 0 else cols[1] for v in L.h]
    fig, ax = plt.subplots(figsize=(NATURE_W1 * 1.5, NATURE_W1 * 1.2))
    ax.axvline(0, color="#BBBBBB", lw=0.8)
    ax.hlines(range(len(L)), 0, L.h, color="#B0B0B0", lw=1.6)
    ax.scatter(L.h, range(len(L)), s=80, c=c, zorder=3, edgecolor="white", lw=.8)
    ax.set_yticks(range(len(L))); ax.set_yticklabels(L.lab, fontsize=8)
    ax.set_xlabel("Relative response h (vs batch control)")
    ax.set_title("Strongest perturbation responses")
    save_fig(fig, path); plt.close(fig)


def fig_network_diagram(J, names, R, path, thres=0.15):
    """网络图:边=自旋耦合(红正蓝负,粗细=强度),节点大小=被扰动幅度。"""
    import networkx as nx

    Z = J / (np.abs(J).max() + 1e-12)
    G = nx.Graph()
    G.add_nodes_from(names)
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            if abs(Z[i, j]) >= thres:
                G.add_edge(names[i], names[j], w=Z[i, j])
    pos = nx.spring_layout(G, seed=SEED, weight=None, k=1.1)
    amp = R.abs().max(axis=0).reindex(names).fillna(0).values
    size = 380 + 900 * amp / (amp.max() + 1e-12)
    cols = pal(3, "npg")

    fig, ax = plt.subplots(figsize=(NATURE_W1 * 1.5, NATURE_W1 * 1.35))
    for u, v, d in G.edges(data=True):
        ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]],
                color=cols[0] if d["w"] > 0 else cols[1],
                lw=0.6 + 4.0 * abs(d["w"]), alpha=.75, zorder=1,
                solid_capstyle="round")
    ax.scatter([pos[n][0] for n in names], [pos[n][1] for n in names],
               s=size, c="white", edgecolor="#333333", lw=1.2, zorder=2)
    for n in names:
        ax.text(pos[n][0], pos[n][1], n, ha="center", va="center",
                fontsize=8, zorder=3)
    ax.plot([], [], color=cols[0], lw=3, label="Positive coupling")
    ax.plot([], [], color=cols[1], lw=3, label="Negative coupling")
    # 图例放画布外下方,避免遮住 spring 布局落在角落的节点
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.01), ncol=2,
              fontsize=7, frameon=False)
    ax.set_axis_off()
    ax.set_title("Inferred spin-network (node size = perturbation amplitude)")
    save_fig(fig, path); plt.close(fig)


# --------------------------------------------------------------------------- #
# 7. main
# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--expr", default=os.path.join(EXAMPLE, "expression.csv"),
                    help="log 归一化表达 csv,行=细胞 列=基因")
    ap.add_argument("--meta", default=os.path.join(EXAMPLE, "cell_meta.csv"),
                    help="细胞元数据 csv,需含 sample_id / batch / if_control")
    ap.add_argument("--truth", default=os.path.join(EXAMPLE, "ground_truth_coupling.csv"),
                    help="真值耦合矩阵(可选;有才做基线评分)")
    ap.add_argument("--num-spin", type=int, default=8, help="基因程序 / 自旋数")
    ap.add_argument("--genes-per-program", type=int, default=6,
                    help="仅用于示例数据的程序↔真值对齐")
    ap.add_argument("--run-dspin", action="store_true",
                    help="尝试真正的 D-SPIN(需 pip install dspin)")
    ap.add_argument("--outdir", default=RESULTS)
    ap.add_argument("--seed", type=int, default=SEED)
    a = ap.parse_args()

    np.random.seed(a.seed)
    os.makedirs(a.outdir, exist_ok=True)
    os.makedirs(ASSETS, exist_ok=True)
    set_pub_style(base_size=9)

    print("[582] Step 1  读入表达与细胞元数据")
    X, meta = load_inputs(a.expr, a.meta)
    print(f"       {X.shape[0]} cells × {X.shape[1]} genes · "
          f"{meta.sample_id.nunique()} samples · {meta.batch.nunique()} batches")

    print(f"[582] Step 2  基因程序发现(NMF,{a.num_spin} 个程序)+ 三态离散")
    A, S, W = discover_programs(X, a.num_spin, a.seed)
    names = list(S.columns)
    tg = top_genes(W)
    pd.DataFrame({p: ",".join(g) for p, g in tg.items()}, index=["top_genes"]).T.to_csv(
        os.path.join(a.outdir, "program_top_genes.csv"))
    S.to_csv(os.path.join(a.outdir, "program_states.csv"))

    print("[582] Step 3  基线 B0 朴素相关网络(混池 Pearson)")
    J_corr = baseline_correlation(S)

    print("[582] Step 4  基线 B1 朴素平均场逆 Ising(样本内协方差 → -C^-1)")
    J_mf, H = baseline_mean_field(S, meta.sample_id)
    pd.DataFrame(J_corr, index=names, columns=names).to_csv(
        os.path.join(a.outdir, "baseline_correlation_network.csv"))
    pd.DataFrame(J_mf, index=names, columns=names).to_csv(
        os.path.join(a.outdir, "baseline_meanfield_network.csv"))

    print("[582] Step 5  扰动响应向量(减同 batch 对照)")
    R = response_relative_to_control(H, meta)
    R.to_csv(os.path.join(a.outdir, "relative_response_vectors.csv"))

    summary = {"n_cells": int(X.shape[0]), "n_genes": int(X.shape[1]),
               "num_spin": a.num_spin, "n_samples": int(meta.sample_id.nunique()),
               "seed": a.seed}

    scores, J_true_ord = None, None
    if a.truth and os.path.exists(a.truth):
        print("[582] Step 6  对齐到真值程序并给两条基线评分")
        J_true = pd.read_csv(a.truth, index_col=0).values
        order, load_r = align_to_truth(W, J_true.shape[0], a.genes_per_program)
        summary["program_matching_loading_r"] = round(load_r, 3)
        Jc = J_corr[np.ix_(order, order)]
        Jm = J_mf[np.ix_(order, order)]
        scores = {"Pooled correlation": score_network(Jc, J_true),
                  "Mean-field Ising": score_network(Jm, J_true)}
        summary["baseline_scores"] = scores
        for k, v in scores.items():
            print(f"       {k:22s} {v}")
        J_true_ord = J_true
        J_corr_ord, J_mf_ord = Jc, Jm
        # 响应向量也按真值程序顺序重排,便于和 ground_truth_response 对照
        R = R.iloc[:, order]
        R.columns = [f"P{i + 1}" for i in range(J_true.shape[0])]
        names = list(R.columns)
        J_corr, J_mf = Jc, Jm
    else:
        print("[582] Step 6  未提供真值,跳过基线评分(真实数据的正常情况)")

    print("[582] Step 7  出图")
    mats = {}
    if J_true_ord is not None:
        mats["True coupling"] = J_true_ord
    mats["B0 pooled correlation"] = J_corr
    mats["B1 mean-field Ising"] = J_mf
    fig_network_heatmaps(mats, names, os.path.join(ASSETS, "582_fig1_network_heatmap"))
    if scores:
        fig_edge_scatter(J_corr, J_mf, J_true_ord, scores,
                         os.path.join(ASSETS, "582_fig2_edge_recovery_scatter"))
        fig_auroc_dumbbell(scores, os.path.join(ASSETS, "582_fig3_baseline_dumbbell"))
    fig_response_heatmap(R, os.path.join(ASSETS, "582_fig4_response_heatmap"))
    fig_response_lollipop(R, os.path.join(ASSETS, "582_fig5_response_lollipop"))
    fig_network_diagram(J_mf, names, R, os.path.join(ASSETS, "582_fig6_network_diagram"))

    if a.run_dspin:
        print("[582] Step 8  D-SPIN 正式路径")
        d = dspin_call(X, meta, a.num_spin, a.outdir)
        summary["dspin"] = d
        for k, v in d.items():
            print(f"       {k}: {v}")
    else:
        print("[582] Step 8  未请求 D-SPIN 正式路径(--run-dspin);只跑基线")
        summary["dspin"] = {"status": "not_requested"}

    with open(os.path.join(a.outdir, "582_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=1, ensure_ascii=False, default=str)
    print(f"[582] 完成 · 结果 {a.outdir} · 展示图 {ASSETS}")


if __name__ == "__main__":
    main()
