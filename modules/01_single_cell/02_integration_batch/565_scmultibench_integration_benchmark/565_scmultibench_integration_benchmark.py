# -*- coding: utf-8 -*-
"""565 · scMultiBench — scoring single-cell multimodal integrations with benchmark metrics.

把 scMultiBench (Liu et al., Nature Methods 2025) 的 **评测层** 封装成可复用模块:
输入 = 整合后的低维 embedding + celltype/batch 标签,输出 = 同一套 scIB 指标
(bio conservation / batch correction / overall)+ 顶刊风对比图。

两条路径:
  1) **本机路径(默认,零依赖门槛)**:用 numpy/sklearn/networkx/scipy 按公开公式
     重算 scIB 指标(ARI/NMI/ASW_label/ASW_batch/graph connectivity/cLISI/iLISI/kBET)。
     这是**近似实现**,不是 scib 包的原始代码,数值可能与官方管线有偏差 —— README 已注明。
  2) **上游路径(--use-scib,守卫式)**:若本机装了 `scib`,直接调用 scMultiBench
     evaluation_pipelines/scib_metrics/scib_metrics.py 中**逐字核对过**的 scib.me.* 调用
     (2026-07-20 又对 theislab/scib 各指标的 def 签名逐个复核了参数名)。
     未装则打印真实安装命令后跳过,绝不静默降级。
     ⚠ **未执行验证**:本机没装 scib(也没有 kBET 所需的 R 包 + rpy2),这条路径只做过
     静态 API 核对,**从未真正跑通过**;每个指标单独 try,失败只记 <指标>_error。

必带基线:PCA-on-concatenated-features(不做任何批次校正)。任何"整合更好"的说法
都必须先赢过它。

上游仓库 : https://github.com/PYangLab/scMultiBench (Apache-2.0)
API 读取自: https://raw.githubusercontent.com/PYangLab/scMultiBench/main/evaluation_pipelines/scib_metrics/scib_metrics.py
论文     : Liu C, Ding S, Kim HJ, Long S, Xiao D, Ghazanfar S, Yang P.
           Multitask benchmarking of single-cell multimodal omics integration methods.
           Nat Methods. 2025;22(11):2449-2460. doi:10.1038/s41592-025-02856-3 · PMID 41083898
指标定义 : Luecken MD, et al. Benchmarking atlas-level data integration in single-cell
           genomics. Nat Methods. 2022;19(1):41-50. PMID 34949812 (scIB)
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import warnings

warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
EXAMPLE = os.path.join(HERE, "example_data")
RESULTS = os.path.join(HERE, "results")
ASSETS = os.path.join(HERE, "assets")
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..", "..", "_framework")))

SEED = 2026

# scIB 分组:哪些指标算生物保留,哪些算批次校正(Luecken 2022 Table 1 的划分)
BIO_METRICS = ["ARI_cluster_label", "NMI_cluster_label", "ASW_label", "cLISI",
               "isolated_label_ASW"]
BATCH_METRICS = ["ASW_batch", "graph_connectivity", "iLISI", "kBET_acceptance"]
W_BIO, W_BATCH = 0.6, 0.4          # scIB 的 overall 加权(0.6 bio / 0.4 batch)


# ---------------------------------------------------------------- 输入
def load_example():
    """读 example_data/ 的合成 RNA+ADT 计数与标签(synthetic, for demo only)。"""
    import pandas as pd
    rna = pd.read_csv(os.path.join(EXAMPLE, "rna_counts.csv"), comment="#", index_col=0)
    adt = pd.read_csv(os.path.join(EXAMPLE, "adt_counts.csv"), comment="#", index_col=0)
    meta = pd.read_csv(os.path.join(EXAMPLE, "metadata.csv"), comment="#", index_col=0)
    meta = meta.loc[rna.index]
    return {"rna": rna, "adt": adt}, meta


def _lognorm(df):
    """CPM-style 归一 + log1p;各模态独立做,再 z-score(防某一模态量纲主导 PCA)。"""
    import numpy as np
    X = df.values.astype(float)
    tot = X.sum(1, keepdims=True)
    tot[tot == 0] = 1.0
    X = np.log1p(X / tot * 1e4)
    sd = X.std(0, keepdims=True)
    sd[sd == 0] = 1.0
    return (X - X.mean(0, keepdims=True)) / sd


# ---------------------------------------------------------------- 被评测的 embedding
def build_embeddings(mods, meta, n_comp=20):
    """产出待评测的 embedding 字典。第一个是必跑的朴素基线。"""
    import numpy as np
    from sklearn.decomposition import PCA

    Z = np.hstack([_lognorm(df) for df in mods.values()])          # 各模态拼接
    embs = {}

    # 基线 ①:直接对拼接特征做 PCA,完全不做批次校正 —— 本模块的对照下限
    embs["PCA_concat (baseline)"] = PCA(n_components=n_comp, random_state=SEED).fit_transform(Z)

    # 基线 ②:批次内中心化(最朴素的线性批次校正),看"最简单的校正能拿多少分"
    Zc = Z.copy()
    for b in meta["batch"].unique():
        m = (meta["batch"] == b).values
        Zc[m] -= Zc[m].mean(0, keepdims=True)
    embs["PCA_batch_centered"] = PCA(n_components=n_comp, random_state=SEED).fit_transform(Zc)

    # 基线 ③:批次内 z-score(中心化 + 方差归一),常见的"过校正"对照
    Zs = Z.copy()
    for b in meta["batch"].unique():
        m = (meta["batch"] == b).values
        sd = Zs[m].std(0, keepdims=True)
        sd[sd == 0] = 1.0
        Zs[m] = (Zs[m] - Zs[m].mean(0, keepdims=True)) / sd
    embs["PCA_batch_zscore"] = PCA(n_components=n_comp, random_state=SEED).fit_transform(Zs)

    # 可选 ④:harmonypy(装了才跑,没装跳过并说明)
    try:
        import harmonypy
        ho = harmonypy.run_harmony(embs["PCA_concat (baseline)"], meta, ["batch"])
        # ⚠ Z_corr 的朝向随 harmonypy 版本变化:旧的纯 Python 版是 (d x N),
        # 现在的 C++ 后端版 Z_corr property 返回 self._cpp.Z_corr.T 即 (N x d)。
        # 不能写死 .T —— 按 cell 数判断朝向。
        Zh = np.asarray(ho.Z_corr)
        if Zh.shape[0] != Z.shape[0]:
            Zh = Zh.T
        if Zh.shape[0] != Z.shape[0]:
            raise ValueError(f"harmonypy 返回形状 {np.asarray(ho.Z_corr).shape} 与 {Z.shape[0]} 细胞对不上")
        embs["Harmony"] = Zh
    except Exception as e:
        print(f"       [skip] Harmony 未评测({type(e).__name__});pip install harmonypy 后可加入对比")

    return embs


# ---------------------------------------------------------------- 指标(本机近似实现)
def _knn_idx(emb, k=15):
    from sklearn.neighbors import NearestNeighbors
    k = min(k, emb.shape[0] - 1)
    nn = NearestNeighbors(n_neighbors=k + 1).fit(emb)
    return nn.kneighbors(emb, return_distance=False)[:, 1:]        # 去掉自身


def _simpson_scores(emb, labels, k=15):
    """kNN 邻域的 inverse Simpson index(LISI 的离散近似,未做 perplexity 加权)。"""
    import numpy as np
    idx = _knn_idx(emb, k)
    lab = np.asarray(labels)
    uniq = np.unique(lab)
    code = np.searchsorted(uniq, lab)
    neigh = code[idx]                                              # cells x k
    out = np.empty(len(lab))
    for i in range(len(lab)):
        p = np.bincount(neigh[i], minlength=len(uniq)) / neigh.shape[1]
        out[i] = 1.0 / np.sum(p ** 2)
    return out, len(uniq)


def compute_metrics_local(emb, celltype, batch, k=15):
    """按 scIB 公开公式在本机重算指标。近似实现,非 scib 原始代码。"""
    import numpy as np
    import networkx as nx
    from sklearn.cluster import KMeans
    from sklearn.metrics import (adjusted_rand_score, normalized_mutual_info_score,
                                 silhouette_score, silhouette_samples)
    from scipy.stats import chisquare

    ct = np.asarray(celltype)
    bt = np.asarray(batch)
    n_ct, n_bt = len(np.unique(ct)), len(np.unique(bt))
    m = {}

    # 聚类:scIB 用 optimised Leiden;本机 leidenalg 不保证可用,这里用 KMeans(k=真实类型数)
    # —— 这是与官方管线的已知偏差,README 已注明。
    clus = KMeans(n_clusters=n_ct, n_init=10, random_state=SEED).fit_predict(emb)
    m["ARI_cluster_label"] = float(adjusted_rand_score(ct, clus))
    m["NMI_cluster_label"] = float(normalized_mutual_info_score(ct, clus))

    # ASW(label):scIB 把 silhouette 缩放到 [0,1]
    m["ASW_label"] = float((silhouette_score(emb, ct) + 1) / 2)

    # ASW(batch):每个 celltype 内,对 batch 求 silhouette,取 1-|s| 的均值(越大=批次越混匀)
    vals = []
    for c in np.unique(ct):
        sel = ct == c
        if len(np.unique(bt[sel])) < 2 or sel.sum() < 10:
            continue
        s = silhouette_samples(emb[sel], bt[sel])
        vals.append(np.mean(1 - np.abs(s)))
    m["ASW_batch"] = float(np.mean(vals)) if vals else float("nan")

    # isolated label ASW:细胞数最少的类型的 ASW(scIB 的 isolated-label 思路的简化版)
    rare = min(np.unique(ct), key=lambda c: (ct == c).sum())
    s_all = silhouette_samples(emb, ct)
    m["isolated_label_ASW"] = float((np.mean(s_all[ct == rare]) + 1) / 2)

    # graph connectivity:每个 celltype 的 kNN 子图中最大连通分量占比
    idx = _knn_idx(emb, k)
    G = nx.Graph()
    G.add_nodes_from(range(emb.shape[0]))
    G.add_edges_from((i, j) for i in range(idx.shape[0]) for j in idx[i])
    gc = []
    for c in np.unique(ct):
        nodes = np.where(ct == c)[0]
        sub = G.subgraph(nodes)
        if sub.number_of_nodes() == 0:
            continue
        gc.append(len(max(nx.connected_components(sub), key=len)) / sub.number_of_nodes())
    m["graph_connectivity"] = float(np.mean(gc))

    # cLISI / iLISI:inverse Simpson 缩放到 [0,1](1 = 好)
    c_raw, nc = _simpson_scores(emb, ct, k)
    i_raw, nb = _simpson_scores(emb, bt, k)
    m["cLISI"] = float((nc - np.median(c_raw)) / (nc - 1)) if nc > 1 else float("nan")
    m["iLISI"] = float((np.median(i_raw) - 1) / (nb - 1)) if nb > 1 else float("nan")

    # kBET(近似):每个细胞的邻域 batch 组成 vs 全局组成做卡方,报"未被拒绝"的比例
    glob = np.array([(bt == b).sum() for b in np.unique(bt)], float) / len(bt)
    order = np.unique(bt)
    code = np.searchsorted(order, bt)
    neigh = code[idx]
    acc = 0
    for i in range(neigh.shape[0]):
        obs = np.bincount(neigh[i], minlength=n_bt).astype(float)
        exp = glob * obs.sum()
        exp[exp <= 0] = 1e-9
        acc += chisquare(obs, exp).pvalue > 0.05
    m["kBET_acceptance"] = float(acc / neigh.shape[0])
    return m


def compute_metrics_scib(emb, celltype, batch, cluster=None, batch_cluster=None,
                         iso_threshold=None):
    """守卫式上游路径:调用真实 scib.me.*,调用形式逐字取自 scMultiBench 上游脚本
    evaluation_pipelines/scib_metrics/scib_metrics.py。未装 scib 则返回 skipped。

    iso_threshold=None 时按上游 scib_metrics.py:72 的 `num = np.max(batch)+1` 复现,
    即 n_batches+1(scib 文档:"max number of batches per label for label to be
    considered as isolated" → 取 n_batches+1 等于把所有 label 都当 isolated)。
    传整数则覆盖。注意 scib 自身的默认是 None(取 label 出现的最小批次数),
    与上游 scMultiBench 的用法不同 —— 这里跟随上游。"""
    try:
        import scib
        import anndata as ad
        import scanpy as sc
    except ImportError as e:
        return {"status": "skipped",
                "reason": f"未安装 {e.name};上游用法见 https://github.com/PYangLab/scMultiBench",
                "install": "pip install scib  # 需 scanpy/scikit-misc 等,建议独立环境"}
    import numpy as np
    import pandas as pd
    from sklearn.cluster import KMeans

    ct = np.asarray(celltype).astype(str)
    bt = np.asarray(batch).astype(str)
    if cluster is None:
        cluster = KMeans(n_clusters=len(np.unique(ct)), n_init=10,
                         random_state=SEED).fit_predict(emb).astype(str)
    if batch_cluster is None:
        batch_cluster = KMeans(n_clusters=len(np.unique(bt)), n_init=10,
                               random_state=SEED).fit_predict(emb).astype(str)

    A = ad.AnnData(np.asarray(emb))
    A.obsm["X_emb"] = np.asarray(emb)
    A.obs = pd.DataFrame({"celltype": pd.Categorical(ct), "batch": pd.Categorical(bt),
                          "cluster": pd.Categorical(cluster),
                          "batch_cluster": pd.Categorical(batch_cluster)},
                         index=A.obs_names)
    # 上游脚本 scib_metrics.py:101 在 silhouette_batch/graph_connectivity/ilisi 之前跑了
    # sc.pp.neighbors(adata, use_rep="X_emb")(位置在 clisi/ari/nmi/asw 之后)。
    # 这里提前到所有指标之前——等价且更安全:ari/nmi/silhouette 不读邻接图,
    # clisi_graph/ilisi_graph 内部经 recompute_knn() 自建 knn(scib/metrics/lisi.py),
    # 而 graph_connectivity 若缺 adata.uns["neighbors"] 会直接
    # raise KeyError(scib/metrics/graph_connectivity.py:44-45)。
    sc.pp.neighbors(A, use_rep="X_emb")

    # 复现上游 scib_metrics.py:72 的 num = np.max(batch)+1
    num = len(np.unique(bt)) + 1 if iso_threshold is None else iso_threshold
    out = {}

    def _try(key, fn):
        """逐指标守卫:某个指标缺外部依赖(如 kBET 需 R 包 + rpy2)或数据不满足前提时,
        只把该指标标成 unavailable,不让整个评测崩掉。"""
        try:
            out[key] = float(fn())
        except Exception as exc:                       # noqa: BLE001 - 指标级隔离
            out[key] = None
            out[key + "_error"] = f"{type(exc).__name__}: {exc}"

    # ↓ 函数名与关键字参数与上游脚本逐字一致(2026-07-20 复核 scMultiBench 上游脚本 +
    #   theislab/scib 各指标 def 签名,全部对得上)
    _try("clisi", lambda: scib.me.clisi_graph(A, label_key="celltype", type_="embed",
                                              use_rep="X_emb"))
    _try("ari", lambda: scib.me.ari(A, cluster_key="cluster", label_key="celltype"))
    _try("nmi", lambda: scib.me.nmi(A, cluster_key="cluster", label_key="celltype"))
    _try("asw", lambda: scib.me.silhouette(A, label_key="cluster", embed="X_emb"))
    _try("iasw", lambda: scib.me.isolated_labels_asw(A, batch_key="batch", label_key="cluster",
                                                     embed="X_emb", iso_threshold=num))
    _try("if1", lambda: scib.me.isolated_labels_f1(A, batch_key="batch", label_key="celltype",
                                                   cluster_key="cluster", embed="X_emb",
                                                   iso_threshold=num))
    _try("asw_batch", lambda: scib.me.silhouette_batch(A, batch_key="batch",
                                                       label_key="celltype", embed="X_emb"))
    _try("gc", lambda: scib.me.graph_connectivity(A, label_key="celltype"))
    _try("ilisi", lambda: scib.me.ilisi_graph(A, batch_key="batch", type_="embed",
                                              use_rep="X_emb"))
    _try("ari_batch", lambda: 1 - abs(scib.me.ari(A, cluster_key="batch_cluster",
                                                  label_key="batch")))
    _try("nmi_batch", lambda: 1 - abs(scib.me.nmi(A, cluster_key="batch_cluster",
                                                  label_key="batch")))
    # kBET 需要 R 的 kBET 包(scib 内部经 rpy2 调用),纯 pip 装 scib 并不会带上它
    _try("kbet", lambda: scib.me.kBET(A, batch_key="batch", label_key="celltype",
                                      type_="embed", embed="X_emb"))
    ok = [k for k in out if not k.endswith("_error") and out[k] is not None]
    out["status"] = "ok" if ok else "all_metrics_failed"
    out["n_metrics_ok"] = len(ok)
    return out


def aggregate(df):
    """按 scIB 加权汇总:overall = 0.6*mean(bio) + 0.4*mean(batch)。"""
    bio = df[[c for c in BIO_METRICS if c in df.columns]].mean(1)
    bat = df[[c for c in BATCH_METRICS if c in df.columns]].mean(1)
    df = df.copy()
    df["Bio conservation"] = bio
    df["Batch correction"] = bat
    df["Overall"] = W_BIO * bio + W_BATCH * bat
    return df.sort_values("Overall", ascending=False)


# ---------------------------------------------------------------- 出图(无条形图)
def make_figures(df, outdir):
    import matplotlib.pyplot as plt
    import numpy as np
    from pubstyle import set_pub_style, save_fig, pal, CMAP_CONT, NATURE_W2

    set_pub_style(base_size=10)
    files = []
    metrics = [c for c in df.columns if c not in ("Bio conservation", "Batch correction", "Overall")]
    methods = list(df.index)

    # ① 指标热图(scIB 风格总览)
    fig, ax = plt.subplots(figsize=(NATURE_W2, 0.55 * len(methods) + 2.2))
    M = df[metrics].values.astype(float)
    im = ax.imshow(M, cmap=CMAP_CONT, vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(metrics)))
    ax.set_xticklabels(metrics, rotation=40, ha="right")
    ax.set_yticks(range(len(methods)))
    ax.set_yticklabels(methods)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            ax.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center", fontsize=7,
                    color="white" if M[i, j] < 0.55 else "black")
    ax.set_title("scIB-style metrics per integration (higher = better)")
    fig.colorbar(im, ax=ax, shrink=0.7, label="Score")
    for s in ("top", "right"):
        ax.spines[s].set_visible(True)
    f = os.path.join(outdir, "fig1_metric_heatmap")
    save_fig(fig, f); plt.close(fig); files.append(f + ".png")

    # ② bio vs batch 散点(权衡图)
    fig, ax = plt.subplots(figsize=(4.6, 4.2))
    cols = pal(len(methods), "npg")
    for i, mname in enumerate(methods):
        ax.scatter(df["Batch correction"][mname], df["Bio conservation"][mname],
                   s=140, color=cols[i], edgecolor="black", linewidth=0.8, zorder=3,
                   label=mname)
    # 标签上下交错,避免点密集时互相压字
    for i, mname in enumerate(methods):
        dy, va = (9, "bottom") if i % 2 == 0 else (-9, "top")
        ax.annotate(mname, (df["Batch correction"][mname], df["Bio conservation"][mname]),
                    textcoords="offset points", xytext=(0, dy), fontsize=8,
                    ha="center", va=va)
    ax.set_xlabel("Batch correction (mean)")
    ax.set_ylabel("Bio conservation (mean)")
    ax.set_title("Integration trade-off")
    lo = min(df["Batch correction"].min(), df["Bio conservation"].min()) - 0.08
    hi = max(df["Batch correction"].max(), df["Bio conservation"].max()) + 0.08
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    f = os.path.join(outdir, "fig2_bio_vs_batch_scatter")
    save_fig(fig, f); plt.close(fig); files.append(f + ".png")

    # ③ overall 排名 lollipop + 基线参考线(明确"有没有赢过朴素基线")
    fig, ax = plt.subplots(figsize=(5.2, 0.55 * len(methods) + 1.8))
    order = df.sort_values("Overall").index
    y = np.arange(len(order))
    base_key = [m for m in methods if "baseline" in m]
    if base_key:
        ax.axvline(df["Overall"][base_key[0]], color="#999999", ls="--", lw=1,
                   zorder=1, label="naive baseline")
    ax.hlines(y, 0, df["Overall"][order], color="#BBBBBB", lw=1.6, zorder=2)
    ax.scatter(df["Overall"][order], y, s=130, color=pal(3, "npg")[1],
               edgecolor="black", linewidth=0.8, zorder=3)
    for yi, v in zip(y, df["Overall"][order]):
        ax.text(v + 0.035, yi, f"{v:.3f}", va="center", fontsize=8)
    ax.set_yticks(y); ax.set_yticklabels(order)
    ax.set_xlabel("Overall score  (0.6 x bio + 0.4 x batch)")
    ax.set_title("Ranking vs naive baseline")
    ax.set_xlim(0, max(1.0, float(df["Overall"].max()) + 0.12))
    ax.legend(loc="lower right")
    f = os.path.join(outdir, "fig3_overall_lollipop")
    save_fig(fig, f); plt.close(fig); files.append(f + ".png")
    return files


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rna", default=os.path.join(EXAMPLE, "rna_counts.csv"))
    ap.add_argument("--adt", default=os.path.join(EXAMPLE, "adt_counts.csv"))
    ap.add_argument("--meta", default=os.path.join(EXAMPLE, "metadata.csv"),
                    help="须含 cell,celltype,batch 三列")
    ap.add_argument("--emb", action="append", default=[],
                    help="待评测的外部 embedding:--emb 名称=path.csv(行=细胞,与 meta 对齐);可多次")
    ap.add_argument("--k", type=int, default=15, help="kNN 邻居数(LISI/kBET/连通性)")
    ap.add_argument("--use-scib", action="store_true",
                    help="额外走上游 scib 路径(需 pip install scib;未装则跳过并说明)")
    ap.add_argument("--outdir", default=RESULTS)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    os.makedirs(ASSETS, exist_ok=True)

    import numpy as np
    import pandas as pd
    np.random.seed(SEED)

    print("[565] Step 1 读入数据")
    if a.rna == os.path.join(EXAMPLE, "rna_counts.csv"):
        print("       使用 example_data/(synthetic, for demo only)")
    mods = {}
    for name, path in (("rna", a.rna), ("adt", a.adt)):
        if path and os.path.exists(path):
            mods[name] = pd.read_csv(path, comment="#", index_col=0)
    if not mods:
        sys.exit("没有可用的模态输入(--rna / --adt)")
    meta = pd.read_csv(a.meta, comment="#", index_col=0)
    first = next(iter(mods.values()))
    meta = meta.loc[first.index]
    for c in ("celltype", "batch"):
        if c not in meta.columns:
            sys.exit(f"metadata 缺少必需列 '{c}'")
    print(f"       {first.shape[0]} cells · modalities={list(mods)} · "
          f"{meta['celltype'].nunique()} cell types · {meta['batch'].nunique()} batches")

    print("[565] Step 2 构建待评测 embedding(含必跑朴素基线)")
    embs = build_embeddings(mods, meta)
    for spec in a.emb:                                    # 用户自己的整合结果
        if "=" not in spec:
            sys.exit(f"--emb 需写成 名称=path.csv,收到:{spec}")
        name, path = spec.split("=", 1)
        E = pd.read_csv(path, comment="#", index_col=0).loc[meta.index]
        embs[name] = E.values.astype(float)
    print(f"       待评测:{list(embs)}")

    print("[565] Step 3 计算指标(本机 scIB 公式近似实现)")
    rows = {}
    for name, E in embs.items():
        rows[name] = compute_metrics_local(E, meta["celltype"], meta["batch"], k=a.k)
        print(f"       {name}: " + " ".join(f"{k}={v:.3f}" for k, v in rows[name].items()))
    df = aggregate(pd.DataFrame(rows).T)
    csv = os.path.join(a.outdir, "metric.csv")
    df.to_csv(csv, index_label="method")
    print(f"       -> {csv}")

    scib_out = None
    if a.use_scib:
        print("[565] Step 4 上游 scib 路径")
        scib_out = {n: compute_metrics_scib(E, meta["celltype"], meta["batch"])
                    for n, E in embs.items()}
        for n, v in scib_out.items():
            print(f"       {n}: {v.get('status')} {v.get('reason','')}")
        pd.DataFrame(scib_out).T.to_csv(os.path.join(a.outdir, "metric_scib.csv"),
                                        index_label="method")
    else:
        print("[565] Step 4 未请求 --use-scib,跳过上游 scib 路径")

    print("[565] Step 5 出图")
    figs = make_figures(df, a.outdir)
    import shutil
    for f in figs:
        shutil.copy(f, os.path.join(ASSETS, os.path.basename(f)))
    print("       " + " ".join(os.path.basename(f) for f in figs))

    best = df.index[0]
    baseline = [m for m in df.index if "baseline" in m]
    summary = {
        "n_cells": int(first.shape[0]), "modalities": list(mods), "k": a.k,
        "ranking": df["Overall"].round(4).to_dict(),
        "best": best,
        "beats_naive_baseline": (bool(df["Overall"][best] > df["Overall"][baseline[0]])
                                 if baseline else None),
        "metric_impl": "local scIB-formula approximation (KMeans clustering, discrete LISI)",
        "scib_path": (scib_out[list(scib_out)[0]].get("status") if scib_out else "not requested"),
        "seed": SEED,
    }
    with open(os.path.join(a.outdir, "565_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=1, ensure_ascii=False)
    print(f"[565] 完成 · best={best} · 结果在 {a.outdir}")


if __name__ == "__main__":
    main()
