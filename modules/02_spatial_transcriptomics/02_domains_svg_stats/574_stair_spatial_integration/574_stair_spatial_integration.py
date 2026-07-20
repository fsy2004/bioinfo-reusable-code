"""574 · STAIR — multi-slice spatial transcriptomics alignment & integration.

STAIR (Yu & Xie, Genome Biology 2025) is a deep-learning method for spatial
transcriptomic **alignment, integration and de-novo 3D reconstruction**. It uses a
heterogeneous graph attention network (HGAT) over intra-slice (homogeneous) and
inter-slice (heterogeneous) neighbourhoods to learn a shared embedding, from which it
derives spatial domains that are consistent across slices, and reconstructs z-order +
x/y rigid alignment from the attention scores.

This module ships:
  * a RUNNABLE BASELINE (CPU, local deps only) — a three-rung integration ladder
    (PCA / ComBat+PCA / spatial-smoothed+ComBat+PCA) scored with the standard
    integration trade-off (bio conservation ARI-NMI vs batch mixing entropy).
    STAIR is a multi-slice *integration* method, so the honest comparator is a plain
    batch-correction + clustering pipeline. Any STAIR claim must beat this floor.
  * a GUARDED STAIR path (--run-stair) — STAIR is not installed locally and needs a
    GPU. Signatures below were read from the upstream source, not guessed.

Verified upstream sources (本地克隆逐行核对,2026-07-21 复核):
  repo      https://github.com/yuyuanyuana/STAIR  (STAIR-tools 1.3.1, MIT, setup.py)
  API src   STAIR/emb_alignment.py · STAIR/utils.py · STAIR/loc_alignment.py
            STAIR/loc_prediction.py · STAIR/embedding/dataset_ae.py
            (逐个符号的文件:行号见 run_stair() docstring)
  tutorial  https://stair-tutorial.readthedocs.io/en/latest/STAIR-Tutorial.html
            (链接出自上游 README「Tutorial」节,本模块未抓取其内容)
  paper     Genome Biol 2025;26(1):427 · PMID 41398698 · doi:10.1186/s13059-025-03895-x
            (esummary 核实:Yu Y, Xie Z;标题/卷期页/DOI 全部对上)
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

# 框架统一出图样式
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..", "..", "_framework")))
from pubstyle import set_pub_style, pal, save_fig, NATURE_W1, NATURE_W2  # noqa: E402

SEED = 2026


# ---------------------------------------------------------------------------
# 合成示例数据:3 张切片、4 个同心空间域、切片特异批次效应 + 中心平移/旋转
# ---------------------------------------------------------------------------
def make_synthetic(out_dir: str, n_side: int = 16, n_slices: int = 3,
                   n_genes: int = 200, n_domains: int = 4, seed: int = SEED,
                   marker_fc: float = 0.9, depth: float = 2.0, p_drop: float = 0.3):
    """写出 synthetic 表达矩阵与 meta。真实域标签保留,用作整合质量的金标准。

    ★ 信噪比是刻意调低的(弱 marker + 浅测序 + dropout)。若把信号调干净,
    朴素 PCA 会直接拿到 ARI=1.0,阶梯失去分辨力、空间平滑反而显得有害 ——
    那样的"基准"什么都证明不了。这里的参数让单点表达本身是模糊的,
    只有借助邻域信息才能稳定判域,这才贴近真实 ST 数据的处境。
    """
    import numpy as np
    import pandas as pd

    rng = np.random.default_rng(seed)
    os.makedirs(out_dir, exist_ok=True)

    gx, gy = np.meshgrid(np.arange(n_side), np.arange(n_side))
    gx, gy = gx.ravel().astype(float), gy.ravel().astype(float)
    c = (n_side - 1) / 2.0
    rmax = np.hypot(c, c)

    # 每个域一组 marker 基因
    gene_domain = rng.integers(0, n_domains, n_genes)
    base = rng.uniform(0.3, 1.2, n_genes)

    rows, meta = [], []
    for s in range(n_slices):
        # 切片间的刚体错位:平移 + 旋转 —— 这正是 STAIR 要对齐的东西
        shift = rng.normal(0, 1.2, 2)
        th = rng.uniform(-0.25, 0.25)
        xr = (gx - c) * np.cos(th) - (gy - c) * np.sin(th) + c + shift[0]
        yr = (gx - c) * np.sin(th) + (gy - c) * np.cos(th) + c + shift[1]

        # 真实域 = 以(错位后的)切片自身中心为准的同心环
        r = np.hypot(xr - c - shift[0], yr - c - shift[1]) / rmax
        dom = np.clip((r * n_domains).astype(int), 0, n_domains - 1)

        # 域程序 + 切片特异批次效应(基因水平乘性 + 文库大小)
        prof = base[None, :] * (1.0 + marker_fc * (gene_domain[None, :] == dom[:, None]))
        batch_gene = np.exp(rng.normal(0, 0.55, n_genes))       # 基因×批次 乘性偏移
        libsize = np.exp(rng.normal(0, 0.20, len(gx)))          # 每点文库差异
        lam = prof * batch_gene[None, :] * libsize[:, None] * depth
        counts = rng.poisson(np.clip(lam, 0.01, None))
        counts = counts * (rng.random(counts.shape) > p_drop)   # dropout(ST 常态)

        rows.append(counts)
        for i in range(len(gx)):
            meta.append((f"S{s}_spot{i:04d}", f"slice{s}", xr[i], yr[i], int(dom[i])))

    X = np.vstack(rows)
    meta = pd.DataFrame(meta, columns=["spot", "slice", "x", "y", "true_domain"])
    expr = pd.DataFrame(X, index=meta["spot"], columns=[f"G{i:03d}" for i in range(n_genes)])

    hdr = "# synthetic, for demo only -- module 574 STAIR baseline\n"
    for path, df, idx in ((os.path.join(out_dir, "slices_expression.csv"), expr, True),
                          (os.path.join(out_dir, "slices_meta.csv"), meta, False)):
        with open(path, "w", newline="") as fh:
            fh.write(hdr)
            df.to_csv(fh, index=idx)
    return expr, meta


def load_example(out_dir: str):
    import pandas as pd
    ex = os.path.join(out_dir, "slices_expression.csv")
    me = os.path.join(out_dir, "slices_meta.csv")
    if not (os.path.exists(ex) and os.path.exists(me)):
        print("[574] example_data missing, generating synthetic slices")
        return make_synthetic(out_dir)
    expr = pd.read_csv(ex, comment="#", index_col=0)
    meta = pd.read_csv(me, comment="#")
    return expr, meta


# ---------------------------------------------------------------------------
# 基线:三级整合阶梯(全部只用本机已装依赖)
# ---------------------------------------------------------------------------
def _preprocess(expr, meta):
    """标准 scanpy 预处理 → AnnData(normalize/log1p/scale)。"""
    import anndata as ad
    import numpy as np
    import scanpy as sc

    A = ad.AnnData(expr.values.astype("float32"))
    A.obs_names = list(expr.index)
    A.var_names = list(expr.columns)
    A.obs["slice"] = meta["slice"].values
    # true_domain 是可选的:真实数据通常没有金标准,此时只报批次混合,不报 ARI/NMI
    if "true_domain" in meta.columns:
        A.obs["true_domain"] = meta["true_domain"].values.astype(str)
    else:
        A.obs["true_domain"] = "NA"
    A.obsm["spatial"] = meta[["x", "y"]].values.astype("float32")
    # ★ 原始 count 必须留一份:STAIR 的 MyDataset 把 adata.X 当 count 用
    #   (自己再做 normalize_total + log1p,且 nb/zinb 似然、seurat_v3 HVG 都要 count),
    #   把 log 过的矩阵喂进去会二次 log1p 并让似然失配。
    A.layers["counts"] = A.X.copy()
    sc.pp.normalize_total(A, target_sum=1e4)
    sc.pp.log1p(A)
    A.raw = A
    return A


def _spatial_smooth(A, k: int = 6):
    """切片内空间 kNN 均值平滑 —— STAIR 同质图(intra-slice)的朴素替身。

    这不是 HGAT:没有注意力、没有跨切片异质边。它存在的意义是把"空间信息"
    这一项单独拆出来,好判断 STAIR 的增益到底来自空间图还是来自注意力。
    """
    import numpy as np
    from sklearn.neighbors import NearestNeighbors

    Xs = np.array(A.X, dtype=float, copy=True)
    out = np.empty_like(Xs)
    for sl in A.obs["slice"].unique():
        m = (A.obs["slice"] == sl).values
        coords = A.obsm["spatial"][m]
        nn = NearestNeighbors(n_neighbors=min(k + 1, m.sum())).fit(coords)
        _, idx = nn.kneighbors(coords)
        out[m] = Xs[m][idx].mean(axis=1)
    B = A.copy()
    B.X = out.astype("float32")
    return B


def _embed(A, n_pcs: int = 20, combat: bool = False, smooth: bool = False):
    import numpy as np
    import scanpy as sc

    B = _spatial_smooth(A) if smooth else A.copy()
    if combat:
        sc.pp.combat(B, key="slice")           # 非空间批次校正(scanpy 自带,无需额外包)
    sc.pp.scale(B, max_value=10)
    sc.tl.pca(B, n_comps=n_pcs, svd_solver="arpack", random_state=SEED)
    return B.obsm["X_pca"]


def _batch_entropy(emb, batches, k: int = 30, seed: int = SEED):
    """局部邻域批次熵(0=完全分离, 1=完全混合)。iLISI 同精神的轻量实现。"""
    import numpy as np
    from sklearn.neighbors import NearestNeighbors

    b = np.asarray(batches)
    labs = np.unique(b)
    code = np.searchsorted(labs, b)
    nn = NearestNeighbors(n_neighbors=min(k + 1, len(b))).fit(emb)
    _, idx = nn.kneighbors(emb)
    idx = idx[:, 1:]                                   # 去掉自身
    ent = []
    for row in code[idx]:
        p = np.bincount(row, minlength=len(labs)) / row.size
        p = p[p > 0]
        ent.append(-(p * np.log(p)).sum())
    return float(np.mean(ent) / np.log(len(labs)))


def run_baseline(A, n_domains: int, outdir: str):
    """三级阶梯 + scIB 式双轴评分(生物保真 vs 批次混合)。"""
    import numpy as np
    import pandas as pd
    from sklearn.cluster import KMeans
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

    truth = A.obs["true_domain"].values
    batches = A.obs["slice"].values
    has_truth = len(np.unique(truth)) > 1          # 无金标准时 ARI/NMI 无意义
    if not has_truth:
        print("       no 'true_domain' column -> reporting batch mixing only")

    ladder = {
        "PCA (no correction)":      dict(combat=False, smooth=False),
        "ComBat + PCA":             dict(combat=True,  smooth=False),
        "Spatial + ComBat + PCA":   dict(combat=True,  smooth=True),
    }

    embs, recs = {}, []
    for name, kw in ladder.items():
        print(f"       fitting: {name}")
        E = _embed(A, **kw)
        lab = KMeans(n_clusters=n_domains, n_init=20, random_state=SEED).fit_predict(E)
        embs[name] = (E, lab)
        recs.append({
            "method": name,
            "ARI": adjusted_rand_score(truth, lab) if has_truth else np.nan,
            "NMI": normalized_mutual_info_score(truth, lab) if has_truth else np.nan,
            "batch_entropy": _batch_entropy(E, batches),
        })
    df = pd.DataFrame(recs)
    # 综合分:与 scIB 一致的加权(生物保真 0.6 / 批次去除 0.4);无金标准时退化为纯批次分
    if has_truth:
        df["overall"] = 0.6 * df[["ARI", "NMI"]].mean(axis=1) + 0.4 * df["batch_entropy"]
    else:
        df["overall"] = df["batch_entropy"]
    df.attrs["has_truth"] = has_truth
    df.to_csv(os.path.join(outdir, "574_baseline_metrics.csv"), index=False)
    return df, embs


# ---------------------------------------------------------------------------
# 出图(顶刊风格;无条形图)
# ---------------------------------------------------------------------------
def fig_spatial_domains(A, embs, figdir):
    import matplotlib.pyplot as plt
    import numpy as np

    slices = list(dict.fromkeys(A.obs["slice"]))
    panels = []
    if len(np.unique(A.obs["true_domain"])) > 1:        # 无金标准就不画这一行
        panels.append(("Ground truth", A.obs["true_domain"].astype("category").cat.codes.values))
    panels += [(n, lab) for n, (_, lab) in embs.items()]

    nr, nc = len(panels), len(slices)
    fig, axes = plt.subplots(nr, nc, figsize=(NATURE_W2, 2.0 * nr), squeeze=False)
    cols = pal(max(4, len(np.unique(panels[0][1]))), "okabe_ito")
    xy = A.obsm["spatial"]
    for i, (name, lab) in enumerate(panels):
        for j, sl in enumerate(slices):
            ax = axes[i][j]
            m = (A.obs["slice"] == sl).values
            ax.scatter(xy[m, 0], xy[m, 1], c=[cols[v % len(cols)] for v in lab[m]],
                       s=9, linewidths=0)
            ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
            for sp in ax.spines.values():
                sp.set_visible(False)
            if i == 0:
                ax.set_title(sl, fontsize=9)
            if j == 0:
                ax.set_ylabel(name, fontsize=8, rotation=0, ha="right", va="center")
    fig.suptitle("Cross-slice spatial domains (synthetic benchmark)", fontsize=10)
    save_fig(fig, os.path.join(figdir, "fig1_spatial_domains"))
    plt.close(fig)


def fig_tradeoff(df, figdir):
    """整合权衡散点:x = 批次混合, y = 生物保真。右上角最优。"""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(NATURE_W1 + 1.1, NATURE_W1))
    cols = pal(len(df), "npg")
    bio = df[["ARI", "NMI"]].mean(axis=1)
    ax.plot(df["batch_entropy"], bio, "-", color="#BBBBBB", lw=1, zorder=1)
    for i, r in df.iterrows():
        ax.scatter(r["batch_entropy"], bio[i], s=110, color=cols[i],
                   edgecolor="black", linewidth=0.8, zorder=3, label=r["method"])
    ax.set_xlabel("Batch mixing entropy  (higher = better mixed)")
    ax.set_ylabel("Bio conservation  (mean of ARI, NMI)")
    ax.set_title("Integration trade-off", fontsize=10)
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=7)
    save_fig(fig, os.path.join(figdir, "fig2_integration_tradeoff"))
    plt.close(fig)


def fig_lollipop(df, figdir):
    """棒棒糖图逐指标对比(库规矩:不用条形图)。"""
    import matplotlib.pyplot as plt
    import numpy as np

    metrics = ["ARI", "NMI", "batch_entropy", "overall"]
    fig, axes = plt.subplots(1, len(metrics), figsize=(NATURE_W2, 1.9), sharey=True)
    cols = pal(len(df), "npg")
    y = np.arange(len(df))[::-1]
    for ax, met in zip(axes, metrics):
        ax.hlines(y, 0, df[met], color="#CCCCCC", lw=1.4, zorder=1)
        ax.scatter(df[met], y, s=70, color=cols, edgecolor="black", linewidth=0.7, zorder=3)
        ax.set_title(met, fontsize=9)
        ax.set_xlim(0, max(1.0, df[met].max() * 1.15))
        ax.set_ylim(-0.6, len(df) - 0.4)       # 收紧纵向留白,避免点被拉散
    axes[0].set_yticks(y); axes[0].set_yticklabels(df["method"], fontsize=7)
    fig.suptitle("Baseline integration ladder", fontsize=10, y=1.12)  # y 抬高,免与子图标题相撞
    save_fig(fig, os.path.join(figdir, "fig3_metric_lollipop"))
    plt.close(fig)


def fig_umap(A, embs, figdir):
    import matplotlib.pyplot as plt
    import numpy as np
    import umap

    slices = list(dict.fromkeys(A.obs["slice"]))
    cols = pal(max(3, len(slices)), "okabe_ito")   # 色盲安全:避开 lancet 的红绿对
    fig, axes = plt.subplots(1, len(embs), figsize=(NATURE_W2, 2.5), squeeze=False)
    for ax, (name, (E, _)) in zip(axes[0], embs.items()):
        U = umap.UMAP(n_neighbors=15, min_dist=0.2, random_state=SEED).fit_transform(E)
        for k, sl in enumerate(slices):
            m = (A.obs["slice"] == sl).values
            ax.scatter(U[m, 0], U[m, 1], s=6, color=cols[k], linewidths=0, label=sl, alpha=0.8)
        ax.set_title(name, fontsize=8)
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():           # 无刻度的坐标框是视觉噪声
            sp.set_visible(False)
    axes[0][0].legend(fontsize=6, markerscale=1.6, loc="best")
    fig.suptitle("Embedding coloured by slice (mixed = batch effect removed)",
                 fontsize=10, y=1.04)
    save_fig(fig, os.path.join(figdir, "fig4_umap_by_slice"))
    plt.close(fig)


# ---------------------------------------------------------------------------
# STAIR 守卫式封装 —— 签名取自上游源码,未在本机执行过
# ---------------------------------------------------------------------------
def run_stair(adata, result_path: str, batch_key: str = "batch",
              slice_order=None, n_domains: int = 9, hvg: int = 3000,
              clustering: str = "kmeans"):
    """STAIR 真实调用序列。需 `pip install` 官方 wheel + CUDA GPU。

    ★ 签名逐个核对自**本地克隆的上游源码**(非教程页面),行号为 STAIR @ main:
      STAIR/emb_alignment.py :
        L23  class Emb_Align(adata, batch_key=None, hvg=False, n_hidden=128,
                 n_latent=32, dropout_rate=0.2, likelihood='nb', device=None,
                 num_workers=4, result_path=None, make_log=True)
        L114 .prepare(count_key=None, lib_size='explog', normalize=True, scale=False)
        L148 .preprocess(lr=0.001, weight_decay=0, epoch_ae=100, batch_size=128, plot=False)
        L226 .latent(batch_size=10000, return_data=False)      → obsm['latent'] (L285)
             ↑ 必须在 prepare_hgat 之前调用:异质图构建读的就是它
               (STAIR/embedding/dataset_hgat.py L107-108 `adata_tmp.obsm['latent']`
                与 `adata_tmp.obsm[spatial_key]`)
        L291 .prepare_hgat(slice_key=None, slice_order=None, spatial_key='spatial',
                 n_neigh_hom=10, c_neigh_het=0.9, kernal_thresh=0.)
        L334 .train_hgat(gamma=0.8, epoch_hgat=150, re_weight=1., si_weight=0.,
                 lr=0.001, weight_decay=0., negative_slope=0.2, dropout_hom=0.5,
                 dropout_het=0.5, mini_batch=False, batch_size=256, batches=100,
                 num_hops=2, plot=False)
        L533 .predict_hgat(mini_batch=False, batches=100, num_hops=2,
                 get_attention=False) → (adata, atte_);  obsm['STAIR'] 被写入
      STAIR/utils.py :
        L85  cluster_func(adata, clustering, use_rep, res=1, cluster_num=None,
                 key_add='cluster')   # 'louvain'/'leiden'/'kmeans'/'mclust'
        L62  mclust_R(...)  ← clustering='mclust' 走这里,**需 rpy2 + R 包 mclust**;
             本机无 rpy2,故本封装默认用纯 sklearn 的 'kmeans' 分支。
      3D 重建路径(同样已核对源码,非仅教程):
        STAIR/loc_prediction.py  L13  sort_slices(atte, start=None, return_tree=False)
        STAIR/loc_alignment.py   L12  class Loc_Align(adata, batch_key,
                 batch_order=None, make_log=True, result_path='.')
          L65  .init_align(emb_key, spatial_key='spatial', num_mnn=1,
                   init_align_key='transform_init', use_scale=False, return_result=False)
          L109 .detect_fine_points(slice_boundary=True, domain_boundary=True,
                   domain_key='layer_cluster', num_domains=1, sep_sort=True,
                   alpha=70, return_result=False)
          L201 .fine_align(fine_align_key='transform_fine', max_iterations=20,
                   tolerance=1e-10, return_result=False) → adata_aligned

    ★ 本函数未在本机执行过(无 wheel、无 GPU),仅保证签名与源码一致。
    """
    try:
        from STAIR.emb_alignment import Emb_Align
        from STAIR.utils import cluster_func
    except ImportError as e:
        return {"status": "skipped",
                "reason": f"STAIR not installed ({e.name}); "
                          "pip install https://github.com/yuyuanyuana/STAIR/releases/"
                          "download/1.3.1/STAIR_tools-1.3.1-py3-none-any.whl"}
    try:
        import torch
        if not torch.cuda.is_available():
            return {"status": "skipped", "reason": "no CUDA GPU; STAIR HGAT training needs one"}
        dev = "cuda:0"
    except ImportError:
        return {"status": "skipped", "reason": "torch not installed"}

    # ★ STAIR 要的是 count,不是 log 后的矩阵。
    #   源码 STAIR/embedding/dataset_ae.py L25-27: count_key=None 时直接
    #   `adata.layers['counts'] = adata.X.copy()`,随后 L63/L66 自己做
    #   normalize_total + log1p;L68-74 的 HVG 又是 flavor='seurat_v3'(要 raw count)。
    #   所以这里必须把 X 换回 counts 层,否则二次 log1p + nb/zinb 似然失配。
    if "counts" not in adata.layers:
        return {"status": "skipped",
                "reason": "no 'counts' layer; STAIR needs raw counts, not log-normalised X"}
    ad_st = adata.copy()
    ad_st.X = ad_st.layers["counts"].copy()

    # HVG 数不能超过基因数(seurat_v3 会报错);合成示例只有 200 基因。
    hvg_use = hvg if (hvg and ad_st.shape[1] > hvg) else False

    emb = Emb_Align(ad_st, batch_key=batch_key, hvg=hvg_use, likelihood="zinb",
                    result_path=result_path, device=dev)
    emb.prepare()                       # count_key=None → 用 X 当 count(已是 count)
    emb.preprocess()
    emb.latent()
    emb.prepare_hgat(slice_key=batch_key, slice_order=slice_order,
                     n_neigh_hom=10, c_neigh_het=0.9)
    emb.train_hgat(gamma=0.8)
    ad_st, atte = emb.predict_hgat()    # 源码返回 (adata, atte_);atte 供 sort_slices 定 z 序
    ad_st = cluster_func(ad_st, clustering=clustering, use_rep="STAIR",
                         cluster_num=n_domains, key_add="STAIR")
    return {"status": "ran", "obsm_key": "STAIR", "domain_key": "STAIR",
            "clustering": clustering, "hvg": hvg_use,
            "n_slices_in_attention": int(atte.shape[0])}


# ---------------------------------------------------------------------------
def main():
    import pandas as pd

    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--expr", default=os.path.join(EXAMPLE, "slices_expression.csv"))
    ap.add_argument("--meta", default=os.path.join(EXAMPLE, "slices_meta.csv"))
    ap.add_argument("--outdir", default=RESULTS)
    ap.add_argument("--figdir", default=ASSETS)
    ap.add_argument("--n-domains", type=int, default=4)
    ap.add_argument("--run-stair", action="store_true",
                    help="attempt the real STAIR path (needs the wheel + a CUDA GPU)")
    a = ap.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    os.makedirs(a.figdir, exist_ok=True)
    set_pub_style(base_size=9)

    print("[574] Step 1 · load slices")
    if a.expr == os.path.join(EXAMPLE, "slices_expression.csv"):
        expr, meta = load_example(EXAMPLE)
    else:
        expr = pd.read_csv(a.expr, comment="#", index_col=0)
        meta = pd.read_csv(a.meta, comment="#")
    print(f"       {expr.shape[0]} spots x {expr.shape[1]} genes, "
          f"{meta['slice'].nunique()} slices")

    print("[574] Step 2 · preprocess")
    A = _preprocess(expr, meta)

    print("[574] Step 3 · baseline integration ladder")
    df, embs = run_baseline(A, a.n_domains, a.outdir)
    print(df.to_string(index=False, float_format=lambda v: f"{v:.3f}"))

    print("[574] Step 4 · figures")
    fig_spatial_domains(A, embs, a.figdir)
    if df.attrs.get("has_truth", True):
        fig_tradeoff(df, a.figdir)      # 权衡图的 y 轴就是 ARI/NMI,无金标准时画不出来
    else:
        print("       trade-off panel skipped (no ground-truth domains)")
    fig_lollipop(df, a.figdir)
    try:
        fig_umap(A, embs, a.figdir)
    except Exception as e:                     # UMAP 是锦上添花,失败不拖垮主流程
        print(f"       UMAP panel skipped: {type(e).__name__}: {e}")

    print("[574] Step 5 · STAIR path")
    if a.run_stair:
        st = run_stair(A, result_path=a.outdir, batch_key="slice",
                       slice_order=list(dict.fromkeys(A.obs["slice"])),
                       n_domains=a.n_domains)
    else:
        st = {"status": "not requested", "reason": "pass --run-stair (needs wheel + GPU)"}
    for k, v in st.items():
        print(f"       {k}: {v}")

    with open(os.path.join(a.outdir, "574_summary.json"), "w") as fh:
        json.dump({"baseline": df.to_dict("records"), "stair": st}, fh, indent=1, default=str)
    print(f"[574] done -> {a.outdir}")


if __name__ == "__main__":
    main()
