"""563 · CONCORD — self-supervised contrastive integration / denoising / dim-reduction.

CONCORD (Zhu et al., *Nature Biotechnology* 2026; PMID 41491253,
doi:10.1038/s41587-025-02950-z) learns a single latent space that simultaneously denoises,
reduces dimension and removes batch effect, using dataset-aware + neighbourhood-aware
contrastive sampling. Per the published abstract, its representations preserve "both local
geometric relationships and global topological structures" without external supervision.

This module ships:
  * a RUNNABLE BASELINE bank that needs nothing beyond scanpy/sklearn — unintegrated PCA,
    per-batch mean-centred PCA, and ComBat+PCA — plus a common evaluation harness;
  * a GUARDED CONCORD path (`--run-concord`) that only fires when `concord-sc` is importable.

The baselines exist because an integration method that claims "better" must be shown against
a naive floor on the same metrics. Never report CONCORD alone.

Upstream API grounded in a local clone of Gartner-Lab/Concord (concord-sc v1.0.13,
src/concord/__init__.py:1), re-verified line-by-line 2026-07-21. Do not invent beyond these:
  Concord               src/concord/concord.py:44   (class), :59 (__init__ signature)
  Concord.default_params src/concord/concord.py:94-143  (every **kwargs key must appear here;
                        setup_config():171 raises ValueError on any unknown key)
  Concord.fit_transform src/concord/concord.py:613
  -> obsm write          src/concord/concord.py:590  (_add_results_to_adata)
  ccd.ul.select_features src/concord/utils/feature_selector.py:182
  ccd.ul.run_umap        src/concord/utils/dim_reduction.py:5
  ccd.pl.plot_embedding  src/concord/plotting/pl_embedding.py:29
  namespaces             src/concord/__init__.py -> `from . import ml, pl, ul, bm, sm`
  README/docs            upstream README.md; https://qinzhu.github.io/Concord_documentation/
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import warnings

warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")
ASSETS = os.path.join(HERE, "assets")
EXAMPLE = os.path.join(HERE, "example_data")
FRAMEWORK = os.path.abspath(os.path.join(HERE, "..", "..", "..", "_framework"))
sys.path.insert(0, FRAMEWORK)

from pubstyle import (  # noqa: E402
    CMAP_DIVERGE, NATURE_W1, NATURE_W2, pal, save_fig, set_pub_style,
)

SEED = 0


# ---------------------------------------------------------------------------
# 合成示例数据:一条连续的分化轨迹 + 3 个细胞类型 + 3 个批次
# 关键:真实低维坐标被保存下来,才能量化"全局几何是否被保住"
# ---------------------------------------------------------------------------
def make_synthetic(n_per_batch: int = 260, n_genes: int = 400, n_batches: int = 3,
                   seed: int = SEED):
    import numpy as np
    import pandas as pd

    rng = np.random.default_rng(seed)
    n_cells = n_per_batch * n_batches

    # 真实生物结构:一条 1D 轨迹(pseudotime)+ 沿轨迹切出的 3 个细胞类型
    t = rng.uniform(0, 1, n_cells)
    ct = np.digitize(t, [1 / 3, 2 / 3])
    ct_name = np.array(["Progenitor", "Intermediate", "Mature"])[ct]

    # 基因程序:每个基因是 pseudotime 的一个平滑高斯响应(峰值位置随机)
    peak = rng.uniform(-0.15, 1.15, n_genes)
    width = rng.uniform(0.12, 0.45, n_genes)
    amp = rng.uniform(1.0, 5.0, n_genes)
    signal = amp[None, :] * np.exp(-((t[:, None] - peak[None, :]) ** 2) / (2 * width[None, :] ** 2))

    # 批次效应:基因水平的加性偏移 + 批次特异的测序深度缩放(真实数据里两者都有)
    batch = np.repeat(np.arange(n_batches), n_per_batch)
    rng.shuffle(batch)
    shift = rng.normal(0, 1.1, (n_batches, n_genes))
    depth = np.array([1.0, 0.62, 1.55])[:n_batches]
    mu = np.clip((signal + shift[batch]) * depth[batch][:, None], 0.02, None)

    counts = rng.poisson(mu).astype("float32")

    import anndata as ad
    A = ad.AnnData(counts)
    A.var_names = [f"G{i:04d}" for i in range(n_genes)]
    A.obs_names = [f"C{i:05d}" for i in range(n_cells)]
    A.obs["batch"] = pd.Categorical([f"batch{b + 1}" for b in batch])
    A.obs["cell_type"] = pd.Categorical(ct_name)
    A.obs["true_time"] = t
    # 真实几何参考:无批次效应的干净信号空间
    A.obsm["X_true"] = signal.astype("float32")
    return A


def write_example(path: str, adata, n_ref: int = 20) -> None:
    """把合成数据落成 CSV(表达矩阵 + 元信息 + 真实几何参考)。

    第三个文件 `true_geometry.csv` 是「无批次效应的干净信号空间」的前 n_ref 个主成分,
    只有合成数据才有。它是评估 trustworthiness / global geometry 的参照系;没有它,
    未整合 PCA 会自己当参照系而拿到平凡的 1.000,比较就废了。
    """
    import pandas as pd
    from sklearn.decomposition import PCA

    os.makedirs(path, exist_ok=True)
    expr = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names)
    meta = adata.obs[["batch", "cell_type", "true_time"]].copy()
    ref = PCA(n_components=n_ref, random_state=SEED).fit_transform(adata.obsm["X_true"])
    refdf = pd.DataFrame(ref, index=adata.obs_names,
                         columns=[f"TruePC{i + 1}" for i in range(n_ref)])
    header = "# synthetic, for demo only -- module 563 CONCORD\n"
    for fn, df in (("counts.csv", expr), ("cell_meta.csv", meta),
                   ("true_geometry.csv", refdf)):
        with open(os.path.join(path, fn), "w", newline="") as fh:
            fh.write(header)
            df.to_csv(fh)


def read_example(path: str):
    """读回 example_data/ 的 CSV。true_geometry.csv 可选,缺失即当作真实数据处理。"""
    import anndata as ad
    import pandas as pd

    expr = pd.read_csv(os.path.join(path, "counts.csv"), comment="#", index_col=0)
    meta = pd.read_csv(os.path.join(path, "cell_meta.csv"), comment="#", index_col=0)
    A = ad.AnnData(expr.values.astype("float32"))
    A.obs_names = list(expr.index)
    A.var_names = list(expr.columns)
    for c in meta.columns:
        A.obs[c] = meta[c].values
    for c in ("batch", "cell_type"):
        if c in A.obs:
            A.obs[c] = A.obs[c].astype("category")
    ref_fp = os.path.join(path, "true_geometry.csv")
    if os.path.exists(ref_fp):
        ref = pd.read_csv(ref_fp, comment="#", index_col=0)
        A.obsm["X_true"] = ref.loc[A.obs_names].values.astype("float32")
    return A


# ---------------------------------------------------------------------------
# 预处理(所有方法共用同一份输入,否则比较不公平)
# ---------------------------------------------------------------------------
def preprocess(adata, n_hvg: int = 300):
    import scanpy as sc

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers["lognorm"] = adata.X.copy()
    try:
        sc.pp.highly_variable_genes(adata, n_top_genes=min(n_hvg, adata.n_vars - 1))
        adata.uns["hvg_list"] = list(adata.var_names[adata.var["highly_variable"]])
    except Exception:
        adata.uns["hvg_list"] = list(adata.var_names)
    return adata


# ---------------------------------------------------------------------------
# 基线库:本机依赖即可跑完的朴素对照
# ---------------------------------------------------------------------------
def embed_pca(adata, n_comps: int = 30):
    """无整合对照:直接 PCA。批次效应会原样留在里面,这是下限。"""
    import numpy as np
    from sklearn.decomposition import PCA

    X = np.asarray(adata.layers["lognorm"])
    return PCA(n_components=n_comps, random_state=SEED).fit_transform(X)


def embed_batch_center(adata, batch_key: str = "batch", n_comps: int = 30):
    """最朴素的整合:按批次做基因级均值中心化后再 PCA(等价于线性去批次截距项)。"""
    import numpy as np
    from sklearn.decomposition import PCA

    X = np.asarray(adata.layers["lognorm"]).copy()
    b = np.asarray(adata.obs[batch_key])
    for lev in np.unique(b):
        m = b == lev
        X[m] -= X[m].mean(axis=0, keepdims=True)
    return PCA(n_components=n_comps, random_state=SEED).fit_transform(X)


def embed_combat(adata, batch_key: str = "batch", n_comps: int = 30):
    """ComBat(Johnson et al. Biostatistics 2007),scanpy 自带,不需要额外装包。"""
    import numpy as np
    import scanpy as sc
    from sklearn.decomposition import PCA

    tmp = adata.copy()
    tmp.X = np.asarray(adata.layers["lognorm"]).copy()
    sc.pp.combat(tmp, key=batch_key)
    return PCA(n_components=n_comps, random_state=SEED).fit_transform(np.asarray(tmp.X))


# ---------------------------------------------------------------------------
# CONCORD 路径:守卫式封装。装不上就诚实退出,绝不静默降级冒充 CONCORD 结果。
# ---------------------------------------------------------------------------
def embed_concord(adata, batch_key: str = "batch", latent_dim: int = 30,
                  n_epochs: int = 15, output_key: str = "Concord"):
    """调用真实 CONCORD。每个符号与形参都对得上本地克隆的上游源码(v1.0.13)。

    逐条核实(2026-07-21,读 C:\\Users\\fsy\\Desktop\\upstream-sources\\563_Concord):
      ccd.Concord(adata, save_dir='save/', copy_adata=False, verbose=False, **kwargs)
          concord.py:59。前四个是显式形参,其余一律走 **kwargs;
          setup_config(concord.py:171)会把不在 default_params 里的键直接 ValueError
          (concord.py:187),所以下面每个键都必须在 default_params(concord.py:94-143)中存在:
            input_feature (:96, 默认 None)   domain_key   (:109, 默认 None)
            latent_dim    (:104, 默认 100)   n_epochs     (:100, 默认 15)
            seed          (:95,  默认 0)     device       (:143, 默认 cuda:0 if available)
          相关但本函数未设的:batch_size(:99,256)/ lr(:101,1e-2)/
            clr_temperature(:116,0.4)/ p_intra_knn(:133,0.0)/ p_intra_domain(:134,1.0)/
            encoder_dims(:105,[1000])/ decoder_dims(:106,[1000])/ preload_dense(:139,False)
      .fit_transform(output_key="Concord", return_decoded=False, decoder_domain=None,
                     return_class=True, return_class_prob=True, save_model=True)
          concord.py:613;经 _add_results_to_adata(concord.py:585)在 :590 处执行
          `adata_to_update.obsm[output_key] = embeddings`,写的是 _adata_original(:649),
          所以我们直接读传进去的那个 adata.obsm[output_key] 是对的。
          注意 fit_transform 无 return 语句(返回 None),必须回读 obsm,不能接返回值。
      输入格式:default_params 的 normalize_total / log1p 均为 False,上游注释写明
          "default adata.X should be normalized"(concord.py:97),_check_input_format
          (:845)只对 raw counts 发 warning。本模块传入的 adata.X 已 normalize+log1p。
      faiss 缺失不致命:knn.py:80 会 warn 后回落 sklearn。
    另外两个上游辅助函数本模块未调用(出图走本库 pubstyle),但接口同样已核对:
      ccd.ul.run_umap(source_key='encoded', result_key='encoded_UMAP', n_components=2,
          n_pc=None, n_neighbors=30, min_dist=0.1, metric='cosine', ...)  dim_reduction.py:5
      ccd.ul.select_features(adata, n_top_features=2000, flavor='seurat_v3', ...)
          feature_selector.py:182

    超参的推荐取值随数据规模变化,上游教程为准;此处只固定我们直接读到的形参名。
    """
    try:
        import concord as ccd
    except ImportError as e:
        return None, {"status": "skipped",
                      "reason": f"concord not importable ({e.name}); install: pip install concord-sc"}

    import numpy as np
    import torch

    feats = adata.uns.get("hvg_list") or list(adata.var_names)
    dev = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    model = ccd.Concord(adata, save_dir=os.path.join(RESULTS, "concord_save"),
                        input_feature=feats, domain_key=batch_key,
                        latent_dim=latent_dim, n_epochs=n_epochs,
                        seed=SEED, device=dev, verbose=False)
    model.fit_transform(output_key=output_key)
    if output_key not in adata.obsm:
        return None, {"status": "error",
                      "reason": f"fit_transform did not populate adata.obsm['{output_key}']"}
    return np.asarray(adata.obsm[output_key]), {
        "status": "ok",
        "device": str(dev),
        "latent_dim": latent_dim,
        "n_epochs": n_epochs,
        "n_input_features": len(feats),
    }


# ---------------------------------------------------------------------------
# 评估:批次混合 vs 生物保真 vs 几何保真(三者都要看,单看一个必然被 gaming)
# ---------------------------------------------------------------------------
def knn_batch_entropy(emb, batch, k: int = 30):
    """每个细胞 kNN 邻域内批次标签的归一化香农熵。1=完全混合,0=完全隔离。

    思路同 Azizi et al. Cell 2018 的 batch mixing entropy;这里做了按全局批次占比
    归一化,避免批次大小不均时天花板不到 1。
    """
    import numpy as np
    from sklearn.neighbors import NearestNeighbors

    b = np.asarray(batch)
    levels, codes = np.unique(b, return_inverse=True)
    nn = NearestNeighbors(n_neighbors=min(k + 1, len(b))).fit(emb)
    idx = nn.kneighbors(return_distance=False)[:, :k]
    nb = codes[idx]
    counts = np.stack([(nb == j).sum(axis=1) for j in range(len(levels))], axis=1).astype(float)
    p = counts / counts.sum(axis=1, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        ent = -np.nansum(np.where(p > 0, p * np.log(p), 0.0), axis=1)
    # 上限 = 全局批次比例下的熵
    glob = np.array([(codes == j).mean() for j in range(len(levels))])
    ceil = -np.sum(glob[glob > 0] * np.log(glob[glob > 0]))
    return ent / ceil if ceil > 0 else ent


def knn_label_purity(emb, labels, k: int = 30):
    """kNN 邻域内与自己同细胞类型的比例——生物结构有没有被整合抹平。"""
    import numpy as np
    from sklearn.neighbors import NearestNeighbors

    y = np.asarray(labels)
    _, codes = np.unique(y, return_inverse=True)
    nn = NearestNeighbors(n_neighbors=min(k + 1, len(y))).fit(emb)
    idx = nn.kneighbors(return_distance=False)[:, :k]
    return (codes[idx] == codes[:, None]).mean(axis=1)


def geometry_rho(emb, ref, n_sub: int = 700, seed: int = SEED):
    """全局几何保真:成对距离在 latent 与参考空间的 Spearman 相关(子采样以控开销)。

    这是 CONCORD 论文的核心主张(global geometry),所以必须显式测,不能只看 UMAP 好看。
    """
    import numpy as np
    from scipy.spatial.distance import pdist
    from scipy.stats import spearmanr

    rng = np.random.default_rng(seed)
    n = emb.shape[0]
    sel = rng.choice(n, size=min(n_sub, n), replace=False)
    d1 = pdist(np.asarray(emb)[sel])
    d2 = pdist(np.asarray(ref)[sel])
    return float(spearmanr(d1, d2).correlation)


def evaluate(emb, adata, ref=None, k: int = 30):
    import numpy as np
    from sklearn.cluster import KMeans
    from sklearn.manifold import trustworthiness
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

    ct = np.asarray(adata.obs["cell_type"])
    bt = np.asarray(adata.obs["batch"])
    ent = knn_batch_entropy(emb, bt, k=k)
    pur = knn_label_purity(emb, ct, k=k)

    km = KMeans(n_clusters=len(np.unique(ct)), n_init=10, random_state=SEED).fit_predict(emb)
    m = {
        "batch_mixing_entropy": float(np.mean(ent)),
        "celltype_knn_purity": float(np.mean(pur)),
        "ARI_celltype": float(adjusted_rand_score(ct, km)),
        "NMI_celltype": float(normalized_mutual_info_score(ct, km)),
    }
    if ref is not None:
        sub = np.random.default_rng(SEED).choice(emb.shape[0], size=min(600, emb.shape[0]),
                                                 replace=False)
        m["trustworthiness"] = float(trustworthiness(np.asarray(ref)[sub], np.asarray(emb)[sub],
                                                     n_neighbors=min(15, len(sub) - 1)))
        m["global_geometry_rho"] = geometry_rho(emb, ref)
    return m, ent, pur


# ---------------------------------------------------------------------------
# 出图(框架样式;库规矩:不用条形图)
# ---------------------------------------------------------------------------
def _umap2d(emb, seed: int = SEED):
    import numpy as np
    import umap

    return umap.UMAP(n_neighbors=15, min_dist=0.3, random_state=seed).fit_transform(
        np.asarray(emb))


def fig_embeddings(embs, adata, outdir):
    """散点矩阵:行=方法,列=按批次上色 / 按细胞类型上色。整合好坏一眼能看。"""
    import matplotlib.pyplot as plt
    import numpy as np

    names = list(embs)
    fig, axes = plt.subplots(len(names), 2, figsize=(NATURE_W2 * 0.72, 2.5 * len(names)),
                             squeeze=False)
    for i, nm in enumerate(names):
        xy = _umap2d(embs[nm])
        for j, key in enumerate(["batch", "cell_type"]):
            ax = axes[i, j]
            lv = list(dict.fromkeys(np.asarray(adata.obs[key])))
            cols = pal(len(lv), "okabe_ito" if key == "batch" else "npg")
            for c, lev in zip(cols, lv):
                m = np.asarray(adata.obs[key]) == lev
                ax.scatter(xy[m, 0], xy[m, 1], s=3, c=c, lw=0, alpha=0.75,
                           label=str(lev), rasterized=True)
            ax.set_xticks([]); ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_visible(True)
            if i == 0:
                ax.set_title("Coloured by " + key.replace("_", " "), fontsize=9)
            if j == 0:
                ax.set_ylabel(nm, fontsize=9, fontweight="bold")
            if i == len(names) - 1:
                ax.legend(markerscale=3, fontsize=6, loc="upper center",
                          bbox_to_anchor=(0.5, -0.04), ncol=len(lv))
    fig.suptitle("UMAP of each latent space", fontsize=11, fontweight="bold")
    fig.tight_layout()
    save_fig(fig, os.path.join(outdir, "563_fig1_embeddings"))
    plt.close(fig)


def fig_metric_heatmap(mdf, outdir):
    """指标热图:方法 × 指标,列内 min-max 归一化(绝对值另见 CSV)。"""
    import matplotlib.pyplot as plt
    import numpy as np

    Z = mdf.copy()
    rng_ = Z.max() - Z.min()
    Zn = (Z - Z.min()) / rng_.replace(0, np.nan)
    Zn = Zn.fillna(0.5)

    fig, ax = plt.subplots(figsize=(NATURE_W1 * 1.45, 0.55 * len(Z) + 1.6))
    im = ax.imshow(Zn.values, cmap="viridis", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(Z.shape[1]))
    ax.set_xticklabels([c.replace("_", "\n") for c in Z.columns], fontsize=7)
    ax.set_yticks(range(Z.shape[0]))
    ax.set_yticklabels(Z.index, fontsize=8)
    for i in range(Z.shape[0]):
        for j in range(Z.shape[1]):
            ax.text(j, i, f"{Z.values[i, j]:.2f}", ha="center", va="center", fontsize=6.5,
                    color="white" if Zn.values[i, j] < 0.55 else "black")
    ax.set_title("Integration metrics (cell = raw value, colour = column-scaled)",
                 fontsize=9)
    fig.colorbar(im, ax=ax, shrink=0.7, label="column-scaled")
    fig.tight_layout()
    save_fig(fig, os.path.join(outdir, "563_fig2_metric_heatmap"))
    plt.close(fig)


def fig_tradeoff(mdf, outdir):
    """权衡散点:横轴批次混合、纵轴生物保真。右上角才是真赢,只往右是过度整合。"""
    import matplotlib.pyplot as plt

    cols = pal(len(mdf), "npg")
    fig, ax = plt.subplots(figsize=(NATURE_W1 * 1.25, NATURE_W1 * 1.05))
    for c, (nm, row) in zip(cols, mdf.iterrows()):
        ax.scatter(row["batch_mixing_entropy"], row["celltype_knn_purity"], s=110,
                   c=c, edgecolor="black", lw=0.8, zorder=3, label=nm)
        ax.annotate(nm, (row["batch_mixing_entropy"], row["celltype_knn_purity"]),
                    textcoords="offset points", xytext=(7, 6), fontsize=7)
    ax.set_xlabel("Batch mixing entropy (kNN, higher = better mixed)")
    ax.set_ylabel("Cell-type kNN purity (higher = bio preserved)")
    ax.set_title("Batch removal vs biological conservation", fontsize=9)
    fig.tight_layout()
    save_fig(fig, os.path.join(outdir, "563_fig3_tradeoff"))
    plt.close(fig)


def fig_percell_violin(ent_map, outdir):
    """每细胞批次混合熵的 violin + 抖动散点(raincloud 风格),看的是分布不是均值。"""
    import matplotlib.pyplot as plt
    import numpy as np

    names = list(ent_map)
    data = [ent_map[n] for n in names]
    cols = pal(len(names), "npg")
    rng = np.random.default_rng(SEED)

    fig, ax = plt.subplots(figsize=(NATURE_W1 * 1.5, NATURE_W1 * 1.0))
    vp = ax.violinplot(data, showextrema=False, widths=0.85)
    for body, c in zip(vp["bodies"], cols):
        body.set_facecolor(c); body.set_alpha(0.35); body.set_edgecolor(c); body.set_lw(1.0)
    for i, (d, c) in enumerate(zip(data, cols), start=1):
        sub = rng.choice(d, size=min(220, len(d)), replace=False)
        ax.scatter(i + rng.normal(0, 0.045, len(sub)), sub, s=4, c=c, alpha=0.5, lw=0,
                   rasterized=True)
        q1, med, q3 = np.percentile(d, [25, 50, 75])
        ax.plot([i, i], [q1, q3], color="black", lw=2.2, zorder=4, solid_capstyle="round")
        ax.scatter([i], [med], s=22, c="white", edgecolor="black", zorder=5, lw=0.9)
    ax.set_xticks(range(1, len(names) + 1))
    ax.set_xticklabels(names, rotation=18, ha="right", fontsize=7.5)
    ax.set_ylabel("Per-cell batch mixing entropy")
    ax.set_title("Distribution of local batch mixing", fontsize=9)
    fig.tight_layout()
    save_fig(fig, os.path.join(outdir, "563_fig4_percell_entropy"))
    plt.close(fig)


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--indir", default=EXAMPLE, help="目录,含 counts.csv + cell_meta.csv")
    ap.add_argument("--h5ad", default=None, help="改用 h5ad 输入(需含 batch/cell_type obs 列)")
    ap.add_argument("--outdir", default=RESULTS)
    ap.add_argument("--batch-key", default="batch")
    ap.add_argument("--label-key", default="cell_type")
    ap.add_argument("--n-comps", type=int, default=30)
    ap.add_argument("--knn", type=int, default=30)
    ap.add_argument("--run-concord", action="store_true",
                    help="尝试真实 CONCORD(需 pip install concord-sc,建议 GPU)")
    ap.add_argument("--concord-epochs", type=int, default=15)
    ap.add_argument("--regen-example", action="store_true", help="重新生成 example_data/")
    a = ap.parse_args()

    import numpy as np
    import pandas as pd

    np.random.seed(SEED)
    set_pub_style(base_size=9)
    os.makedirs(a.outdir, exist_ok=True)
    os.makedirs(ASSETS, exist_ok=True)

    # --- 输入 ---
    if a.h5ad:
        import anndata as ad
        adata = ad.read_h5ad(a.h5ad)
    else:
        if a.regen_example or not os.path.exists(os.path.join(a.indir, "counts.csv")):
            print("[563] generating synthetic example_data/ ...")
            write_example(a.indir, make_synthetic())
        adata = read_example(a.indir)
    for key in (a.batch_key, a.label_key):
        if key not in adata.obs:
            sys.exit(f"[563] input lacks obs['{key}']")
    adata.obs["batch"] = adata.obs[a.batch_key].astype("category")
    adata.obs["cell_type"] = adata.obs[a.label_key].astype("category")
    print(f"[563] {adata.n_obs} cells x {adata.n_vars} genes · "
          f"{adata.obs['batch'].nunique()} batches · {adata.obs['cell_type'].nunique()} labels")

    # --- 共用预处理 ---
    adata = preprocess(adata)

    # 几何参考空间:合成数据用无批次的干净信号;真实数据没有 ground truth,
    # 退回用未整合 PCA 作参考(此时 trustworthiness/geometry 只是"相对未整合的形变量")
    ref = np.asarray(adata.obsm["X_true"]) if "X_true" in adata.obsm else None

    # --- 基线 ---
    print("[563] Step 1 baselines: PCA (no integration) / batch-centred PCA / ComBat+PCA")
    embs = {
        "PCA (unintegrated)": embed_pca(adata, a.n_comps),
        "Batch-centred PCA": embed_batch_center(adata, "batch", a.n_comps),
    }
    try:
        embs["ComBat + PCA"] = embed_combat(adata, "batch", a.n_comps)
    except Exception as e:
        print(f"       ComBat skipped: {type(e).__name__}: {e}")
    if ref is None:
        # 真实数据没有 ground truth。此时**不报** trustworthiness / geometry_rho ——
        # 若拿未整合 PCA 当参照,未整合方法会自比得 1.000,是自证的假赢。
        print("       no ground-truth geometry in input -> trustworthiness / "
              "global_geometry_rho are omitted (no honest reference exists)")

    # --- CONCORD(守卫) ---
    concord_info = {"status": "not requested (--run-concord)"}
    if a.run_concord:
        print("[563] Step 2 CONCORD (guarded)")
        emb, concord_info = embed_concord(adata, "batch", latent_dim=a.n_comps,
                                          n_epochs=a.concord_epochs)
        if emb is not None:
            embs["CONCORD"] = emb
        print(f"       {concord_info}")
    else:
        print("[563] Step 2 CONCORD not requested; baselines only "
              "(pip install concord-sc, then --run-concord)")

    # --- 评估 ---
    print("[563] Step 3 evaluate: batch mixing / bio conservation / geometry")
    rows, ent_map = {}, {}
    for nm, E in embs.items():
        m, ent, _ = evaluate(E, adata, ref=ref, k=a.knn)
        rows[nm] = m
        ent_map[nm] = ent
        print(f"       {nm:22s} " + "  ".join(f"{k}={v:.3f}" for k, v in m.items()))
    mdf = pd.DataFrame(rows).T
    mdf.index.name = "method"
    mdf.to_csv(os.path.join(a.outdir, "563_integration_metrics.csv"))
    pd.DataFrame(ent_map, index=adata.obs_names).to_csv(
        os.path.join(a.outdir, "563_percell_batch_entropy.csv"))
    for nm, E in embs.items():
        np.save(os.path.join(a.outdir, "embedding_" + nm.split()[0].lower() + ".npy"), E)

    # --- 出图 ---
    print("[563] Step 4 figures -> assets/")
    fig_embeddings(embs, adata, ASSETS)
    fig_metric_heatmap(mdf, ASSETS)
    fig_tradeoff(mdf, ASSETS)
    fig_percell_violin(ent_map, ASSETS)

    # 依赖版本快照(铁律6:可复现)
    import importlib.metadata as im
    versions = {}
    for p in ("scanpy", "anndata", "scikit-learn", "umap-learn", "numpy", "pandas",
              "scipy", "matplotlib", "concord-sc"):
        try:
            versions[p] = im.version(p)
        except Exception:
            versions[p] = "not installed"
    versions["python"] = sys.version.split()[0]

    with open(os.path.join(a.outdir, "563_summary.json"), "w") as fh:
        json.dump({"n_cells": int(adata.n_obs), "n_genes": int(adata.n_vars),
                   "seed": SEED, "versions": versions,
                   "methods": list(embs), "metrics": mdf.to_dict(orient="index"),
                   "concord": concord_info,
                   # ref is None 时几何指标是被**省略**的,不存在任何替代参照系,
                   # 早先写 "unintegrated PCA" 是错的(代码从未这样回退)。
                   "geometry_reference": "X_true (synthetic)" if ref is not None
                                         else "none (geometry metrics omitted)"},
                  fh, indent=1, default=str)
    print(f"[563] done -> {a.outdir}")


if __name__ == "__main__":
    main()
