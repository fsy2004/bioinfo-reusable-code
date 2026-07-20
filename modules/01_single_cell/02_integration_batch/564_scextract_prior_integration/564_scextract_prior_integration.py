"""564 · scExtract prior-informed integration (scanorama_prior / cellhint_prior).

上游: Wu Y & Tang F, Genome Biology 2025, doi:10.1186/s13059-025-03639-x (PMID 40537825)
仓库: https://github.com/yxwucq/scExtract

本模块做两件事:
  1. **可跑基线 (always, CPU, 本机依赖即可)** —— 未校正 PCA vs ComBat 的批次整合对照,
     用「批次混合度 (batch entropy)」×「生物学保真度 (cell-type kNN purity)」双轴评估,
     并单独追踪 **batch-specific 稀有细胞类型** 的保真度。这是整合类方法的最低对照:
     任何"更好"的说法都必须先赢过它,尤其不能靠抹平真实差异来换混合度。
  2. **scExtract 守卫式封装** —— 未安装/缺 embedding 字典时优雅退出并打印真实安装命令,
     不臆造 API。函数签名来自实际读取的上游源码(URL 见 README)。

★ 只封装 scanorama_prior / cellhint_prior 两个先验整合算法。
  **不构建 LLM 自动注释链路**(需外部 API key,且细胞类型判定不应交给 LLM)。
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
FRAMEWORK = os.path.abspath(os.path.join(HERE, "..", "..", "..", "_framework"))
sys.path.insert(0, FRAMEWORK)

from pubstyle import (  # noqa: E402
    CMAP_CONT, NATURE_W1, NATURE_W2, pal, save_fig, set_pub_style,
)

SEED = 0


# ---------------------------------------------------------------------------
# 合成示例数据: 3 个批次 × 4 个共享细胞类型 + 1 个只存在于 batch3 的稀有类型。
# 稀有类型是关键设计: 它专门用来检验"整合有没有把真实差异一起抹平"。
# ---------------------------------------------------------------------------
def make_synthetic(n_genes: int = 300, seed: int = SEED):
    import anndata as ad
    import numpy as np
    import pandas as pd

    rng = np.random.default_rng(seed)
    shared = ["Tcell", "Bcell", "Myeloid", "Fibroblast"]
    # (batch, celltype, n_cells)
    spec = []
    for b, nb in zip(["batch1", "batch2", "batch3"], [140, 160, 150]):
        for ct in shared:
            spec.append((b, ct, nb))
    spec.append(("batch3", "RareStromal", 90))   # batch-specific,不可被抹掉

    # 每个细胞类型有一组标志基因(程序性上调),批次有整体乘性+加性漂移。
    # 关键: 生物学信号幅度 (~1.1) 明显小于批次漂移 —— 这才是真实 scRNA 整合的处境;
    # 信号若压过噪声,任何方法都满分,基线就失去了鉴别力。
    ct_names = shared + ["RareStromal"]
    marker_idx = {ct: rng.choice(n_genes, 30, replace=False) for ct in shared}
    # RareStromal 与 Fibroblast 共享 2/3 标志基因 —— 过度校正最容易把它并进 Fibroblast,
    # 这正是 "aligning without flattening real differences" 要考的点。
    fib = marker_idx["Fibroblast"]
    marker_idx["RareStromal"] = np.concatenate(
        [fib[:20], rng.choice(np.setdiff1d(np.arange(n_genes), fib), 10, replace=False)])

    batch_scale = {"batch1": 1.0, "batch2": 1.9, "batch3": 0.55}
    batch_shift = {"batch1": 0.0, "batch2": 1.2, "batch3": -0.4}
    batch_gene = {b: rng.normal(0, 1.1, n_genes) for b in batch_scale}

    X, obs_b, obs_c = [], [], []
    for b, ct, n in spec:
        base = rng.normal(2.0, 1.0, (n, n_genes))                  # 细胞级噪声
        base[:, marker_idx[ct]] += rng.normal(1.1, 0.45, (n, 30))  # 生物学信号(弱)
        base = base * batch_scale[b] + batch_shift[b] + batch_gene[b][None, :]
        X.append(base)
        obs_b += [b] * n
        obs_c += [ct] * n
    X = np.vstack(X).astype("float32")
    X -= X.min()                                   # 保持非负,像 log-normalized 表达

    A = ad.AnnData(X)
    A.var_names = [f"G{i:04d}" for i in range(n_genes)]
    A.obs_names = [f"C{i:05d}" for i in range(A.n_obs)]
    A.obs["batch"] = pd.Categorical(obs_b)
    A.obs["cell_type"] = pd.Categorical(obs_c)
    A.uns["synthetic"] = "synthetic, for demo only"
    return A


# ---------------------------------------------------------------------------
# 评估指标
# ---------------------------------------------------------------------------
def _knn_scores(emb, batch, celltype, k: int = 30):
    """返回 (批次混合熵, 细胞类型 kNN 保真度, 每类型保真度 dict)。

    批次混合熵: 每个细胞 k 近邻中批次标签的 Shannon 熵 / log(n_batch),1=完全混合。
    类型保真度: k 近邻中与自身同类型的比例,1=生物学结构完好。
    两者是权衡关系 —— 单看混合度会奖励"把所有细胞搅成一团"的坏整合。
    """
    import numpy as np
    from sklearn.neighbors import NearestNeighbors

    nn = NearestNeighbors(n_neighbors=k + 1).fit(emb)
    idx = nn.kneighbors(emb, return_distance=False)[:, 1:]

    b = np.asarray(batch)
    c = np.asarray(celltype)
    ub = np.unique(b)

    ent = []
    for row in b[idx]:
        p = np.array([(row == u).mean() for u in ub])
        p = p[p > 0]
        ent.append(-(p * np.log(p)).sum())
    ent = np.array(ent) / np.log(len(ub))

    pur = (c[idx] == c[:, None]).mean(axis=1)
    per_ct = {ct: float(pur[c == ct].mean()) for ct in np.unique(c)}
    return float(ent.mean()), float(pur.mean()), per_ct


def _ari(emb, celltype, seed: int = SEED):
    """KMeans(k=真类别数) 对真标签的 ARI。用 KMeans 而非 leiden,避免 leidenalg 依赖。"""
    import numpy as np
    from sklearn.cluster import KMeans
    from sklearn.metrics import adjusted_rand_score

    k = len(np.unique(celltype))
    lab = KMeans(n_clusters=k, n_init=10, random_state=seed).fit_predict(emb)
    return float(adjusted_rand_score(celltype, lab))


# ---------------------------------------------------------------------------
# 基线: 未校正 PCA vs ComBat (+ Harmony 若本机恰好有 harmonypy)
# ---------------------------------------------------------------------------
def run_baseline(adata, n_pcs: int = 30, k: int = 30):
    import numpy as np
    import scanpy as sc

    embeds = {}

    un = adata.copy()
    sc.pp.scale(un, max_value=10)
    sc.tl.pca(un, n_comps=n_pcs, svd_solver="arpack", random_state=SEED)
    embeds["Uncorrected PCA"] = np.asarray(un.obsm["X_pca"])

    cb = adata.copy()
    try:
        sc.pp.combat(cb, key="batch")
        sc.pp.scale(cb, max_value=10)
        sc.tl.pca(cb, n_comps=n_pcs, svd_solver="arpack", random_state=SEED)
        embeds["ComBat"] = np.asarray(cb.obsm["X_pca"])
    except Exception as e:                     # 基线不因单个方法失败而整体失败
        print(f"       ComBat skipped: {type(e).__name__}: {e}")

    # Harmony 是 scExtract 的依赖之一;本机若没有就跳过,绝不安装
    try:
        import harmonypy  # noqa: F401
        hm = adata.copy()
        hm.obsm["X_pca"] = embeds["Uncorrected PCA"]
        sc.external.pp.harmony_integrate(hm, key="batch", random_state=SEED)
        embeds["Harmony"] = np.asarray(hm.obsm["X_pca_harmony"])
    except Exception as e:
        print(f"       Harmony skipped ({type(e).__name__}): harmonypy not installed locally")

    rows = []
    for name, emb in embeds.items():
        ent, pur, per_ct = _knn_scores(emb, adata.obs["batch"], adata.obs["cell_type"], k=k)
        rows.append({
            "method": name,
            "batch_mixing_entropy": round(ent, 4),
            "celltype_knn_purity": round(pur, 4),
            "rare_type_purity": round(per_ct.get("RareStromal", float("nan")), 4),
            "kmeans_ARI_vs_truth": round(_ari(emb, adata.obs["cell_type"]), 4),
        })
    return embeds, rows


# ---------------------------------------------------------------------------
# scExtract 守卫式封装 —— 只走 scanorama_prior / cellhint_prior
# ---------------------------------------------------------------------------
def _check_upstream_preconditions(file_list, method: str):
    """上游对输入 h5ad 的硬性要求(实读自 integration/integrate.py,非推测)。

    这些要求上游没有做友好检查,不满足会在包内部抛 AttributeError / KeyError,
    所以在调用前先自查并给出可读的原因。逐条源码依据(本地克隆 commit f3ef7fb / v0.2.0,
    行号逐条 grep 复核过):
      · integrate.py:258 / :310  `assert len(file_list) == 1` —— scanorama(_prior) 只吃合并后的单文件
      · integrate.py:263 / :315 / :353  `adata_all.raw.to_adata()` —— 单文件路径必须带 .raw
      · integrate.py:58  `csr_matrix(adata.raw.X)` —— 多文件 merge_datasets 同样必须带 .raw
      · integrate.py:265 / :317 / :355  `highly_variable_genes(batch_key='Dataset')` —— 单文件需 obs['Dataset']
      · integrate.py:42-43  多文件路径 obs['Dataset'] 缺失时由 merge_datasets 用文件名自动补,故不强制
      · integrate.py:359  `adata_all.obs['cell_type']` —— 单文件路径需 obs['cell_type']
      · integrate.py:45-52 (merge_datasets) —— 多文件路径需 cell_type / leiden / louvain 之一
    """
    import anndata as ad

    probs = []
    single = len(file_list) == 1
    if method in ("scanorama_prior", "scanorama") and not single:
        probs.append("scanorama(_prior) asserts len(file_list) == 1 (integrate.py:258/:310); "
                     "merge with --method cellhint_prior first")
    for f in file_list:
        if not os.path.exists(f):
            probs.append(f"{f}: not found")
            continue
        A = ad.read_h5ad(f, backed="r")
        if A.raw is None:
            probs.append(f"{os.path.basename(f)}: adata.raw is None; upstream calls "
                         "adata_all.raw.to_adata() (integrate.py:263) / "
                         "csr_matrix(adata.raw.X) (integrate.py:58)")
        if single and "Dataset" not in A.obs:
            probs.append(f"{os.path.basename(f)}: obs['Dataset'] missing; upstream uses it as "
                         "the batch key in highly_variable_genes (integrate.py:265). "
                         "(Only the single-file path needs it — merge_datasets auto-fills it "
                         "from the file name, integrate.py:42-43)")
        if single and "cell_type" not in A.obs:
            probs.append(f"{os.path.basename(f)}: obs['cell_type'] missing (integrate.py:359)")
        if not single and not ({"cell_type", "leiden", "louvain"} & set(A.obs.columns)):
            probs.append(f"{os.path.basename(f)}: needs one of cell_type/leiden/louvain "
                         "(merge_datasets, integrate.py:45-52)")
    return probs


def run_scextract(file_list, method: str, output_path: str,
                  embedding_dict_path=None, config_path="config.ini", **kw):
    """调用上游真实函数 scextract.integration.integrate.integrate_processed_datasets。

    签名实读自源码 src/scextract/integration/integrate.py:160-174 (commit f3ef7fb, v0.2.0):
        integrate_processed_datasets(file_list: List[str], method: str, output_path: str,
                                     config_path: str = 'config.ini',
                                     alignment_path: Optional[str] = None,
                                     embedding_dict_path: Optional[str] = None,
                                     downsample: Optional[bool] = False,
                                     downsample_cells_per_label: Optional[int] = 1000,
                                     search_factor: int = 5, approx: bool = False,
                                     use_gpu: bool = False, batch_size: int = 5000,
                                     dimred: int = 100, use_pct: bool = False, **kwargs) -> None
    method 的合法取值(实读自 utils/parse_args.py:101):
        ['scExtract', 'scanorama_prior', 'scanorama', 'cellhint_prior', 'cellhint']
    本模块只允许 prior 两支;'scExtract' 走 LLM 注释链路,刻意不暴露。
    """
    allowed = {"scanorama_prior", "cellhint_prior", "scanorama", "cellhint"}
    if method not in allowed:
        return {"status": "refused",
                "reason": f"method must be one of {sorted(allowed)}; "
                          "'scExtract' (LLM annotation path) is intentionally not wrapped here"}
    try:
        from scextract.integration.integrate import integrate_processed_datasets
    except ImportError as e:
        return {"status": "skipped",
                "reason": f"scExtract not installed ({e.name}). "
                          "git clone https://github.com/yxwucq/scExtract && pip install -e . "
                          "(upstream README Step1). The integration backends are NOT in its "
                          "pyproject dependencies; integrate.py imports the module names "
                          "`scanorama_prior` / `cellhint_prior` (separate forks at "
                          "github.com/yxwucq/scanorama_prior and github.com/yxwucq/cellhint_prior), "
                          "or plain `scanorama` / `cellhint` for the no-prior methods; "
                          "install whichever matches --method separately"}
    if method.endswith("_prior") and not embedding_dict_path:
        return {"status": "skipped",
                "reason": "*_prior needs --embedding-dict (cell-type embedding pickle). "
                          "Upstream generates it with `scExtract extract_celltype_embedding "
                          "--file_list <h5ad> --cell_type_column cell_type "
                          "--output_embedding_pkl <pkl>`, which calls an LLM API. "
                          "Supply your own pickle, or run that step yourself — this module "
                          "does not build the LLM path."}
    probs = _check_upstream_preconditions(file_list, method)
    if probs:
        return {"status": "skipped",
                "reason": "input does not meet upstream preconditions: " + "; ".join(probs)}
    integrate_processed_datasets(file_list=list(file_list), method=method,
                                 output_path=output_path, config_path=config_path,
                                 embedding_dict_path=embedding_dict_path, **kw)
    return {"status": "done", "method": method, "output_path": output_path}


# ---------------------------------------------------------------------------
# 出图 (无条形图: 散点 / dumbbell / heatmap)
# ---------------------------------------------------------------------------
def plot_all(adata, embeds, rows, outdir):
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import umap

    set_pub_style(base_size=9)
    df = pd.DataFrame(rows)
    cts = list(pd.unique(adata.obs["cell_type"]))
    bts = list(pd.unique(adata.obs["batch"]))
    c_ct = dict(zip(cts, pal(len(cts), "npg")))
    c_bt = dict(zip(bts, pal(len(bts), "okabe_ito")))

    # --- Fig 1: UMAP 网格 (行=方法, 列=按批次/按类型着色) ---
    names = list(embeds)
    fig, axes = plt.subplots(len(names), 2, figsize=(NATURE_W2 * 0.72, 3.0 * len(names)),
                             squeeze=False)
    for i, nm in enumerate(names):
        U = umap.UMAP(random_state=SEED, n_neighbors=25, min_dist=0.3).fit_transform(embeds[nm])
        for j, (key, cmapd) in enumerate([("batch", c_bt), ("cell_type", c_ct)]):
            ax = axes[i, j]
            for lv in (bts if key == "batch" else cts):
                m = (adata.obs[key] == lv).values
                ax.scatter(U[m, 0], U[m, 1], s=3, lw=0, c=cmapd[lv], label=lv, alpha=0.75)
            ax.set_xticks([]); ax.set_yticks([])
            ax.set_title(f"{nm} — by {key.replace('_', ' ')}", fontsize=9)
            if i == 0:
                ax.legend(markerscale=3, fontsize=6, loc="upper right", handletextpad=0.2)
    fig.tight_layout()
    save_fig(fig, os.path.join(outdir, "fig1_umap_batch_vs_celltype"))
    plt.close(fig)

    # --- Fig 2: 权衡散点 (x=批次混合, y=类型保真; 点大小=稀有类型保真) ---
    fig, ax = plt.subplots(figsize=(NATURE_W1 * 1.25, NATURE_W1 * 1.05))
    cols = pal(len(df), "lancet")
    for (_, r), c in zip(df.iterrows(), cols):
        ax.scatter(r.batch_mixing_entropy, r.celltype_knn_purity,
                   s=40 + 260 * r.rare_type_purity, c=c, lw=1.1, edgecolor="black", zorder=3)
        ax.annotate(r.method, (r.batch_mixing_entropy, r.celltype_knn_purity),
                    textcoords="offset points", xytext=(9, 9), fontsize=8)
    ax.set_xlabel("Batch mixing entropy  (higher = better mixed)")
    ax.set_ylabel("Cell-type kNN purity  (higher = biology kept)")
    ax.set_title("Integration trade-off\npoint size = rare-type purity", fontsize=9)
    ax.margins(0.22)
    fig.tight_layout()
    save_fig(fig, os.path.join(outdir, "fig2_mixing_vs_purity_tradeoff"))
    plt.close(fig)

    # --- Fig 3: dumbbell —— 每个方法从 uncorrected 出发的移动量 ---
    metrics = ["batch_mixing_entropy", "celltype_knn_purity", "rare_type_purity"]
    labels = ["Batch mixing", "Cell-type purity", "Rare-type purity"]
    # 高度按"实际会画的行数"算(方法数 -1 个参照 × 指标数),否则单方法时图会被拉得极空
    n_rows = max(1, len(df) - 1) * len(metrics)
    fig, ax = plt.subplots(figsize=(NATURE_W1 * 1.5, max(2.4, 0.5 * n_rows + 1.3)))
    ref = df[df.method == "Uncorrected PCA"].iloc[0]
    ypos, yticks, ylabs = 0, [], []
    mcol = dict(zip(df.method, pal(len(df), "lancet")))
    for met, lab in zip(metrics, labels):
        for _, r in df.iterrows():
            if r.method == "Uncorrected PCA":
                continue
            ax.plot([ref[met], r[met]], [ypos, ypos], color="#BBBBBB", lw=1.6, zorder=1)
            ax.scatter(ref[met], ypos, s=42, c="#999999", lw=0.8, edgecolor="black", zorder=3)
            ax.scatter(r[met], ypos, s=52, c=mcol[r.method], lw=0.8, edgecolor="black", zorder=3)
            yticks.append(ypos); ylabs.append(f"{lab}\n{r.method}")
            ypos += 1
        ypos += 0.6
    ax.set_yticks(yticks); ax.set_yticklabels(ylabs, fontsize=7)
    ax.set_xlabel("Score (grey = uncorrected PCA reference)")
    ax.set_title("Shift from the uncorrected baseline", fontsize=9)
    ax.margins(x=0.14)
    fig.tight_layout()
    save_fig(fig, os.path.join(outdir, "fig3_shift_from_baseline"))
    plt.close(fig)

    # --- Fig 4: 跨批次细胞类型质心相似度 heatmap (先验整合权重的经验类比物) ---
    best = df.sort_values("celltype_knn_purity", ascending=False).iloc[0].method
    emb = embeds[best]
    keys, cents = [], []
    for b in bts:
        for ct in cts:
            m = ((adata.obs["batch"] == b) & (adata.obs["cell_type"] == ct)).values
            if m.sum() >= 5:
                keys.append(f"{b}|{ct}")
                cents.append(emb[m].mean(0))
    C = np.vstack(cents)
    Cn = C / np.linalg.norm(C, axis=1, keepdims=True)
    S = Cn @ Cn.T
    fig, ax = plt.subplots(figsize=(NATURE_W1 * 1.7, NATURE_W1 * 1.55))
    im = ax.imshow(S, cmap=CMAP_CONT, vmin=np.percentile(S, 2), vmax=1)
    ax.set_xticks(range(len(keys))); ax.set_xticklabels(keys, rotation=90, fontsize=5.5)
    ax.set_yticks(range(len(keys))); ax.set_yticklabels(keys, fontsize=5.5)
    ax.set_title(f"Cross-batch cell-type centroid similarity ({best})", fontsize=9)
    fig.colorbar(im, ax=ax, shrink=0.7, label="cosine similarity")
    fig.tight_layout()
    save_fig(fig, os.path.join(outdir, "fig4_crossbatch_celltype_similarity"))
    plt.close(fig)

    return ["fig1_umap_batch_vs_celltype", "fig2_mixing_vs_purity_tradeoff",
            "fig3_shift_from_baseline", "fig4_crossbatch_celltype_similarity"]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--h5ad", default=os.path.join(EXAMPLE, "synthetic_3batch.h5ad"),
                    help="AnnData;需含 obs['batch'] 与 obs['cell_type'];缺失则自动生成合成示例")
    ap.add_argument("--batch-key", default="batch")
    ap.add_argument("--label-key", default="cell_type")
    ap.add_argument("--n-pcs", type=int, default=30)
    ap.add_argument("--k", type=int, default=30, help="kNN 指标的邻居数")
    ap.add_argument("--run-scextract", action="store_true",
                    help="尝试真实 scExtract 先验整合(需安装 + embedding 字典)")
    ap.add_argument("--method", default="scanorama_prior",
                    choices=["scanorama_prior", "cellhint_prior", "scanorama", "cellhint"])
    ap.add_argument("--embedding-dict", default=None, help="细胞类型 embedding pickle 路径")
    ap.add_argument("--outdir", default=RESULTS)
    a = ap.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    os.makedirs(ASSETS, exist_ok=True)
    os.makedirs(EXAMPLE, exist_ok=True)

    import anndata as ad

    if not os.path.exists(a.h5ad):
        print(f"[564] {os.path.basename(a.h5ad)} not found -> generating synthetic demo (seed={SEED})")
        adata = make_synthetic()
        adata.write_h5ad(a.h5ad, compression="gzip")   # 压缩后 ~1MB,适合入库
    else:
        adata = ad.read_h5ad(a.h5ad)
    for kcol in (a.batch_key, a.label_key):
        if kcol not in adata.obs:
            sys.exit(f"input lacks obs['{kcol}'] — prior-informed integration needs both "
                     "a batch key and a cell-type label key")
    adata.obs["batch"] = adata.obs[a.batch_key]
    adata.obs["cell_type"] = adata.obs[a.label_key]
    print(f"[564] {adata.n_obs} cells x {adata.n_vars} genes; "
          f"batches={list(adata.obs['batch'].unique())}")

    print("[564] baseline: uncorrected PCA vs ComBat (+Harmony if available)")
    embeds, rows = run_baseline(adata, n_pcs=a.n_pcs, k=a.k)
    import pandas as pd
    df = pd.DataFrame(rows)
    print(df.to_string(index=False))
    df.to_csv(os.path.join(a.outdir, "564_baseline_metrics.csv"), index=False)

    figs = plot_all(adata, embeds, rows, a.outdir)
    import shutil
    for f in figs:
        shutil.copy(os.path.join(a.outdir, f + ".png"), os.path.join(ASSETS, f + ".png"))

    sx = {"status": "not requested (--run-scextract)"}
    if a.run_scextract:
        print(f"[564] scExtract path: method={a.method}")
        sx = run_scextract([a.h5ad], a.method,
                           os.path.join(a.outdir, "scextract_integrated.h5ad"),
                           embedding_dict_path=a.embedding_dict)
        for k_, v in sx.items():
            print(f"       {k_}: {v}")

    with open(os.path.join(a.outdir, "564_summary.json"), "w") as fh:
        json.dump({"baseline_metrics": rows, "scextract": sx, "figures": figs},
                  fh, indent=1, default=str)
    print(f"[564] wrote {a.outdir} + {len(figs)} figures -> assets/")


if __name__ == "__main__":
    main()
