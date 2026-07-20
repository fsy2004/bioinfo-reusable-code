"""575 · SCALE —— 空间组学「多尺度」空间域识别(带朴素基线)。

SCALE (Yousefi et al., Nucleic Acids Research 2026, doi:10.1093/nar/gkaf1456, PMID 41495880)
用图神经表征学习 + 基于熵的搜索,在多个尺度上找「稳定」的空间域。
仓库 https://github.com/imsb-uke/scale

本模块两条路:
  A. 朴素基线(默认,本机 scanpy/sklearn/leidenalg 即可跑完):
     空间近邻平滑 → PCA → (邻域尺度 k × Leiden 分辨率) 网格 → 多随机种子 ARI 稳定性
     → 选稳定且簇数分层的两个尺度。这是「多尺度稳定域」这一想法的无 GNN 下限。
  B. SCALE 真身(--run-scale):守卫式封装。scale 包未安装/无 torch-geometric 时优雅退出
     并打印真实安装命令,绝不静默降级、绝不假装跑过。

API 来源(实际读过的文件,2026-07-20):
  https://raw.githubusercontent.com/imsb-uke/scale/HEAD/scale/__init__.py
  https://raw.githubusercontent.com/imsb-uke/scale/HEAD/scale/scale.py
  https://raw.githubusercontent.com/imsb-uke/scale/HEAD/scale/config.py
  https://raw.githubusercontent.com/imsb-uke/scale/HEAD/scale/training.py
  https://raw.githubusercontent.com/imsb-uke/scale/HEAD/scale/search/_stability.py
  https://raw.githubusercontent.com/imsb-uke/scale/HEAD/notebooks/vignette.ipynb
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
RESULTS = os.path.join(HERE, "results")
ASSETS = os.path.join(HERE, "assets")
EXAMPLE = os.path.join(HERE, "example_data", "spatial_counts_synthetic.csv")

sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..", "..", "_framework")))
from pubstyle import set_pub_style, save_fig, pal, CMAP_CONT  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402

SEED = 200  # 与 SCALE 默认 cfg.seed 一致


# ---------------------------------------------------------------- 输入
def load_table(path: str):
    """读 csv:行=spot/cell,列 = x, y, [domain_coarse, domain_fine], 其余为基因。"""
    df = pd.read_csv(path, index_col=0)
    for c in ("x", "y"):
        if c not in df.columns:
            sys.exit(f"输入缺少坐标列 '{c}'")
    meta_cols = [c for c in ("x", "y", "domain_coarse", "domain_fine") if c in df.columns]
    genes = [c for c in df.columns if c not in meta_cols]
    coords = df[["x", "y"]].to_numpy(float)
    X = df[genes].to_numpy(float)
    truth = {c: df[c].to_numpy() for c in ("domain_coarse", "domain_fine") if c in df.columns}
    return df.index.to_numpy(), coords, X, np.array(genes), truth


# ---------------------------------------------------------------- 基线
def spatial_smooth(X: np.ndarray, coords: np.ndarray, k: int) -> np.ndarray:
    """把每个 spot 的表达替换为其 k 个空间近邻(含自身)的均值。

    这是「让表征带上空间上下文」的最朴素做法 —— SCALE 用 GNN 在空间图上学表征,
    这里用一次固定权重的邻域平均代替。k 就是这条基线的「尺度」旋钮。
    """
    from sklearn.neighbors import NearestNeighbors

    nn = NearestNeighbors(n_neighbors=min(k, len(coords))).fit(coords)
    _, idx = nn.kneighbors(coords)
    return X[idx].mean(axis=1)


def leiden_labels(emb: np.ndarray, res: float, seed: int, n_neighbors: int = 15) -> np.ndarray:
    """在给定表征上跑 Leiden(scanpy 的 igraph flavor,与 SCALE 的聚类后端一致)。"""
    import scanpy as sc
    import anndata as ad

    A = ad.AnnData(emb.astype("float32"))
    sc.pp.neighbors(A, use_rep="X", n_neighbors=n_neighbors, random_state=seed)
    sc.tl.leiden(A, resolution=res, key_added="L", random_state=seed,
                 flavor="igraph", n_iterations=2, directed=False)
    return A.obs["L"].to_numpy()


def stability_grid(X: np.ndarray, coords: np.ndarray, k_set, res_set,
                   n_repeats: int = 3, n_pcs: int = 15):
    """(k × resolution) 网格,每格跑 n_repeats 个随机种子,稳定性 = 种子间平均成对 ARI。

    与 SCALE 的 calc_stability 同一思路(它也是用 adjusted_rand_score 衡量重复间一致性),
    但这里换成的是朴素平滑表征而非 GNN 嵌入,且不做熵搜索。
    """
    from itertools import combinations
    from sklearn.decomposition import PCA
    from sklearn.metrics import adjusted_rand_score

    stab = pd.DataFrame(index=pd.Index(k_set, name="k_neighbors"),
                        columns=[round(r, 3) for r in res_set], dtype=float)
    nclu = stab.copy()
    labels_store: dict[tuple, np.ndarray] = {}

    for k in k_set:
        emb = spatial_smooth(X, coords, k)
        emb = PCA(n_components=min(n_pcs, emb.shape[1] - 1), random_state=SEED).fit_transform(emb)
        for res in res_set:
            reps = [leiden_labels(emb, res, seed=s) for s in range(n_repeats)]
            aris = [adjusted_rand_score(a, b) for a, b in combinations(reps, 2)]
            stab.loc[k, round(res, 3)] = float(np.mean(aris)) if aris else 1.0
            nclu.loc[k, round(res, 3)] = float(np.mean([len(set(r)) for r in reps]))
            labels_store[(k, round(res, 3))] = reps[0]
    return stab, nclu, labels_store


def pick_levels(stab: pd.DataFrame, nclu: pd.DataFrame, top_frac: float = 0.15,
                min_nclusters_start: int = 3):
    """从稳定性高的候选里挑出两个「尺度层」:粗层(簇少)与细层(簇多)。

    诚实说明:SCALE 用的是基于熵的层级搜索(scale.search.calc_entropy),会检查跨层的
    嵌套/单调性。这里只做「取 top 稳定候选 → 按簇数取两端」的朴素代理,不是熵搜索。
    """
    flat = stab.stack().sort_values(ascending=False)
    n_top = max(2, int(round(top_frac * flat.size)))
    cand = flat.head(n_top).index.tolist()
    cand = [c for c in cand if nclu.loc[c[0], c[1]] >= min_nclusters_start] or cand
    cand_sorted = sorted(cand, key=lambda c: nclu.loc[c[0], c[1]])
    return {"level_0_coarse": cand_sorted[0], "level_1_fine": cand_sorted[-1]}


def nonspatial_control(X: np.ndarray, res_set, n_pcs: int = 15):
    """无空间信息对照:直接对原始表达做 PCA + Leiden。用来看空间平滑到底加了多少。"""
    from sklearn.decomposition import PCA

    emb = PCA(n_components=min(n_pcs, X.shape[1] - 1), random_state=SEED).fit_transform(X)
    return {round(r, 3): leiden_labels(emb, r, seed=0) for r in res_set}


# ---------------------------------------------------------------- SCALE 守卫路径
def run_scale_guarded(csv_path: str):
    """SCALE 真身。签名取自上游源码,但仍以官方 vignette 为准。

    上游真实入口(读自 scale/__init__.py 与 scale/scale.py):
        from scale import run_scale
        from scale.config import load_config, Config
        cfg = load_config()                 # 可改 cfg.distance_set / cfg.resolution_set / cfg.n_repeats
        run_scale(adata, cfg, use_svgs=True, use_hvgs=False,
                  sample_key=None, integration_method=None, layer=None, celltype_key=None,
                  spatial_key="spatial", n_levels=2, top_n=0.15, ...)
    结果落在 adata.obs['scale_l*'](各尺度域标签)、adata.obsm['scale_clusterings']、adata.uns['scale']。
    也可按 vignette 手动分步:train → select_best_lambdas → calc_clusterings → calc_stability → calc_entropy。
    注意:上游 vignette 里 calc_clusterings(adata, flavor=..., n_iterations=...) 少传了必需的 cfg
    形参,与当前 scale/clustering.py 签名不一致 —— 生产运行前请以官方最新教程为准,此处不固定。
    """
    try:
        import scale as scale_pkg  # noqa: F401
        from scale import run_scale  # noqa: F401
        from scale.config import load_config  # noqa: F401
    except ImportError as e:
        return {"status": "skipped",
                "reason": f"未安装 scale ({e.name})",
                "install": "git clone https://github.com/imsb-uke/scale.git && cd scale && poetry install"}
    try:
        import torch_geometric  # noqa: F401
    except ImportError:
        return {"status": "skipped",
                "reason": "缺 torch-geometric(SCALE 的 GNN 依赖)",
                "install": "见 https://pytorch-geometric.readthedocs.io 按本机 torch/CUDA 版本安装"}
    return {"status": "ready",
            "note": "scale 与 torch-geometric 均可导入;请按官方 notebooks/vignette.ipynb 组织 AnnData"
                    "(需 adata.obsm['spatial'])后调用 run_scale(adata, load_config())。",
            "input_csv": csv_path}


# ---------------------------------------------------------------- 出图
def fig_stability(stab: pd.DataFrame, picks: dict, out: str):
    """稳定性热图:x=Leiden resolution, y=空间邻域 k,色=种子间平均 ARI。"""
    M = stab.to_numpy(float)
    fig, ax = plt.subplots(figsize=(7.0, 3.4))
    # 色标下界跟着数据走:全部落在 0.8-1.0 时用 0-1 会把差异压平
    im = ax.imshow(M, aspect="auto", cmap=CMAP_CONT,
                   vmin=float(np.nanmin(M)), vmax=1.0, origin="lower")
    ax.set_xticks(range(stab.shape[1]))
    ax.set_xticklabels([f"{c:.2f}" for c in stab.columns], rotation=90, fontsize=7)
    ax.set_yticks(range(stab.shape[0]))
    ax.set_yticklabels(stab.index, fontsize=8)
    ax.set_xlabel("Leiden resolution")
    ax.set_ylabel("Spatial neighbourhood k")
    ax.set_title("Stability across scales (mean pairwise ARI)", fontsize=10)
    for name, (k, r) in picks.items():
        xi, yi = list(stab.columns).index(r), list(stab.index).index(k)
        ax.scatter(xi, yi, s=90, facecolors="none", edgecolors="#E64B35", linewidths=1.8)
        # 靠右的标注往左放,免得被画布裁掉
        right = xi > stab.shape[1] * 0.6
        ax.annotate(name.replace("level_", "L").replace("_", " "), (xi, yi),
                    textcoords="offset points", xytext=(-8 if right else 8, 7),
                    ha="right" if right else "left", fontsize=7, color="#E64B35")
    fig.colorbar(im, ax=ax, label="stability (ARI)", fraction=0.03, pad=0.02)
    fig.tight_layout()
    save_fig(fig, out)


def fig_domains(coords, truth, labels_by_level, out: str):
    """空间散点:真值(粗/细)与基线在两个尺度上给出的域。"""
    panels = []
    if "domain_coarse" in truth:
        panels.append(("Ground truth · coarse", truth["domain_coarse"]))
    if "domain_fine" in truth:
        panels.append(("Ground truth · fine", truth["domain_fine"]))
    for name, lab in labels_by_level.items():
        panels.append((f"Baseline · {name.replace('_', ' ')}", lab))

    n = len(panels)
    fig, axes = plt.subplots(1, n, figsize=(3.0 * n, 3.1), squeeze=False)
    for ax, (title, lab) in zip(axes[0], panels):
        cats = sorted(set(map(str, lab)))
        colors = pal(len(cats))
        cmap = {c: colors[i % len(colors)] for i, c in enumerate(cats)}
        ax.scatter(coords[:, 0], coords[:, 1], s=7, linewidths=0,
                   c=[cmap[str(v)] for v in lab])
        ax.set_title(f"{title}\n({len(cats)} domains)", fontsize=8)
        ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
        for s in ax.spines.values():
            s.set_visible(False)
    fig.tight_layout()
    save_fig(fig, out)


def fig_spatial_gain(gain_df: pd.DataFrame, out: str):
    """Dumbbell:每个分辨率下 无空间对照 → 空间平滑基线 的真值 ARI 变化(不用条形图)。"""
    fig, ax = plt.subplots(figsize=(5.6, 0.24 * len(gain_df) + 1.6))
    y = np.arange(len(gain_df))
    ax.hlines(y, gain_df["ari_nonspatial"], gain_df["ari_spatial"],
              color="#BBBBBB", linewidth=1.6, zorder=1)
    ax.scatter(gain_df["ari_spatial"], y, s=34, color="#E64B35",
               label="spatially smoothed baseline", zorder=2)
    # 对照画成空心圈,两端重合时也看得见(低分辨率下两者常常一样)
    ax.scatter(gain_df["ari_nonspatial"], y, s=42, facecolors="none",
               edgecolors="#4DBBD5", linewidths=1.6,
               label="non-spatial control", zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels([f"res {r}" for r in gain_df.index], fontsize=7)
    ax.set_xlabel("ARI vs. ground-truth fine domains")
    ax.set_title("What spatial context buys, per resolution")
    ax.legend(frameon=False, fontsize=7, loc="lower right")
    fig.tight_layout()
    save_fig(fig, out)


# ---------------------------------------------------------------- 主流程
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", default=EXAMPLE, help="csv:行=spot,列 x,y,[domain_*],基因…")
    ap.add_argument("--outdir", default=RESULTS)
    ap.add_argument("--k-set", default="6,10,16,24", help="空间邻域尺度网格")
    ap.add_argument("--res-min", type=float, default=0.1)
    ap.add_argument("--res-max", type=float, default=1.2)
    ap.add_argument("--res-step", type=float, default=0.1)
    ap.add_argument("--n-repeats", type=int, default=3, help="每格随机种子数(稳定性用)")
    ap.add_argument("--top-frac", type=float, default=0.15, help="取稳定性前多少比例做候选")
    ap.add_argument("--run-scale", action="store_true", help="尝试调用真正的 SCALE(需安装)")
    a = ap.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    os.makedirs(ASSETS, exist_ok=True)
    np.random.seed(SEED)
    set_pub_style()

    print("[575] Step 1 读取输入")
    _, coords, X, genes, truth = load_table(a.input)
    print(f"       {X.shape[0]} spots × {X.shape[1]} genes;真值列 {list(truth)}")

    k_set = [int(v) for v in a.k_set.split(",")]
    res_set = list(np.round(np.arange(a.res_min, a.res_max + 1e-9, a.res_step), 3))

    print("[575] Step 2 基线:空间平滑 × 多分辨率 Leiden × 多种子稳定性网格")
    stab, nclu, store = stability_grid(X, coords, k_set, res_set, n_repeats=a.n_repeats)
    stab.to_csv(os.path.join(a.outdir, "stability_grid.csv"))
    nclu.to_csv(os.path.join(a.outdir, "n_clusters_grid.csv"))

    print("[575] Step 3 挑多尺度层(朴素代理,非 SCALE 熵搜索)")
    picks = pick_levels(stab, nclu, top_frac=a.top_frac)
    labels_by_level = {name: store[key] for name, key in picks.items()}
    for name, (k, r) in picks.items():
        print(f"       {name}: k={k} res={r} stability={stab.loc[k, r]:.3f} "
              f"n_clusters={nclu.loc[k, r]:.0f}")

    print("[575] Step 4 无空间信息对照 + 真值评估")
    from sklearn.metrics import adjusted_rand_score
    ctrl = nonspatial_control(X, res_set)
    rows = []
    fine = truth.get("domain_fine")
    for r in res_set:
        k_best = stab[r].astype(float).idxmax()
        lab_sp = store[(k_best, r)]
        rows.append({
            "resolution": r,
            "k_best": k_best,
            "ari_spatial": adjusted_rand_score(fine, lab_sp) if fine is not None else np.nan,
            "ari_nonspatial": adjusted_rand_score(fine, ctrl[r]) if fine is not None else np.nan,
            "n_clusters_spatial": len(set(lab_sp)),
            "n_clusters_nonspatial": len(set(ctrl[r])),
        })
    gain = pd.DataFrame(rows).set_index("resolution")
    gain.to_csv(os.path.join(a.outdir, "spatial_vs_nonspatial_ari.csv"))
    print(f"       spatial 中位 ARI={gain['ari_spatial'].median():.3f} | "
          f"non-spatial 中位 ARI={gain['ari_nonspatial'].median():.3f}")

    lab_out = pd.DataFrame({n: v for n, v in labels_by_level.items()})
    lab_out.insert(0, "y", coords[:, 1]); lab_out.insert(0, "x", coords[:, 0])
    lab_out.to_csv(os.path.join(a.outdir, "baseline_domain_labels.csv"), index=False)

    print("[575] Step 5 出图")
    for d in (a.outdir, ASSETS):
        fig_stability(stab, picks, os.path.join(d, "575_stability_grid"))
        fig_domains(coords, truth, labels_by_level, os.path.join(d, "575_domains_multiscale"))
        if fine is not None:
            fig_spatial_gain(gain, os.path.join(d, "575_spatial_gain_dumbbell"))

    scale_status = run_scale_guarded(a.input) if a.run_scale else {
        "status": "not requested", "hint": "加 --run-scale 尝试真正的 SCALE"}
    print(f"[575] SCALE 路径: {scale_status['status']}"
          + (f" — {scale_status.get('reason', scale_status.get('note',''))}" if scale_status.get('reason') or scale_status.get('note') else ""))
    if scale_status.get("install"):
        print(f"       安装: {scale_status['install']}")

    summary = {
        "n_spots": int(X.shape[0]), "n_genes": int(X.shape[1]),
        "k_set": k_set, "n_resolutions": len(res_set), "n_repeats": a.n_repeats,
        "picked_levels": {n: {"k": int(k), "resolution": float(r),
                              "stability": float(stab.loc[k, r]),
                              "n_clusters": int(nclu.loc[k, r])}
                          for n, (k, r) in picks.items()},
        "median_ari_spatial": float(gain["ari_spatial"].median()),
        "median_ari_nonspatial": float(gain["ari_nonspatial"].median()),
        "scale_upstream": scale_status,
        "seed": SEED,
    }
    with open(os.path.join(a.outdir, "575_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=1, ensure_ascii=False, default=str)
    print(f"[575] 完成 → {a.outdir}")


if __name__ == "__main__":
    main()
