"""579 · SIMO — 单细胞多组学到空间的概率对齐(spatial integration of multi-omics).

上游方法 SIMO(Yang et al., Nat Commun 2025;16:1265,doi:10.1038/s41467-025-56523-4,
PMID 39893194,repo https://github.com/ZJUFanLab/SIMO)用 fused Gromov-Wasserstein
最优传输把 scRNA 细胞概率性地分配到空间 spot,再用 unbalanced OT 做跨模态标签传递,
把第二个非转录组模态(scATAC / 甲基化)也拉到同一张空间图上。

本模块跑三条路线,结果并排对比:
  A. naive 基线    —— 最大表达相关性贪心分配,不用 OT、不用空间图(本库要求的朴素对照)
  B. OT 参照       —— 直接用 POT 的 fused Gromov-Wasserstein 求解(SIMO 同一算法族)
  C. SIMO 正牌路线 —— 守卫式:未装 simo-omics 就跳过并打印真实安装命令

★ B 是我们自己用 POT 写的传输求解,**不是** SIMO 的 alignment_1_batch;
  SIMO 的自适应权重、aware-label 约束、坐标 refine 与跨模态 UOT 都不在 B 里。
  B 的作用是给"OT 比贪心好多少"一个可跑的量化下界,不是复现 SIMO。
"""
from __future__ import annotations
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
EXAMPLE = os.path.join(HERE, "example_data")
RESULTS = os.path.join(HERE, "results")
ASSETS = os.path.join(HERE, "assets")
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..", "..", "_framework")))

SEED = 2026


# ---------------------------------------------------------------- 数据与预处理
def load_inputs(datadir: str):
    """读 example_data/ 的 5 个 CSV;换自己的数据保持同样列名即可。"""
    st = pd.read_csv(os.path.join(datadir, "st_expression.csv"), index_col=0)
    coords = pd.read_csv(os.path.join(datadir, "st_coords.csv"))
    rna = pd.read_csv(os.path.join(datadir, "sc_rna_expression.csv"), index_col=0)
    rna_meta = pd.read_csv(os.path.join(datadir, "sc_rna_meta.csv"))
    m2 = pd.read_csv(os.path.join(datadir, "sc_mod2_gene_activity.csv"), index_col=0)
    m2_meta = pd.read_csv(os.path.join(datadir, "sc_mod2_meta.csv"))
    genes = [g for g in st.columns if g in rna.columns]        # 只用共享基因
    return st[genes], coords, rna[genes], rna_meta, m2[genes], m2_meta


def _align_meta(meta: pd.DataFrame, key: str, order) -> pd.DataFrame:
    """把 meta 表按表达矩阵的行序对齐。显式 reindex,避免 .loc 把索引名弄丢。"""
    out = meta.set_index(key).reindex(order)
    if out.isna().any().any():
        missing = out[out.isna().any(1)].index.tolist()[:5]
        sys.exit(f"meta 表缺少这些 {key}: {missing} …(表达矩阵与 meta 行不匹配)")
    out.index.name = key
    return out.reset_index()


def _cpm_log(df: pd.DataFrame) -> np.ndarray:
    """CPM + log1p,消掉细胞/spot 间深度差异(否则相关性被文库大小主导)。"""
    X = df.to_numpy(dtype=float)
    tot = X.sum(1, keepdims=True)
    tot[tot == 0] = 1.0
    return np.log1p(X / tot * 1e4)


def _zscore_cols(X: np.ndarray) -> np.ndarray:
    sd = X.std(0, keepdims=True)
    sd[sd == 0] = 1.0
    return (X - X.mean(0, keepdims=True)) / sd


def _corr(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """行×行 Pearson 相关矩阵,返回 (nA, nB)。"""
    Az = _zscore_cols(A.T).T
    Bz = _zscore_cols(B.T).T
    Az /= np.sqrt((Az ** 2).sum(1, keepdims=True)) + 1e-12
    Bz /= np.sqrt((Bz ** 2).sum(1, keepdims=True)) + 1e-12
    return Az @ Bz.T


# ---------------------------------------------------------------- A. naive 基线
def baseline_greedy(corr: np.ndarray) -> np.ndarray:
    """朴素对照:每个细胞直接选相关性最高的 spot。无 OT、无空间结构、无配额约束。

    这是空间映射类方法的地板:多个细胞会挤到同一个 spot,
    而且完全不用空间坐标信息。
    """
    return corr.argmax(1)


# ---------------------------------------------------------------- B. OT 参照
def ot_fgw(corr: np.ndarray, cell_emb: np.ndarray, spot_xy: np.ndarray,
           alpha: float = 0.3, verbose: bool = True):
    """fused Gromov-Wasserstein:表达代价 + 两侧图结构一致性,POT 直接求解。

    M  = 1 - corr           细胞×spot 的表达不相似度(fused 项)
    C1 = 细胞嵌入距离矩阵    源侧结构
    C2 = spot 空间距离矩阵   目标侧结构(这一项让相邻细胞落到相邻 spot)
    alpha 权衡两者,和 SIMO 论文里的 alpha 是同一个含义,但取值未与上游对齐。

    返回 (分配的 spot 下标, 传输计划 T)。装不上 POT 时返回 (None, None)。
    """
    try:
        import ot as pot
    except ImportError:
        if verbose:
            print("       POT 未安装,跳过 OT 参照:pip install POT")
        return None, None

    M = (1.0 - corr).astype(np.float64)
    M /= M.max() + 1e-12
    C1 = _pdist(cell_emb)
    C2 = _pdist(spot_xy)
    C1 /= C1.max() + 1e-12
    C2 /= C2.max() + 1e-12
    p = np.ones(M.shape[0]) / M.shape[0]          # 细胞均匀边缘
    q = np.ones(M.shape[1]) / M.shape[1]          # spot 均匀边缘 -> 隐含配额约束
    T = pot.gromov.fused_gromov_wasserstein(
        M, C1, C2, p, q, loss_fun="square_loss", alpha=alpha, max_iter=200)
    return np.asarray(T).argmax(1), np.asarray(T)


def _pdist(X: np.ndarray) -> np.ndarray:
    from scipy.spatial.distance import cdist
    return cdist(X, X, metric="euclidean")


# ---------------------------------------------------------------- C. SIMO 守卫路径
def run_simo(verbose: bool = True) -> dict:
    """SIMO 正牌路线。未安装则优雅退出并给真实安装命令,绝不静默降级。

    下面这串调用签名**逐条核对自上游源码**(ZJUFanLab/SIMO @ main,函数定义位置见括注),
    不是抄教程 notebook —— 官方 tutorial/mouse_brain.ipynb 只在上游 README 里被链接,
    本地克隆未包含该 notebook,故不声称"读过教程"。
      load_data          simo/helper.py:23
      process_anndata    simo/helper.py:33     (neighbors/umap 默认 True,n_comps 默认 100)
      find_marker        simo/helper.py:138    (return_anndata=False 时返回 gene list,helper.py:202)
      alignment_1_batch  simo/simo.py:1014     (内部循环调 alignment_1[def simo.py:13],
                                                alignment_1 在 simo.py:151 调 fgw_ot[def simo.py:169])
      assign_coord_1     simo/simo.py:357      (layer 默认 'data',top_num 默认 None)
      alignment_2        simo/simo.py:203      (返回 4 元组,simo.py:354;内部 label_transfer 走
                                                ot.unbalanced.mm_unbalanced,helper.py:321)
      assign_coord_2     simo/simo.py:662
      regulation_analysis / spatial_regulation  simo/regulation.py:64 / :312
        (regulation 未被 simo/__init__.py 星号导入,必须写 `from simo.regulation import ...`)
    完整流程(参数名与默认值以上述源码为准,本模块不代跑、不固定默认值):
      from simo import (load_data, process_anndata, find_marker,
                        alignment_1_batch, assign_coord_1, alignment_2, assign_coord_2)
      from simo.regulation import regulation_analysis, spatial_regulation

      st  = process_anndata(st,  neighbors=True, umap=True, n_comps=100)
      rna = process_anndata(rna, neighbors=True, umap=True, n_comps=100)
      gene_list1 = find_marker(st, rna, gene_selection_method='deg', deg_num=100,
                               marker1_by='seurat_clusters', marker2_by='cell_type')
      out1 = alignment_1_batch(adata1=st[:, gene_list1], adata2=rna[:, gene_list1],
                               alpha=0.1, aware_st_label='seurat_clusters',
                               aware_sc_label='cell_type')
      map1 = assign_coord_1(adata1=..., adata2=..., out_data=out1,
                            no_repeated_cells=True, top_num=5, layer='data')
      out2, transfer_df, obs_df1, obs_df2 = alignment_2(
          adata1=rna_mapped[:, gene_list2], adata2=mod2[:, gene_list2],
          coor_df=map1, reg=1, adata1_avg_by='leiden', adata2_avg_by='cell_type',
          modality2_type='neg')
      map2 = assign_coord_2(adata1=st, adata2=mod2, out_data=out2, top_num=3)

    注意:PyPI 包名是 **simo-omics**(setup.py:name='simo-omics'),import 名是 `simo`;
    PyPI 上另有一个完全无关的包 `simo`(自述为 "Smart Home Supremacy",智能家居),
    `pip install simo` 会装错东西。
    """
    try:
        import simo
    except ImportError:
        return {
            "status": "skipped",
            "reason": "simo 未安装",
            "install": "pip install simo-omics   # 注意不是 `pip install simo`(那是无关包)",
            "upstream_repo": "https://github.com/ZJUFanLab/SIMO",
            "tutorial": "上游 README 链接 tutorial/mouse_brain.ipynb 与 tutorial/human_heart.ipynb"
                        "(未随本地克隆取回,未实读)",
        }
    have = [n for n in ("load_data", "process_anndata", "find_marker", "alignment_1_batch",
                        "assign_coord_1", "alignment_2", "assign_coord_2")
            if hasattr(simo, n)]
    return {
        "status": "ready",
        "simo_version": getattr(simo, "__version__", "?"),
        "exported_api_found": have,
        "next": "按官方 mouse_brain.ipynb 逐步执行;本模块不代跑,以免固化未经核对的参数",
    }


# ---------------------------------------------------------------- 打分
def score(assign: np.ndarray, cell_layers: np.ndarray, spot_layers: np.ndarray,
          spot_xy: np.ndarray) -> dict:
    """ground-truth 打分:分配到的 spot 层标签是否等于细胞真实层。

    合成数据里层标签就是真答案,所以这个准确率是可信的外部指标,
    不是方法自己算出来的内部一致性。
    """
    got = spot_layers[assign]
    ok = got == cell_layers
    # 空间偏移:分配位置 与 该细胞真实层所有 spot 的质心 的距离
    disp = []
    for lay in np.unique(cell_layers):
        cen = spot_xy[spot_layers == lay].mean(0)
        m = cell_layers == lay
        disp.append(np.linalg.norm(spot_xy[assign[m]] - cen, axis=1))
    counts = np.bincount(assign, minlength=len(spot_layers))
    return {
        "layer_accuracy": round(float(ok.mean()), 4),
        "median_displacement": round(float(np.median(np.concatenate(disp))), 3),
        "n_spots_used": int((counts > 0).sum()),
        # 占用率/最大堆叠:贪心会把大量细胞堆到少数 spot,
        # 单看准确率看不出这个问题,而它直接毁掉下游的空间邻域分析
        "spot_occupancy_frac": round(float((counts > 0).mean()), 4),
        "max_cells_per_spot": int(counts.max()),
        "per_layer": {str(l): round(float(ok[cell_layers == l].mean()), 4)
                      for l in np.unique(cell_layers)},
    }


# ---------------------------------------------------------------- 出图
def make_figures(coords, rna_meta, assign_g, assign_ot, T, sc_g, sc_ot,
                 m2_meta, m2_g_s, m2_ot_s, layers, outdir):
    """出 4 张图。sc_*/m2_*_s 是 score() 返回的打分字典,不是分配向量。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from pubstyle import set_pub_style, pal, save_fig

    set_pub_style(base_size=10)
    cols = dict(zip(layers, pal(len(layers), "npg")))
    xy = coords[["x", "y"]].to_numpy(float)

    # --- 图1:空间映射三联(scatter,不用条形图) ---
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6))
    panels = [("Ground-truth spots", coords["layer"].to_numpy(), xy, 26),
              ("A. Greedy correlation", rna_meta["layer"].to_numpy(), xy[assign_g], 14),
              ("B. Fused GW optimal transport", rna_meta["layer"].to_numpy(), xy[assign_ot], 14)]
    rng = np.random.default_rng(SEED)
    for ax, (title, labs, pos, s) in zip(axes, panels):
        jit = rng.normal(0, 0.16, pos.shape) if s < 20 else 0.0   # 抖动避免点重叠
        for lay in layers:
            m = labs == lay
            ax.scatter((pos + jit)[m, 0], (pos + jit)[m, 1], s=s, c=cols[lay],
                       label=lay, edgecolors="none", alpha=0.85)
        ax.set_title(title)
        ax.set_xlabel("Spatial x"); ax.set_ylabel("Spatial y")
        ax.set_aspect("equal")
    axes[-1].legend(title="Layer", loc="center left", bbox_to_anchor=(1.02, 0.5), markerscale=1.6)
    fig.tight_layout()
    save_fig(fig, os.path.join(outdir, "fig1_spatial_mapping"))
    plt.close(fig)

    # --- 图2:每层准确率 slopegraph(greedy -> OT),两个模态 ---
    fig, axes = plt.subplots(1, 2, figsize=(8, 4.2), sharey=True)
    for ax, (gs, os_, ttl) in zip(axes, [(sc_g, sc_ot, "Modality 1: scRNA-seq"),
                                         (m2_g_s, m2_ot_s, "Modality 2: gene activity")]):
        placed = []                                   # 记录已放标签的 y,避免重叠
        for lay in layers:
            a, b = gs["per_layer"][lay], os_["per_layer"][lay]
            ax.plot([0, 1], [a, b], "-o", color=cols[lay], lw=2, ms=7, label=lay)
            span = max(1e-6, max(os_["per_layer"].values()) - min(os_["per_layer"].values()))
            ty = b
            while any(abs(ty - p) < span * 0.09 for p in placed):
                ty += span * 0.09
            placed.append(ty)
            ax.annotate(lay, (1.03, ty), fontsize=8, va="center", color=cols[lay])
        ax.plot([0, 1], [gs["layer_accuracy"], os_["layer_accuracy"]], "--o",
                color="0.35", lw=2.4, ms=8, zorder=5)
        ax.annotate("overall", (1.03, os_["layer_accuracy"]), fontsize=8,
                    va="center", color="0.35", fontweight="bold")
        ax.set_xticks([0, 1]); ax.set_xticklabels(["Greedy\n(baseline)", "Fused GW\n(OT)"])
        ax.set_xlim(-0.25, 1.45); ax.set_title(ttl)
    axes[0].set_ylabel("Layer assignment accuracy")
    fig.suptitle("Per-layer accuracy: naive baseline vs optimal transport", y=1.02)
    fig.tight_layout()
    save_fig(fig, os.path.join(outdir, "fig2_accuracy_slopegraph"))
    plt.close(fig)

    # --- 图4:空间占用/堆叠 dumbbell(OT 的边缘约束到底买到了什么) ---
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.2))
    rows = [("Modality 1\n(scRNA)", sc_g, sc_ot), ("Modality 2\n(gene activity)", m2_g_s, m2_ot_s)]
    for ax, (metric, lab) in zip(axes, [("spot_occupancy_frac", "Fraction of spots occupied"),
                                        ("max_cells_per_spot", "Max cells piled on one spot")]):
        for i, (name, g, o) in enumerate(rows):
            ax.plot([g[metric], o[metric]], [i, i], "-", color="0.75", lw=3, zorder=1)
            ax.scatter(g[metric], i, s=95, color=pal(2, "npg")[0], zorder=3,
                       label="Greedy (baseline)" if i == 0 else None)
            ax.scatter(o[metric], i, s=95, color=pal(2, "npg")[1], zorder=3,
                       label="Fused GW (OT)" if i == 0 else None)
        ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows])
        ax.set_xlabel(lab); ax.set_ylim(-0.6, len(rows) - 0.4)
    axes[0].legend(loc="upper center", bbox_to_anchor=(1.05, 1.32), ncol=2)
    fig.tight_layout()
    save_fig(fig, os.path.join(outdir, "fig4_spot_occupancy_dumbbell"))
    plt.close(fig)

    # --- 图3:传输质量热图(细胞真实层 × spot 层) ---
    if T is not None:
        cl = rna_meta["layer"].to_numpy()
        sl = coords["layer"].to_numpy()
        Mm = np.zeros((len(layers), len(layers)))
        for i, a in enumerate(layers):
            for j, b in enumerate(layers):
                Mm[i, j] = T[np.ix_(cl == a, sl == b)].sum()
        Mm /= Mm.sum(1, keepdims=True)
        fig, ax = plt.subplots(figsize=(4.6, 4.0))
        im = ax.imshow(Mm, cmap="viridis", vmin=0, vmax=Mm.max())
        ax.set_xticks(range(len(layers))); ax.set_xticklabels(layers, rotation=45, ha="right")
        ax.set_yticks(range(len(layers))); ax.set_yticklabels(layers)
        ax.set_xlabel("Spot layer (target)"); ax.set_ylabel("Cell layer (source)")
        ax.set_title("Transport mass, row-normalised")
        for i in range(len(layers)):
            for j in range(len(layers)):
                ax.text(j, i, f"{Mm[i, j]:.2f}", ha="center", va="center", fontsize=8,
                        color="white" if Mm[i, j] < Mm.max() * 0.6 else "black")
        fig.colorbar(im, ax=ax, shrink=0.8, label="Fraction of mass")
        fig.tight_layout()
        save_fig(fig, os.path.join(outdir, "fig3_transport_heatmap"))
        plt.close(fig)


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--datadir", default=EXAMPLE, help="输入目录(5 个 CSV)")
    ap.add_argument("--outdir", default=RESULTS)
    ap.add_argument("--alpha", type=float, default=0.3, help="FGW 中结构项权重")
    ap.add_argument("--run-simo", action="store_true", help="尝试 SIMO 正牌路线(需 pip install simo-omics)")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    os.makedirs(ASSETS, exist_ok=True)
    np.random.seed(SEED)

    print("[579] Step 1 读入数据")
    st, coords, rna, rna_meta, m2, m2_meta = load_inputs(a.datadir)
    coords = _align_meta(coords, "spot", st.index)
    rna_meta = _align_meta(rna_meta, "cell", rna.index)
    m2_meta = _align_meta(m2_meta, "cell", m2.index)
    layers = sorted(coords["layer"].unique())
    print(f"       spots={st.shape[0]}  rna={rna.shape[0]}  mod2={m2.shape[0]}  genes={st.shape[1]}")

    Xst, Xrna = _cpm_log(st), _cpm_log(rna)
    Xm2 = m2.to_numpy(float)                     # 第二模态已是 [0,1] 活性分,不做 CPM
    spot_xy = coords[["x", "y"]].to_numpy(float)
    spot_lay = coords["layer"].to_numpy()

    from sklearn.decomposition import PCA
    emb_rna = PCA(n_components=15, random_state=SEED).fit_transform(Xrna)
    emb_m2 = PCA(n_components=15, random_state=SEED).fit_transform(Xm2)

    out = {}
    print("[579] Step 2 模态1 scRNA:A 贪心基线 vs B 融合 GW 最优传输")
    corr = _corr(Xrna, Xst)
    ag = baseline_greedy(corr)
    sc_g = score(ag, rna_meta["layer"].to_numpy(), spot_lay, spot_xy)
    print(f"       A greedy : acc={sc_g['layer_accuracy']}  spots_used={sc_g['n_spots_used']}")
    ao, T = ot_fgw(corr, emb_rna, spot_xy, alpha=a.alpha)
    if ao is None:
        sys.exit("POT 不可用,B 路线无法运行(A 已完成)。pip install POT")
    sc_ot = score(ao, rna_meta["layer"].to_numpy(), spot_lay, spot_xy)
    print(f"       B FGW-OT : acc={sc_ot['layer_accuracy']}  spots_used={sc_ot['n_spots_used']}")
    out["modality1_greedy"], out["modality1_fgw_ot"] = sc_g, sc_ot

    print("[579] Step 3 模态2(非转录组):同样两条路线")
    # 第二模态与表达负相关 -> 取负号后再算相关,对应 SIMO 的 modality2_type='neg'
    corr2 = _corr(-Xm2, Xst)
    mg = baseline_greedy(corr2)
    m2_g = score(mg, m2_meta["layer"].to_numpy(), spot_lay, spot_xy)
    mo, _ = ot_fgw(corr2, emb_m2, spot_xy, alpha=a.alpha, verbose=False)
    m2_ot = score(mo, m2_meta["layer"].to_numpy(), spot_lay, spot_xy)
    print(f"       A greedy : acc={m2_g['layer_accuracy']}   B FGW-OT : acc={m2_ot['layer_accuracy']}")
    out["modality2_greedy"], out["modality2_fgw_ot"] = m2_g, m2_ot

    print("[579] Step 4 SIMO 正牌路线")
    out["simo"] = run_simo() if a.run_simo else {"status": "not requested (--run-simo)"}
    for k, v in out["simo"].items():
        print(f"       {k}: {v}")

    print("[579] Step 5 出图 + 落盘")
    pd.DataFrame({"cell": rna_meta["cell"], "true_layer": rna_meta["layer"],
                  "spot_greedy": coords["spot"].to_numpy()[ag],
                  "spot_fgw_ot": coords["spot"].to_numpy()[ao],
                  "x_fgw_ot": spot_xy[ao, 0], "y_fgw_ot": spot_xy[ao, 1]}
                 ).to_csv(os.path.join(a.outdir, "modality1_cell_to_spot.csv"), index=False)
    pd.DataFrame({"cell": m2_meta["cell"], "true_layer": m2_meta["layer"],
                  "spot_greedy": coords["spot"].to_numpy()[mg],
                  "spot_fgw_ot": coords["spot"].to_numpy()[mo]}
                 ).to_csv(os.path.join(a.outdir, "modality2_cell_to_spot.csv"), index=False)

    make_figures(coords, rna_meta, ag, ao, T, sc_g, sc_ot, m2_meta, m2_g, m2_ot,
                 layers, a.outdir)
    for f in ("fig1_spatial_mapping", "fig2_accuracy_slopegraph",
              "fig3_transport_heatmap", "fig4_spot_occupancy_dumbbell"):
        src = os.path.join(a.outdir, f + ".png")
        if os.path.exists(src):
            import shutil
            shutil.copyfile(src, os.path.join(ASSETS, f + ".png"))

    # 依赖快照(铁律6:结果要可复现,版本必须随结果落盘)
    import importlib
    out["_env"] = {"python": sys.version.split()[0], "seed": SEED, "alpha": a.alpha}
    for m in ("numpy", "pandas", "scipy", "sklearn", "ot", "matplotlib"):
        try:
            out["_env"][m] = getattr(importlib.import_module(m), "__version__", "?")
        except ImportError:
            out["_env"][m] = "not installed"

    with open(os.path.join(a.outdir, "579_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1, ensure_ascii=False)
    print(f"[579] done -> {a.outdir}")


if __name__ == "__main__":
    main()
