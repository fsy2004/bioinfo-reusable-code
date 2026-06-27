# ============================================================================
# 模块 543 · squidpy 空间统计 (squidpy spatial statistics)
# ----------------------------------------------------------------------------
# 做什么:
#   用 squidpy 对空间转录组(spot/cell × 坐标)做一站式空间统计工具箱——
#   空间邻接图 → 邻域富集(谁挨着谁) → 空间自相关 Moran's I(哪些基因有
#   空间结构) → 细胞类型共现概率(距离-依赖) → Ripley's L(聚集/分散)。
#
# ★诚实基线(本模块的灵魂):
#   Moran's I 容易"看起来都显著"。本脚本对每个基因额外计算一个【坐标置换
#   null】:把 spatial 坐标随机打乱、保持表达不变,重算 Moran's I。真实空间
#   结构化基因的 observed I 应远高于 null(I→0);若 observed≈null 则信号是
#   偶然。脚本输出 observed-vs-null 散点 + 经验 p,证明信号非随机产物。
#
# Turnkey 用法 (零改动即跑,自动合成 ~500 spot 演示数据):
#   python 543_squidpy_spatial_statistics.py
#
# 换数据用法 (传你自己的 .h5ad,需含 adata.obsm['spatial'] 和一个类型标签列):
#   python 543_squidpy_spatial_statistics.py --input my.h5ad --cluster_key cell_type
#
# 真实 API (squidpy 1.8.2 实测核对):
#   sq.gr.spatial_neighbors(adata, coord_type='generic', n_neighs=6)
#   sq.gr.nhood_enrichment(adata, cluster_key, seed=, n_perms=)  -> uns['<k>_nhood_enrichment']['zscore']
#   sq.gr.spatial_autocorr(adata, mode='moran', n_perms=, seed=) -> uns['moranI'] (cols: I, pval_sim, ...)
#   sq.gr.co_occurrence(adata, cluster_key, interval=) -> uns['<k>_co_occurrence']['occ','interval']
#   sq.gr.ripley(adata, cluster_key, mode='L', n_simulations=, seed=) -> uns['<k>_ripley_L']
# ============================================================================

import sys
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# --- 复用顶刊绘图框架 (照抄定位逻辑) -----------------------------------------
HERE = Path(__file__).resolve().parent
for up in [HERE, *HERE.parents]:
    if (up / "_framework" / "pubstyle.py").exists():
        sys.path.insert(0, str(up / "_framework")); break
try:
    from pubstyle import (set_pub_style, save_fig, pal, panel_labels,
                          CMAP_CONT, CMAP_DIVERGE, NATURE_W1, NATURE_W2)
except Exception:
    def set_pub_style(*a, **k): pass
    def save_fig(fig, f, dpi=300):
        from pathlib import Path as _P
        _P(f).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(str(f) + ".pdf"); fig.savefig(str(f) + ".png", dpi=dpi)
    def pal(n=None, name="npg"):
        import matplotlib; return list(matplotlib.cm.tab10.colors)
    def panel_labels(*a, **k): pass
    CMAP_CONT = "viridis"; CMAP_DIVERGE = "RdBu_r"; NATURE_W1 = 3.5; NATURE_W2 = 7.0

SEED = 42

# --- 路径 (全部从 __file__ 派生, 绝不 hardcode / 绝不 chdir) -----------------
DIR_DATA = HERE / "example_data"
DIR_RES = HERE / "results"
DIR_ASSETS = HERE / "assets"
for d in (DIR_DATA, DIR_RES, DIR_ASSETS):
    d.mkdir(parents=True, exist_ok=True)


# ============================================================================
# 1) 合成空间 AnnData: ~500 spot, 真实空间结构信号
# ============================================================================
def make_synthetic(path: Path, n_per_side: int = 23):
    """生成一个带真实空间结构的小 AnnData 并保存为 .h5ad。

    设计:
      - spot 排成略带抖动的网格 (~n^2 个)。
      - 4 个空间域 (domain) 按 2D 高斯 niche 分布 -> 邻域富集非随机。
      - 部分基因表达由域/坐标驱动 (空间结构化) -> Moran's I 应显著;
        另一部分纯噪声 -> Moran's I 应≈0, 作为内部阴性对照。
    """
    import anndata as ad
    rng = np.random.default_rng(SEED)

    side = np.linspace(0, 30, n_per_side)
    gx, gy = np.meshgrid(side, side)
    coord = np.column_stack([gx.ravel(), gy.ravel()]).astype(float)
    coord += rng.normal(0, 0.35, coord.shape)   # 抖动, 显式种子
    n = coord.shape[0]

    # --- 4 个域: 每个域一个 niche 中心, spot 归属到最近的软分配中心 ---
    centers = np.array([[8, 8], [22, 8], [8, 22], [22, 22]], float)
    dom_names = ["Tumor", "Stroma", "Immune", "Vessel"]
    d2 = ((coord[:, None, :] - centers[None, :, :]) ** 2).sum(-1)
    # 软分配 + 噪声, 让边界自然
    logits = -d2 / 40.0 + rng.normal(0, 0.6, d2.shape)
    domain_idx = logits.argmax(1)
    domain = np.array(dom_names)[domain_idx]

    # --- 基因 ---
    n_struct = 8     # 空间结构化 marker (每域 2 个)
    n_grad = 2       # 沿坐标梯度
    n_noise = 30     # 纯噪声 (阴性对照)
    genes = []
    cols = []

    # 域 marker: 在所属域高表达
    for di, dn in enumerate(dom_names):
        for k in range(2):
            base = rng.poisson(1.0, n).astype(float)
            base[domain_idx == di] += rng.normal(6, 1.0, (domain_idx == di).sum())
            base = np.clip(base, 0, None)
            cols.append(base); genes.append(f"{dn}_mk{k+1}")

    # 坐标梯度基因
    for k in range(n_grad):
        axis = coord[:, k % 2]
        base = 0.4 * axis + rng.normal(0, 1.5, n)
        base = np.clip(base, 0, None)
        cols.append(base); genes.append(f"Gradient{k+1}")

    # 噪声基因
    for k in range(n_noise):
        cols.append(rng.poisson(2.0, n).astype(float)); genes.append(f"Noise{k+1}")

    X = np.column_stack(cols)
    adata = ad.AnnData(X=X.astype(np.float32))
    adata.var_names = genes
    adata.obs_names = [f"spot{i:04d}" for i in range(n)]
    adata.obsm["spatial"] = coord
    adata.obs["domain"] = pd.Categorical(domain, categories=dom_names)
    # 标记哪些基因是"真信号", 供基线评估
    adata.var["is_structured"] = [g.endswith(("mk1", "mk2")) or g.startswith("Gradient")
                                  for g in genes]
    adata.write_h5ad(path)
    print(f"[synth] 合成 {n} spot × {len(genes)} 基因 -> {path.name} "
          f"(结构化 {int(adata.var['is_structured'].sum())} / 噪声 {n_noise})")
    return adata


# ============================================================================
# 主流程
# ============================================================================
def main():
    ap = argparse.ArgumentParser(description="squidpy 空间统计工具箱 (turnkey)")
    ap.add_argument("--input", type=str, default=None,
                    help=".h5ad, 需含 obsm['spatial'] + 类型标签列; 缺省用合成数据")
    ap.add_argument("--cluster_key", type=str, default="domain",
                    help="adata.obs 中的离散类型/域列名")
    ap.add_argument("--n_neighs", type=int, default=6, help="KNN 空间邻居数")
    ap.add_argument("--n_perms", type=int, default=200, help="置换次数 (Moran/nhood)")
    args = ap.parse_args()

    import squidpy as sq
    import anndata as ad

    set_pub_style(base_size=11)

    # --- 载入数据 ---
    if args.input:
        print(f"[load] 读入 {args.input}")
        adata = ad.read_h5ad(args.input)
    else:
        h5 = DIR_DATA / "spatial_demo.h5ad"
        if not h5.exists():
            make_synthetic(h5)
        adata = ad.read_h5ad(h5)
    ckey = args.cluster_key
    assert "spatial" in adata.obsm, "adata.obsm['spatial'] 缺失"
    assert ckey in adata.obs, f"adata.obs 无列 '{ckey}'"
    if not isinstance(adata.obs[ckey].dtype, pd.CategoricalDtype):
        adata.obs[ckey] = adata.obs[ckey].astype("category")
    coord = np.asarray(adata.obsm["spatial"])
    print(f"[data] {adata.n_obs} spot × {adata.n_vars} 基因; "
          f"'{ckey}' 类别 = {list(adata.obs[ckey].cat.categories)}")

    # ------------------------------------------------------------------
    # step 1: 空间邻接图
    # ------------------------------------------------------------------
    print("[step] 1/6 spatial_neighbors (KNN 空间图)")
    sq.gr.spatial_neighbors(adata, coord_type="generic", n_neighs=args.n_neighs)

    # ------------------------------------------------------------------
    # step 2: 邻域富集 (z-score)
    # ------------------------------------------------------------------
    print("[step] 2/6 nhood_enrichment (邻域富集 z-score)")
    sq.gr.nhood_enrichment(adata, ckey, seed=SEED, n_perms=args.n_perms,
                           show_progress_bar=False)
    z = adata.uns[f"{ckey}_nhood_enrichment"]["zscore"]
    cats = list(adata.obs[ckey].cat.categories)
    z_df = pd.DataFrame(z, index=cats, columns=cats)
    z_df.to_csv(DIR_RES / "nhood_enrichment_zscore.csv")

    # ------------------------------------------------------------------
    # step 3: Moran's I (空间自相关) + ★诚实基线: 坐标置换 null
    # ------------------------------------------------------------------
    print("[step] 3/6 spatial_autocorr Moran's I (+ 坐标置换 null 基线)")
    sq.gr.spatial_autocorr(adata, mode="moran", n_perms=args.n_perms, seed=SEED)
    moran = adata.uns["moranI"].copy()
    moran["is_structured"] = adata.var["is_structured"].reindex(moran.index).values

    # ★ 诚实基线: 打乱坐标 -> 空间结构被破坏 -> Moran's I 应趋近 0
    rng = np.random.default_rng(SEED)
    null_I = {}
    n_null = 30
    perm_adata = adata.copy()
    for g in moran.index:
        vals = []
        gi = np.asarray(perm_adata[:, g].X).ravel()
        for _ in range(n_null):
            order = rng.permutation(perm_adata.n_obs)
            tmp = adata.copy()
            tmp.obsm["spatial"] = coord[order]
            sq.gr.spatial_neighbors(tmp, coord_type="generic", n_neighs=args.n_neighs)
            sq.gr.spatial_autocorr(tmp, mode="moran", genes=[g], n_perms=None, seed=SEED)
            vals.append(float(tmp.uns["moranI"].loc[g, "I"]))
        null_I[g] = vals
    null_mean = {g: np.mean(v) for g, v in null_I.items()}
    null_sd = {g: np.std(v) + 1e-9 for g, v in null_I.items()}
    moran["null_I_mean"] = [null_mean[g] for g in moran.index]
    # 经验单尾 p: observed 比多少比例的 null 还小 (信号强 -> p 小)
    moran["perm_p"] = [(np.sum(np.asarray(null_I[g]) >= moran.loc[g, "I"]) + 1) / (n_null + 1)
                       for g in moran.index]
    moran.to_csv(DIR_RES / "moranI_with_null.csv")
    n_sig_struct = int(((moran["perm_p"] < 0.05) & (moran["is_structured"])).sum())
    n_struct = int(moran["is_structured"].sum())
    print(f"[base] 结构化基因经置换检验显著: {n_sig_struct}/{n_struct} "
          f"(噪声基因 observed I 应≈null)")

    # ------------------------------------------------------------------
    # step 4: 共现概率 (距离依赖)
    # ------------------------------------------------------------------
    print("[step] 4/6 co_occurrence (共现概率曲线)")
    sq.gr.co_occurrence(adata, ckey, interval=20)
    occ = adata.uns[f"{ckey}_co_occurrence"]["occ"]      # (n_cat, n_cat, n_int)
    occ_int = adata.uns[f"{ckey}_co_occurrence"]["interval"]

    # ------------------------------------------------------------------
    # step 5: Ripley's L
    # ------------------------------------------------------------------
    print("[step] 5/6 ripley (Ripley's L 聚集/分散)")
    sq.gr.ripley(adata, ckey, mode="L", n_simulations=50, seed=SEED)
    rip = adata.uns[f"{ckey}_ripley_L"]
    rip_stat = rip["L_stat"]    # cols: bins, <ckey>, stats

    # ------------------------------------------------------------------
    # step 6: 绘图 (非条形图)
    # ------------------------------------------------------------------
    print("[step] 6/6 绘图 -> assets/")
    colors = pal(len(cats))
    cat_color = dict(zip(cats, colors))

    # --- Fig A: 邻域富集 z-score heatmap ---
    figA, axA = plt.subplots(figsize=(NATURE_W1, NATURE_W1 * 0.92))
    vmax = np.nanmax(np.abs(z_df.values))
    im = axA.imshow(z_df.values, cmap=CMAP_DIVERGE, vmin=-vmax, vmax=vmax)
    axA.set_xticks(range(len(cats))); axA.set_xticklabels(cats, rotation=45, ha="right")
    axA.set_yticks(range(len(cats))); axA.set_yticklabels(cats)
    for i in range(len(cats)):
        for j in range(len(cats)):
            axA.text(j, i, f"{z_df.values[i, j]:.0f}", ha="center", va="center",
                     fontsize=7, color="black")
    axA.set_title("Neighborhood enrichment (z-score)")
    figA.colorbar(im, ax=axA, fraction=0.046, pad=0.04, label="z-score")
    figA.tight_layout()
    save_fig(figA, str(DIR_ASSETS / "A_nhood_enrichment_heatmap"))
    plt.close(figA)

    # --- Fig B: Moran's I 排序 lollipop ---
    md = moran.sort_values("I", ascending=True)
    figB, axB = plt.subplots(figsize=(NATURE_W1, NATURE_W1 * 1.15))
    ypos = np.arange(len(md))
    pt_col = ["#C0392B" if s else "#95A5A6" for s in md["is_structured"]]
    axB.hlines(ypos, 0, md["I"], color="#BDC3C7", lw=1.0, zorder=1)
    axB.scatter(md["I"], ypos, c=pt_col, s=22, zorder=2, edgecolor="white", linewidth=0.4)
    axB.set_yticks(ypos); axB.set_yticklabels(md.index, fontsize=5.5)
    axB.axvline(0, color="black", lw=0.6)
    axB.set_xlabel("Moran's I"); axB.set_title("Spatial autocorrelation (ranked)")
    from matplotlib.lines import Line2D
    axB.legend(handles=[Line2D([0], [0], marker="o", color="w", markerfacecolor="#C0392B",
                               markersize=6, label="structured"),
                        Line2D([0], [0], marker="o", color="w", markerfacecolor="#95A5A6",
                               markersize=6, label="noise")],
               fontsize=6, loc="lower right", frameon=False)
    figB.tight_layout()
    save_fig(figB, str(DIR_ASSETS / "B_moranI_lollipop"))
    plt.close(figB)

    # --- Fig C: ★诚实基线 observed Moran's I vs 坐标置换 null ---
    figC, axC = plt.subplots(figsize=(NATURE_W1, NATURE_W1 * 0.92))
    for struct, lab, col in [(True, "structured", "#C0392B"), (False, "noise", "#95A5A6")]:
        sub = moran[moran["is_structured"] == struct]
        axC.scatter(sub["null_I_mean"], sub["I"], c=col, s=26, alpha=0.85,
                    edgecolor="white", linewidth=0.4, label=lab)
    lim = [min(moran["null_I_mean"].min(), moran["I"].min()) - 0.05,
           moran["I"].max() + 0.05]
    axC.plot(lim, lim, "--", color="black", lw=0.8, label="y = x (no signal)")
    axC.set_xlim(lim); axC.set_ylim(lim)
    axC.set_xlabel("Null Moran's I (shuffled coordinates)")
    axC.set_ylabel("Observed Moran's I")
    axC.set_title("Honest baseline: signal vs permutation null")
    axC.legend(fontsize=6, frameon=False, loc="upper left")
    figC.tight_layout()
    save_fig(figC, str(DIR_ASSETS / "C_moran_null_baseline"))
    plt.close(figC)

    # --- Fig D: 共现概率曲线 (以第一个域为 condition) ---
    figD, axD = plt.subplots(figsize=(NATURE_W1, NATURE_W1 * 0.85))
    cond_i = 0
    mids = (occ_int[:-1] + occ_int[1:]) / 2.0
    for j, cj in enumerate(cats):
        axD.plot(mids, occ[cond_i, j, :], "-", color=cat_color[cj], lw=1.6, label=cj)
    axD.axhline(1.0, color="black", lw=0.6, ls=":")
    axD.set_xlabel("Distance"); axD.set_ylabel(f"P(x | {cats[cond_i]}) / P(x)")
    axD.set_title(f"Co-occurrence around '{cats[cond_i]}'")
    axD.legend(fontsize=6, frameon=False, title=None)
    figD.tight_layout()
    save_fig(figD, str(DIR_ASSETS / "D_co_occurrence_curve"))
    plt.close(figD)

    # --- Fig E: Ripley's L 曲线 ---
    figE, axE = plt.subplots(figsize=(NATURE_W1, NATURE_W1 * 0.85))
    for cj in cats:
        sub = rip_stat[rip_stat[ckey] == cj]
        axE.plot(sub["bins"], sub["stats"], "-", color=cat_color[cj], lw=1.6, label=cj)
    axE.set_xlabel("Distance r"); axE.set_ylabel("Ripley's L(r)")
    axE.set_title("Ripley's L (clustering vs dispersion)")
    axE.legend(fontsize=6, frameon=False)
    figE.tight_layout()
    save_fig(figE, str(DIR_ASSETS / "E_ripley_L_curve"))
    plt.close(figE)

    # --- Fig F: 空间散点叠表达 (域 + top Moran 基因) ---
    top_gene = moran.sort_values("I", ascending=False).index[0]
    expr = np.asarray(adata[:, top_gene].X).ravel()
    figF, (axF1, axF2) = plt.subplots(1, 2, figsize=(NATURE_W2, NATURE_W1 * 0.95))
    for cj in cats:
        m = adata.obs[ckey].values == cj
        axF1.scatter(coord[m, 0], coord[m, 1], c=cat_color[cj], s=10, label=cj,
                     edgecolor="none")
    axF1.set_title(f"Spatial domains ('{ckey}')"); axF1.set_aspect("equal")
    axF1.legend(fontsize=6, frameon=False, markerscale=1.5, loc="upper right")
    axF1.set_xlabel("x"); axF1.set_ylabel("y")
    sc2 = axF2.scatter(coord[:, 0], coord[:, 1], c=expr, cmap=CMAP_CONT, s=12,
                       edgecolor="none")
    axF2.set_title(f"Top spatial gene: {top_gene}"); axF2.set_aspect("equal")
    axF2.set_xlabel("x"); axF2.set_ylabel("y")
    figF.colorbar(sc2, ax=axF2, fraction=0.046, pad=0.04, label="expression")
    panel_labels([axF1, axF2])
    figF.tight_layout()
    save_fig(figF, str(DIR_ASSETS / "F_spatial_scatter_expr"))
    plt.close(figF)

    # ------------------------------------------------------------------
    # 版本记录
    # ------------------------------------------------------------------
    # 记录依赖版本快照 (session_info 等价物, 铁律6)
    import importlib.metadata as ilm
    with open(DIR_RES / "versions.txt", "w", encoding="utf-8") as fh:
        for p in ["squidpy", "anndata", "scanpy", "numpy", "pandas", "matplotlib"]:
            try:
                fh.write(f"{p}=={ilm.version(p)}\n")
            except Exception:
                fh.write(f"{p}==NA\n")
    print("[done] 图 -> assets/  表 -> results/  版本 -> results/versions.txt")


if __name__ == "__main__":
    main()
