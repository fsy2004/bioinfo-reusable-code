# -*- coding: utf-8 -*-
# =============================================================================
# 531 · LIANA+ 多方法共识细胞-细胞通讯 (consensus cell-cell communication, CCI)
# -----------------------------------------------------------------------------
# 分类:  16_spatial_communication
# 语言:  Python (scanpy + LIANA+ 1.x 底座, 纯 CPU 可跑通)
#
# 用途:
#   用 LIANA+ 的 rank-aggregate (RRA) 把 6 个配体-受体 (L-R) 打分方法
#     CellPhoneDB / CellChat / Connectome / NATMI / log2FC / SingleCellSignalR
#   统一成一份 *共识* 通讯排名 (specificity_rank × magnitude_rank), 不押单一工具。
#
# ★诚实基线 (本模块的灵魂, 别只报好看指标):
#   单方法的 L-R 排名彼此分歧很大 —— "押单一工具"会系统性漏检。这里内置两道对照:
#     (1) 跨方法一致性 heatmap: 直接展示各单方法之间 Spearman 排名相关有多低 ——
#         分歧真实存在, 所以才需要共识 (诚实地暴露, 不掩盖)。
#     (2) 共识 vs 单一方法 top-K Jaccard 重叠 (lollipop): 量化"若只信一个方法,
#         你会和共识差多少", 以 % 报出, 不夸大共识的稳健性。
#
# 依赖:  liana>=1.7  scanpy  anndata  numpy  pandas  scipy  matplotlib  plotnine
#         conda/pip 安装见 README;本模块在 liana 1.7.3 实跑验证。
#
# Turnkey 用法:  python 531_liana_consensus_cci.py        (合成数据→results/+assets/)
# 换数据用法:    python 531_liana_consensus_cci.py --input my.h5ad --groupby cell_type
#                  (--input 为 AnnData .h5ad, 已 normalize_total + log1p;
#                   --groupby 为 obs 中的细胞类型列名)
#
# 输入:  AnnData (.h5ad), .X = log-normalized 表达, obs[groupby] = 细胞类型分组。
# 路径:  全部脚本相对 (ROOT = 脚本所在目录), 无 setwd / 无绝对路径。
# =============================================================================

import sys
from pathlib import Path

# ---- 定位脚本目录 + 载入顶刊绘图框架 (向上搜 _framework) ---------------------
HERE = Path(__file__).resolve().parent
for up in [HERE, *HERE.parents]:
    if (up / "_framework" / "pubstyle.py").exists():
        sys.path.insert(0, str(up / "_framework"))
        break
try:
    from pubstyle import (set_pub_style, save_fig, pal,
                          CMAP_CONT, CMAP_DIVERGE, NATURE_W1, NATURE_W2)
except Exception:  # 框架缺失时最小降级, 不影响真实分析
    def set_pub_style(*a, **k): pass
    def save_fig(fig, f, dpi=300):
        from pathlib import Path as _P
        _P(f).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(str(f) + ".pdf"); fig.savefig(str(f) + ".png", dpi=dpi)
    def pal(n=None, name="npg"):
        import matplotlib; return list(matplotlib.cm.tab10.colors)
    CMAP_CONT = "viridis"; CMAP_DIVERGE = "RdBu_r"; NATURE_W1 = 3.5; NATURE_W2 = 7.0

import argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

SEED = 42

DIR_EX = HERE / "example_data"
DIR_RES = HERE / "results"
DIR_ASSETS = HERE / "assets"
for d in (DIR_EX, DIR_RES, DIR_ASSETS):
    d.mkdir(parents=True, exist_ok=True)


# -----------------------------------------------------------------------------
# 合成示例数据 (synthetic, for demo only):
#   小 AnnData, 4 个细胞类型, 基因符号取自 LIANA "consensus" L-R 资源 (保证资源能匹配)。
#   设计意图 = 让不同打分方法 *真的* 产生分歧, 这样诚实基线才有意义:
#     · 一组 L-R 走"高幅度低特异" (全细胞广泛表达) → magnitude 类方法 (Connectome/
#       SingleCellSignalR) 会偏爱;
#     · 另一组走"低幅度高特异" (只在一对细胞专一表达) → specificity 类方法
#       (CellPhoneDB/NATMI) 会偏爱。
#   于是单方法排名彼此分歧, 共识 (RRA) 才有真实价值。
# -----------------------------------------------------------------------------
def make_synthetic(path_h5ad):
    import scanpy as sc
    import anndata as ad
    import liana as li

    rng = np.random.default_rng(SEED)

    # 从 LIANA consensus 资源取真实 L-R 基因符号 (拆解 complex 的下划线子单元)
    res = li.rs.select_resource("consensus")
    flat = []
    seen = set()
    for col in ("ligand", "receptor"):
        for g in res[col]:
            for sub in str(g).split("_"):
                if sub not in seen:
                    seen.add(sub); flat.append(sub)
    genes = flat[:160]
    g2i = {g: i for i, g in enumerate(genes)}

    n = 400
    cts = rng.choice(["Tcell", "Myeloid", "Fibroblast", "Epithelial"], size=n)
    # 低底噪计数, 后面在特定细胞类型上叠加信号
    X = rng.poisson(0.4, size=(n, len(genes))).astype(float)

    # (A) 细胞类型专一程序 → 制造"高特异"的 L-R (specificity 类方法偏爱)
    ct_specific = {
        "Myeloid":    genes[0:25],
        "Fibroblast": genes[25:50],
        "Tcell":      genes[50:75],
        "Epithelial": genes[75:100],
    }
    for ct, gl in ct_specific.items():
        mask = (cts == ct)
        for g in gl:
            X[mask, g2i[g]] += rng.poisson(4.0, mask.sum())   # 仅该型高表达

    # (B) 一批"广泛高表达"基因 (所有细胞都中高水平) → 制造"高幅度低特异"信号,
    #     magnitude 类方法会把它们排很前, specificity 类则不会 → 制造方法分歧。
    broad = genes[100:135]
    for g in broad:
        X[:, g2i[g]] += rng.poisson(3.0, n)                    # 全细胞普遍升高

    A = ad.AnnData(
        X=X.astype(np.float32),
        obs=pd.DataFrame({"cell_type": pd.Categorical(cts)},
                         index=[f"cell_{i}" for i in range(n)]),
        var=pd.DataFrame(index=genes),
    )
    A.layers["counts"] = A.X.copy()
    sc.pp.normalize_total(A, target_sum=1e4)
    sc.pp.log1p(A)
    A.write_h5ad(path_h5ad)
    return A


# -----------------------------------------------------------------------------
# 主流程
# -----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="LIANA+ consensus cell-cell communication (CCI)")
    ap.add_argument("--input", default=None, help="AnnData .h5ad (normalize_total + log1p)")
    ap.add_argument("--groupby", default="cell_type", help="obs 中细胞类型列名")
    ap.add_argument("--expr_prop", type=float, default=0.1, help="最小表达比例阈值")
    ap.add_argument("--n_perms", type=int, default=100, help="置换次数 (示例用 100, 真实建议 1000)")
    ap.add_argument("--top_n", type=int, default=15, help="dotplot/tileplot 取的 top L-R 数")
    ap.add_argument("--outdir", default=str(DIR_RES))
    args = ap.parse_args()

    import scanpy as sc
    import liana as li

    set_pub_style(base_size=11, palette="npg")
    res_dir = Path(args.outdir); res_dir.mkdir(parents=True, exist_ok=True)

    # ---- 1. 载入 / 合成数据 --------------------------------------------------
    if args.input:
        print(f"[step 1] loading user AnnData: {args.input}")
        adata = sc.read_h5ad(args.input)
    else:
        h5 = DIR_EX / "synthetic_lr.h5ad"
        if not h5.exists():
            print("[step 1] generating synthetic AnnData (4 cell types, real L-R genes; demo only) ...")
            adata = make_synthetic(h5)
        else:
            print(f"[step 1] loading cached synthetic AnnData: {h5}")
            adata = sc.read_h5ad(h5)

    groupby = args.groupby
    assert groupby in adata.obs, f"groupby '{groupby}' not in adata.obs"
    print(f"[step 1] {adata.n_obs} cells x {adata.n_vars} genes; "
          f"{adata.obs[groupby].nunique()} cell types: "
          f"{list(adata.obs[groupby].unique())}")

    # ---- 2. 共识 rank-aggregate (RRA over 6 methods) -------------------------
    print("[step 2] LIANA rank_aggregate (RRA consensus over 6 L-R methods) ...")
    li.mt.rank_aggregate(
        adata, groupby=groupby, expr_prop=args.expr_prop,
        use_raw=False, n_perms=args.n_perms, seed=SEED, verbose=False,
    )
    cons = adata.uns["liana_res"].copy()
    cons.to_csv(res_dir / "consensus_liana_res.csv", index=False)
    print(f"[step 2]   consensus table: {len(cons)} (L-R × source-target) records; "
          f"cols={[c for c in cons.columns if 'rank' in c or c in ('source','target')]}")

    # ---- 3. 诚实基线: 逐个跑单方法, 取每方法对 L-R 的排名 --------------------
    print("[step 3] honest baseline — running 6 individual methods ...")
    # (method_fn, score_col, ascending)  ascending=True 表示该分数"越小越好"
    method_specs = {
        "CellPhoneDB":       (li.mt.cellphonedb,        "cellphone_pvals", True),
        "CellChat":          (li.mt.cellchat,           "cellchat_pvals",  True),
        "Connectome":        (li.mt.connectome,         "scaled_weight",   False),
        "NATMI":             (li.mt.natmi,              "spec_weight",     False),
        "log2FC":            (li.mt.logfc,              "lr_logfc",        False),
        "SingleCellSignalR": (li.mt.singlecellsignalr,  "lrscore",         False),
    }
    key_cols = ["source", "target", "ligand_complex", "receptor_complex"]

    def _pairkey(df):
        return (df["source"] + "|" + df["target"] + "|"
                + df["ligand_complex"] + "|" + df["receptor_complex"])

    rank_tbl = pd.DataFrame({"pair": _pairkey(cons).drop_duplicates().values})

    # 共识自身的排名 (specificity_rank: 越小越特异 → rank=1 最好)
    cons_key = cons.copy(); cons_key["pair"] = _pairkey(cons_key)
    cons_rank = cons_key.set_index("pair")["specificity_rank"]
    rank_tbl["Consensus"] = rank_tbl["pair"].map(cons_rank).rank(method="average")

    for name, (fn, col, ascending) in method_specs.items():
        df = fn(adata, groupby=groupby, expr_prop=args.expr_prop,
                use_raw=False, n_perms=args.n_perms, seed=SEED,
                verbose=False, inplace=False)
        df = df.copy(); df["pair"] = _pairkey(df)
        s = df.groupby("pair")[col].mean()                 # 同 pair 去重
        r = s.rank(method="average", ascending=ascending)  # rank=1 = 该方法最看好
        rank_tbl[name] = rank_tbl["pair"].map(r)
        print(f"[step 3]   {name:18s} ranked {r.notna().sum()} L-R pairs (score={col})")

    rank_tbl.to_csv(res_dir / "method_rank_table.csv", index=False)

    # 跨方法 Spearman 一致性矩阵 (诚实对照 1 的数据)
    methods = ["Consensus"] + list(method_specs.keys())
    M = rank_tbl[methods].dropna()
    corr = np.eye(len(methods))
    for i, a in enumerate(methods):
        for j, b in enumerate(methods):
            corr[i, j] = spearmanr(M[a], M[b]).statistic
    corr_df = pd.DataFrame(corr, index=methods, columns=methods)
    corr_df.to_csv(res_dir / "cross_method_spearman.csv")
    # 报告单方法之间 (排除 Consensus 行列与对角) 的最低/平均一致性
    single = list(method_specs.keys())
    off = corr_df.loc[single, single].values
    off_vals = off[~np.eye(len(single), dtype=bool)]
    print(f"[step 3]   cross-method Spearman among singles: "
          f"mean={off_vals.mean():.2f}, min={off_vals.min():.2f} "
          f"(低 → 方法分歧大 → 需要共识)")

    # 共识 vs 单方法 top-K Jaccard 重叠 (诚实对照 2 的数据)
    K = min(30, len(M))
    Mp = rank_tbl[["pair"] + methods].dropna()
    top_cons = set(Mp.sort_values("Consensus").head(K)["pair"])
    overlap = {}
    for name in method_specs:
        top_m = set(Mp.sort_values(name).head(K)["pair"])
        overlap[name] = len(top_cons & top_m) / len(top_cons | top_m)
    ov = pd.Series(overlap).sort_values()
    ov.to_csv(res_dir / "consensus_vs_single_topK_jaccard.csv", header=["jaccard"])
    print(f"[step 3]   top-{K} Jaccard(consensus, single):")
    for k, v in ov.items():
        print(f"             {k:18s} {v:.3f}")

    # =========================================================================
    # 绘图 (全部非平凡条形图: dotplot / 网络 / heatmap / lollipop / tile)
    # =========================================================================
    levels = (list(adata.obs[groupby].cat.categories)
              if hasattr(adata.obs[groupby], "cat")
              else sorted(adata.obs[groupby].unique()))

    # ---- Fig 1: 共识 L-R dotplot (size=specificity, colour=magnitude) -------
    print("[fig 1] consensus L-R dotplot ...")
    try:
        from plotnine import labs
        gg = li.pl.dotplot(
            adata=adata, uns_key="liana_res",
            colour="magnitude_rank", size="specificity_rank",
            inverse_colour=True, inverse_size=True,   # rank 越小越好 → 反转使"好"=大/亮
            top_n=args.top_n, orderby="magnitude_rank", orderby_ascending=True,
            cmap=CMAP_CONT, figure_size=(7.5, 6.5), return_fig=True,
        )
        gg = gg + labs(title="Consensus ligand-receptor interactions (LIANA rank-aggregate)")
        gg.save(str(DIR_ASSETS / "fig1_consensus_dotplot.png"), dpi=300, verbose=False)
        gg.save(str(DIR_ASSETS / "fig1_consensus_dotplot.pdf"), verbose=False)
    except Exception as e:
        print(f"[warn] dotplot failed ({e}); skipping fig1")

    # ---- Fig 2: source→target 通讯网络 (圆形布局有向图, 边=显著 L-R 计数) ----
    print("[fig 2] source-target communication network ...")
    _manual_network(cons, levels)

    # ---- Fig 3: 跨方法一致性 heatmap (★诚实对照 1) -------------------------
    print("[fig 3] cross-method consistency heatmap (honest baseline 1) ...")
    fig3, ax3 = plt.subplots(figsize=(NATURE_W1 + 1.6, NATURE_W1 + 1.3))
    im = ax3.imshow(corr_df.values, cmap=CMAP_DIVERGE, vmin=-1, vmax=1, aspect="equal")
    ax3.set_xticks(range(len(methods))); ax3.set_xticklabels(methods, rotation=45, ha="right", fontsize=8)
    ax3.set_yticks(range(len(methods))); ax3.set_yticklabels(methods, fontsize=8)
    for i in range(len(methods)):
        for j in range(len(methods)):
            v = corr_df.values[i, j]
            ax3.text(j, i, f"{v:.2f}", ha="center", va="center",
                     fontsize=7, color="black" if abs(v) < 0.55 else "white")
    cb = fig3.colorbar(im, ax=ax3, fraction=0.046, pad=0.04)
    cb.set_label("Spearman rho (interaction ranks)", fontsize=8)
    ax3.set_title("Cross-method ranking agreement\n(low off-diagonal -> consensus needed)", fontsize=9)
    fig3.tight_layout()
    save_fig(fig3, str(DIR_ASSETS / "fig3_cross_method_heatmap")); plt.close(fig3)

    # ---- Fig 4: 共识 vs 单方法 top-K 重叠 (★诚实对照 2, lollipop) -----------
    print("[fig 4] consensus vs single-method overlap (honest baseline 2, lollipop) ...")
    fig4, ax4 = plt.subplots(figsize=(NATURE_W1 + 1.8, NATURE_W1 - 0.1))
    cols = pal(len(ov), name="npg")
    yy = np.arange(len(ov))
    ax4.hlines(yy, 0, ov.values, color="#bbbbbb", lw=2, zorder=1)
    ax4.scatter(ov.values, yy, s=95, c=cols, zorder=3, edgecolor="white", linewidth=0.8)
    for y, v in zip(yy, ov.values):
        ax4.text(v + 0.02, y, f"{v:.2f}", va="center", fontsize=8)
    ax4.set_yticks(yy); ax4.set_yticklabels(ov.index, fontsize=9)
    ax4.set_xlim(0, max(0.1, ov.max() * 1.25))
    ax4.set_xlabel(f"top-{K} Jaccard overlap with consensus", fontsize=9)
    ax4.set_title("If you trusted ONE method, how far from consensus?\n"
                  "(lower = single method diverges from consensus)", fontsize=9)
    ax4.spines[["top", "right"]].set_visible(False)
    fig4.tight_layout()
    save_fig(fig4, str(DIR_ASSETS / "fig4_consensus_vs_single_lollipop")); plt.close(fig4)

    # ---- Fig 5: source × target tileplot (跨细胞对的幅度/特异 tile) ---------
    #   fill=magnitude_rank (颜色), label=specificity_rank (格内数字), 顶刊偏好的
    #   信息密集 tile, 一眼看出"哪对细胞间、哪条 L-R 既强又特异"。
    print("[fig 5] source x target tileplot (magnitude fill) ...")
    try:
        # fill=magnitude_rank(填色越亮越强);label 同填 magnitude_rank 但用空白格式函数
        # 隐藏密集数字,保持 tile 干净(默认把数值印满每格,过载难读)。
        fig5 = li.pl.tileplot(
            adata=adata, uns_key="liana_res",
            fill="magnitude_rank", label="specificity_rank",
            label_fun=lambda x: "",                    # 不在格内打数字,fill 已表达强度
            top_n=args.top_n, orderby="magnitude_rank", orderby_ascending=True,
            cmap=CMAP_CONT, figure_size=(8.5, 6.5), return_fig=True,
        )
        from plotnine import labs as _labs
        fig5 = fig5 + _labs(title="Consensus L-R across source/target (tileplot)")
        fig5.save(str(DIR_ASSETS / "fig5_source_target_tileplot.png"), dpi=300, verbose=False)
        fig5.save(str(DIR_ASSETS / "fig5_source_target_tileplot.pdf"), verbose=False)
    except Exception as e:
        print(f"[warn] tileplot failed ({e}); skipping fig5")

    # ---- 版本记录 (依赖快照, 铁律 6) ---------------------------------------
    import importlib.metadata as ilm
    with open(res_dir / "versions.txt", "w", encoding="utf-8") as f:
        for pkg in ("liana", "scanpy", "anndata", "numpy", "pandas",
                    "scipy", "matplotlib", "plotnine"):
            try:
                f.write(f"{pkg}=={ilm.version(pkg)}\n")
            except Exception:
                f.write(f"{pkg}==NA\n")
    print("[done] results/ + assets/ written.")


def _manual_network(cons, levels):
    """source->target 通讯网络 (圆形布局有向图)。
    节点大小 = 该 source 外向显著 L-R 总数;边宽 = 该 source->target 显著 L-R 计数。
    比 li.pl.circle_plot 更可控 (标签不裁切、纯 matplotlib 矢量导出稳)。"""
    counts = cons.groupby(["source", "target"]).size().reset_index(name="n")
    levels = list(levels)
    ang = {ct: 2 * np.pi * i / len(levels) for i, ct in enumerate(levels)}
    pos = {ct: (np.cos(a), np.sin(a)) for ct, a in ang.items()}
    cols = pal(len(levels), name="npg")
    cmap = {ct: cols[i] for i, ct in enumerate(levels)}

    fig, ax = plt.subplots(figsize=(6.4, 6.0))
    nmax = max(counts["n"].max(), 1)
    for _, row in counts.iterrows():
        s, t, n = row["source"], row["target"], row["n"]
        if s not in pos or t not in pos or s == t:
            continue
        x0, y0 = pos[s]; x1, y1 = pos[t]
        ax.annotate("", xy=(x1 * 0.8, y1 * 0.8), xytext=(x0 * 0.8, y0 * 0.8),
                    arrowprops=dict(arrowstyle="-|>", color=cmap[s],
                                    lw=0.5 + 4 * n / nmax, alpha=0.6,
                                    connectionstyle="arc3,rad=0.15"))
    deg = counts.groupby("source")["n"].sum()
    for ct, (x, y) in pos.items():
        ax.scatter(x, y, s=200 + 600 * deg.get(ct, 0) / max(deg.max(), 1),
                   c=[cmap[ct]], edgecolor="white", linewidth=1.2, zorder=5)
        ha = "left" if x > 0.25 else ("right" if x < -0.25 else "center")
        va = "bottom" if y > 0.25 else ("top" if y < -0.25 else "center")
        ax.text(x * 1.28, y * 1.28, ct, ha=ha, va=va, fontsize=10)
    ax.set_xlim(-2.0, 2.0); ax.set_ylim(-1.8, 1.8)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title("Cell-cell communication network\n(edge width = # significant L-R, source-coloured)",
                 fontsize=11)
    fig.tight_layout()
    save_fig(fig, str(DIR_ASSETS / "fig2_communication_network")); plt.close(fig)


if __name__ == "__main__":
    main()
