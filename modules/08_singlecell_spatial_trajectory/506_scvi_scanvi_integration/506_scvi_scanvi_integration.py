# -*- coding: utf-8 -*-
# =============================================================================
# 506 · scVI / scANVI 单细胞整合与标签转移 (deep generative integration)
# -----------------------------------------------------------------------------
# 用 scvi-tools 的深度生成模型做批次整合(scVI)+ 半监督标签转移(scANVI),
# 并与「未校正 PCA」基线对照,同时报告【批次混合 ↑】与【生物保留(细胞类型纯度)】
# 两个指标 —— 诚实地展示整合方法的取舍(过度混合会抹掉生物信号)。
#
# ★诚实基线(两把屠刀之一,Genome Biol 2025 10.1186/s13059-025-03574-x):
#   单细胞嵌入/整合务必与 PCA / Harmony 等简单基线对照,且同时看「批次混合」
#   与「生物保留」两面,不可只报一个好看的指标。本模块内置 PCA 基线 + 双指标。
#
# Turnkey: python 506_scvi_scanvi_integration.py   (默认合成数据→results/+assets/)
#          换数据: --input your.h5ad --batch_key batch --label_key celltype
# 复用 _framework/pubstyle.py;图中文字英文,注释中文。CPU 即可跑通小示例。
# =============================================================================
from __future__ import annotations
import argparse
import sys
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")

# ---- 定位脚本目录 + 载入顶刊绘图框架(向上搜 _framework) --------------------
HERE = Path(__file__).resolve().parent
for up in [HERE, *HERE.parents]:
    if (up / "_framework" / "pubstyle.py").exists():
        sys.path.insert(0, str(up / "_framework"))
        break
try:
    from pubstyle import set_pub_style, save_fig, pal, panel_labels, NATURE_W2
except Exception:  # 框架缺失时最小降级,不影响分析
    def set_pub_style(*a, **k): pass
    def save_fig(fig, f, dpi=300):
        from pathlib import Path as _P
        _P(f).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(str(f) + ".pdf"); fig.savefig(str(f) + ".png", dpi=dpi)
    def pal(n=None, name="npg"): import matplotlib; return list(matplotlib.cm.tab10.colors)
    def panel_labels(*a, **k): pass
    NATURE_W2 = 7.0

SEED = 42
np.random.seed(SEED)


# ---- 合成数据:2 批次 × 3 细胞类型,真实 scRNA 样计数(NB 过离散 + dropout)----
def make_synthetic(n_genes=300, n_per=150):
    """3 细胞类型(各有 marker 程序)× 2 批次。生成更贴近真实 scRNA 的计数:
    负二项(gamma-poisson 过离散)+ dropout 置零,批次效应=适中基因偏移 + 文库缩放。
    设计意图:未校正(PCA)按【批次】分离;scVI 条件于 batch 后应按【细胞类型】聚拢、
    批次混合 —— 既能体现整合收益,又不至于强到无法校正(诚实展示方法的适用边界)。"""
    types = ["Tcell", "Bcell", "Myeloid"]
    rng = np.random.default_rng(SEED)
    # 细胞类型 marker 程序(各 40 个基因升高);底噪较低
    base = rng.normal(0.1, 0.25, size=(len(types), n_genes))
    for i in range(len(types)):
        mk = slice(i * 40, i * 40 + 40)
        base[i, mk] += rng.normal(1.6, 0.3, size=40)
    # 批次 2:适中技术偏移(100 基因 +0.6)+ 文库缩放(整体 ×~1.25)
    batch_off = np.zeros((2, n_genes))
    bsel = rng.choice(n_genes, size=100, replace=False)
    batch_off[1, bsel] += 0.6
    lib_scale = np.array([1.0, 1.25])
    disp = 2.0   # NB 离散度(越小越过离散)

    X, ct, bt = [], [], []
    for b in range(2):
        for ti, t in enumerate(types):
            mu = np.exp(base[ti] + batch_off[b] + rng.normal(0, 0.2, size=(n_per, n_genes))) * lib_scale[b]
            # gamma-poisson = 负二项(过离散),再做 dropout
            lam = rng.gamma(shape=disp, scale=mu / disp)
            counts = rng.poisson(lam)
            drop = rng.random(counts.shape) < 0.2          # 20% dropout 置零
            counts[drop] = 0
            X.append(counts)
            ct += [t] * n_per
            bt += [f"batch{b+1}"] * n_per
    X = np.vstack(X).astype(np.float32)
    genes = [f"G{i:03d}" for i in range(n_genes)]
    return X, np.array(ct), np.array(bt), genes


# ---- 简单 kNN 指标:批次混合熵 + 细胞类型 kNN 纯度 --------------------------
def knn_metrics(emb, batch, celltype, k=30):
    """批次混合熵(越高越混合,理想整合↑)与细胞类型纯度(越高生物越保留)。"""
    from sklearn.neighbors import NearestNeighbors
    nn = NearestNeighbors(n_neighbors=k + 1).fit(emb)
    idx = nn.kneighbors(emb, return_distance=False)[:, 1:]
    # 批次混合熵(对每个细胞,邻居批次分布的香农熵,按批次数归一)
    ub = np.unique(batch); ent = []
    for i in range(emb.shape[0]):
        p = np.array([(batch[idx[i]] == b).mean() for b in ub]); p = p[p > 0]
        ent.append(-(p * np.log(p)).sum() / np.log(len(ub)))
    # 细胞类型纯度(邻居中与自身同型的比例)
    pur = np.mean([(celltype[idx[i]] == celltype[i]).mean() for i in range(emb.shape[0])])
    return float(np.mean(ent)), float(pur)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default=None, help="h5ad;留空用合成数据")
    ap.add_argument("--batch_key", default="batch")
    ap.add_argument("--label_key", default="celltype")
    ap.add_argument("--epochs", type=int, default=200)
    ap.add_argument("--outdir", default=str(HERE / "results"))
    args = ap.parse_args()

    set_pub_style(base_size=11)
    DRES = Path(args.outdir); DAST = HERE / "assets"; DDAT = HERE / "example_data"
    for d in (DRES, DAST, DDAT):
        d.mkdir(parents=True, exist_ok=True)

    import scanpy as sc
    import anndata as ad

    # ---- 1. 数据 ------------------------------------------------------------
    if args.input:
        adata = sc.read_h5ad(args.input)
        adata.layers["counts"] = adata.X.copy()
    else:
        fcsv = DDAT / "synthetic_counts.npz"
        X, ct, bt, genes = make_synthetic()
        np.savez(fcsv, X=X, ct=ct, bt=bt, genes=genes)
        adata = ad.AnnData(X.astype(np.float32))
        adata.var_names = genes
        adata.obs[args.label_key] = ct
        adata.obs[args.batch_key] = bt
        adata.layers["counts"] = adata.X.copy()
        print(f"[gen] synthetic: {adata.n_obs} cells x {adata.n_vars} genes, "
              f"2 batches x 3 types (for demo only)")

    # ---- 2. 基线:未校正 PCA(log-normalize → PCA → UMAP) ------------------
    base = adata.copy()
    sc.pp.normalize_total(base, target_sum=1e4); sc.pp.log1p(base)
    sc.pp.scale(base, max_value=10)
    sc.tl.pca(base, n_comps=20, random_state=SEED)
    sc.pp.neighbors(base, use_rep="X_pca", random_state=SEED)
    sc.tl.umap(base, random_state=SEED)
    adata.obsm["X_umap_pca"] = base.obsm["X_umap"]
    ent_pca, pur_pca = knn_metrics(base.obsm["X_pca"],
                                   adata.obs[args.batch_key].values,
                                   adata.obs[args.label_key].values)

    # ---- 3. scVI 整合 -------------------------------------------------------
    import scvi
    scvi.settings.seed = SEED
    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key=args.batch_key)
    model = scvi.model.SCVI(adata, n_latent=10, n_layers=1, gene_likelihood="nb")
    model.train(max_epochs=args.epochs, accelerator="cpu", enable_progress_bar=False)
    adata.obsm["X_scVI"] = model.get_latent_representation()
    sci = adata.copy()
    sc.pp.neighbors(sci, use_rep="X_scVI", random_state=SEED)
    sc.tl.umap(sci, random_state=SEED)
    adata.obsm["X_umap_scvi"] = sci.obsm["X_umap"]
    ent_scvi, pur_scvi = knn_metrics(adata.obsm["X_scVI"],
                                     adata.obs[args.batch_key].values,
                                     adata.obs[args.label_key].values)

    # ---- 4. scANVI 半监督标签转移(隐藏一部分标签当作未知,预测后对真值评估)----
    rng = np.random.default_rng(SEED)
    partial = adata.obs[args.label_key].astype(str).copy().values
    held = rng.choice(adata.n_obs, size=int(0.4 * adata.n_obs), replace=False)
    partial[held] = "Unknown"
    adata.obs["labels_partial"] = partial
    scanvi = scvi.model.SCANVI.from_scvi_model(
        model, adata=adata, labels_key="labels_partial", unlabeled_category="Unknown")
    scanvi.train(max_epochs=max(20, args.epochs // 3), accelerator="cpu",
                 enable_progress_bar=False)
    adata.obs["scanvi_pred"] = scanvi.predict()
    truth = adata.obs[args.label_key].values
    acc = float((adata.obs["scanvi_pred"].values[held] == truth[held]).mean())

    # ---- 5. 落盘关键统计 ----------------------------------------------------
    import pandas as pd
    metrics = pd.DataFrame({
        "method": ["PCA (uncorrected)", "scVI (integrated)"],
        "batch_mixing_entropy": [ent_pca, ent_scvi],
        "celltype_purity": [pur_pca, pur_scvi],
    })
    metrics.to_csv(DRES / "integration_metrics.csv", index=False)
    pd.DataFrame({"cell": adata.obs_names, "true": truth,
                  "scanvi_pred": adata.obs["scanvi_pred"].values,
                  "held_out": np.isin(np.arange(adata.n_obs), held)}
                 ).to_csv(DRES / "scanvi_predictions.csv", index=False)
    print(f"[metric] batch-mixing entropy  PCA={ent_pca:.3f} -> scVI={ent_scvi:.3f}  (higher=better mixed)")
    print(f"[metric] cell-type purity       PCA={pur_pca:.3f} -> scVI={pur_scvi:.3f}  (higher=biology kept)")
    print(f"[metric] scANVI label-transfer accuracy on 40% held-out = {acc:.3f}")

    # ---- 6. 出图(Nature 主题;矢量 PDF + 300dpi PNG)-----------------------
    import matplotlib.pyplot as plt
    types = sorted(np.unique(truth)); batches = sorted(np.unique(adata.obs[args.batch_key]))
    cmap_t = {t: pal(len(types), "okabe_ito")[i] for i, t in enumerate(types)}
    cmap_b = {b: pal(len(batches), "npg")[i] for i, b in enumerate(batches)}

    # Fig 1: UMAP 着色批次  PCA vs scVI(整合效果一眼可见)
    fig, ax = plt.subplots(1, 2, figsize=(NATURE_W2, 3.3))
    for a, (key, ttl) in zip(ax, [("X_umap_pca", "PCA (uncorrected)"),
                                  ("X_umap_scvi", "scVI (integrated)")]):
        U = adata.obsm[key]
        for b in batches:
            m = adata.obs[args.batch_key].values == b
            a.scatter(U[m, 0], U[m, 1], s=8, c=cmap_b[b], label=b, alpha=0.75, linewidths=0)
        a.set_title(ttl); a.set_xlabel("UMAP1"); a.set_ylabel("UMAP2")
        a.set_xticks([]); a.set_yticks([])
    ax[1].legend(title="Batch", markerscale=1.6, loc="upper right")
    panel_labels(ax)
    fig.tight_layout(); save_fig(fig, DAST / "umap_batch"); plt.close(fig)

    # Fig 2: scVI UMAP 着色细胞类型(生物结构保留)
    fig, a = plt.subplots(figsize=(4.2, 3.6))
    U = adata.obsm["X_umap_scvi"]
    for t in types:
        m = truth == t
        a.scatter(U[m, 0], U[m, 1], s=9, c=cmap_t[t], label=t, alpha=0.8, linewidths=0)
    a.set_title("scVI latent — cell types"); a.set_xlabel("UMAP1"); a.set_ylabel("UMAP2")
    a.set_xticks([]); a.set_yticks([]); a.legend(title="Cell type", markerscale=1.5)
    fig.tight_layout(); save_fig(fig, DAST / "umap_celltype"); plt.close(fig)

    # Fig 3: 整合 vs 生物保留 散点(scIB 风格)—— 理想在右上角,箭头示 PCA→scVI 的改善。
    #        顶刊偏好这种二维权衡散点而非条形图(信息更密、直接读出"既混合又保真")。
    fig, a = plt.subplots(figsize=(4.3, 3.9))
    mx = metrics["batch_mixing_entropy"].values; my = metrics["celltype_purity"].values
    cols = pal(2, "okabe_ito")
    a.annotate("", xy=(mx[1], my[1]), xytext=(mx[0], my[0]),
               arrowprops=dict(arrowstyle="-|>", color="grey", lw=1.6, alpha=0.7))
    for i, m in enumerate(metrics["method"]):
        a.scatter(mx[i], my[i], s=170, c=cols[i], edgecolor="black", linewidth=0.9, zorder=3)
        a.annotate(m, (mx[i], my[i]), textcoords="offset points",
                   xytext=(9, -4 if i == 0 else 6), fontsize=9)
    a.set_xlabel("Batch-mixing entropy  (→ better integrated)")
    a.set_ylabel("Cell-type purity  (↑ biology kept)")
    a.set_title("Integration vs biology conservation")
    a.set_xlim(min(mx) - 0.12, 1.03); a.set_ylim(min(my) - 0.04, 1.02)
    a.text(1.02, 1.0, "ideal", ha="right", va="top", fontsize=8, color="grey", style="italic")
    fig.tight_layout(); save_fig(fig, DAST / "metrics_scatter"); plt.close(fig)

    # Fig 4: scANVI 标签转移混淆矩阵(held-out)
    from sklearn.metrics import confusion_matrix
    cm = confusion_matrix(truth[held], adata.obs["scanvi_pred"].values[held], labels=types)
    cmn = cm / cm.sum(1, keepdims=True).clip(min=1)
    fig, a = plt.subplots(figsize=(4.0, 3.4))
    im = a.imshow(cmn, cmap="viridis", vmin=0, vmax=1)
    a.set_xticks(range(len(types))); a.set_xticklabels(types, rotation=30, ha="right")
    a.set_yticks(range(len(types))); a.set_yticklabels(types)
    for i in range(len(types)):
        for j in range(len(types)):
            a.text(j, i, f"{cmn[i,j]:.2f}", ha="center", va="center",
                   color="white" if cmn[i, j] < 0.6 else "black", fontsize=9)
    a.set_xlabel("scANVI predicted"); a.set_ylabel("True")
    a.set_title(f"scANVI label transfer (acc={acc:.2f})")
    fig.colorbar(im, ax=a, fraction=0.046, label="Row-normalized")
    fig.tight_layout(); save_fig(fig, DAST / "scanvi_confusion"); plt.close(fig)

    print(f"[fig] assets/: umap_batch, umap_celltype, metrics_scatter, scanvi_confusion (.pdf+.png)")
    # 依赖快照(铁律6)
    try:
        import session_info
        session_info.show(write_req_file=False)
    except Exception:
        print(f"[env] scvi-tools={scvi.__version__}, scanpy={sc.__version__}, numpy={np.__version__}")


if __name__ == "__main__":
    main()
