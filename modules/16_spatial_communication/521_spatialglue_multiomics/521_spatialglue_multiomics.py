# -*- coding: utf-8 -*-
# =============================================================================
# 521 · SpatialGlue 空间多组学整合 (spatial multi-omics domain identification)
# -----------------------------------------------------------------------------
# SpatialGlue 用带跨模态注意力的图神经网络整合空间多组学(RNA+蛋白/ATAC)识别空间域。
# 本 turnkey 提供:① 可本地跑的【诚实基线】(各模态单独 + 简单拼接整合,kmeans→空间域,
# 用 ARI 对真值打分,证明"整合 > 单模态");② 真 SpatialGlue 的调用代码(需 GPU/装包,
# 见 README 服务器说明)。合成互补模态:RNA 分 {域1}vs{2,3}、蛋白分 {域3}vs{1,2},
# 单模态都分不全 3 域,整合才行 → 体现多组学价值。
#
# ★两把屠刀:GNN 整合务必与简单基线(拼接+PCA)对照(本模块内置 ARI 对照),别只跑 GNN。
# Turnkey: python 521_spatialglue_multiomics.py   (合成空间多组学→results/+assets/)
# 复用 _framework/pubstyle.py;无条形图(空间散点 + lollipop)。
# =============================================================================
from __future__ import annotations
import sys, warnings
from pathlib import Path
import numpy as np
warnings.filterwarnings("ignore")

HERE = Path(__file__).resolve().parent
for up in [HERE, *HERE.parents]:
    if (up / "_framework" / "pubstyle.py").exists():
        sys.path.insert(0, str(up / "_framework")); break
from pubstyle import set_pub_style, save_fig, pal, panel_labels, NATURE_W2  # noqa: E402

SEED = 42
np.random.seed(SEED)


def make_spatial_multiomics(side=30):
    """side×side spots,3 个竖条空间域;RNA 与蛋白模态信息互补。"""
    xs, ys = np.meshgrid(np.arange(side), np.arange(side))
    coord = np.c_[xs.ravel(), ys.ravel()].astype(float)
    domain = (coord[:, 0] // (side / 3)).astype(int)          # 0/1/2 三竖条域
    n = coord.shape[0]
    rng = np.random.default_rng(SEED)
    # RNA:60 基因;前 20 基因在【域0】强高 → RNA 只能分 {0} vs {1,2}
    rna = rng.poisson(2, size=(n, 60)).astype(float)
    rna[domain == 0, :20] += rng.poisson(10, size=((domain == 0).sum(), 20))
    # 蛋白:20 ADT;前 8 在【域2】强高 → 蛋白只能分 {2} vs {0,1}
    adt = rng.normal(1, 0.5, size=(n, 20))
    adt[domain == 2, :8] += rng.normal(6, 0.5, size=((domain == 2).sum(), 8))
    return coord, domain, rna, adt


def cluster_rep(rep, k=3, npc=10):
    from sklearn.cluster import KMeans
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    z = StandardScaler().fit_transform(rep)
    pcs = PCA(n_components=min(npc, rep.shape[1] - 1), random_state=SEED).fit_transform(z)
    return KMeans(k, n_init=10, random_state=SEED).fit_predict(pcs), pcs


def main():
    set_pub_style(base_size=11)
    DRES = HERE / "results"; DAST = HERE / "assets"; DDAT = HERE / "example_data"
    for d in (DRES, DAST, DDAT):
        d.mkdir(parents=True, exist_ok=True)
    from sklearn.metrics import adjusted_rand_score
    from sklearn.preprocessing import StandardScaler

    coord, domain, rna, adt = make_spatial_multiomics()
    np.savez(DDAT / "spatial_multiomics.npz", coord=coord, domain=domain, rna=rna, adt=adt)
    print(f"[gen] synthetic spatial multi-omics: {coord.shape[0]} spots, 3 domains, "
          f"RNA(60) + ADT(20), complementary (demo only)")

    # ---- 诚实基线:单模态 vs 拼接整合 -------------------------------------
    rna_lab, rna_pc = cluster_rep(rna)
    adt_lab, adt_pc = cluster_rep(adt)
    # 整合:各模态按自身 PC1 标准差归一(信号轴≈1、噪声轴<1),直接拼接 kmeans。
    # 不再对 joint 重 z-score(那会把噪声维抬到与信号同方差,淹没互补信号)。
    from sklearn.cluster import KMeans
    rna_s = rna_pc[:, :5] / rna_pc[:, 0].std()
    adt_s = adt_pc[:, :5] / adt_pc[:, 0].std()
    joint = np.c_[rna_s, adt_s]
    int_lab = KMeans(3, n_init=10, random_state=SEED).fit_predict(joint)
    aris = {"RNA only": adjusted_rand_score(domain, rna_lab),
            "Protein only": adjusted_rand_score(domain, adt_lab),
            "Integrated (concat+PCA)": adjusted_rand_score(domain, int_lab)}
    print("[baseline] domain-recovery ARI:  " +
          "  ".join(f"{k}={v:.2f}" for k, v in aris.items()))

    # ---- 真 SpatialGlue(需装包/GPU;缺则跳过,仅基线)--------------------
    glue_ari = None
    try:
        import torch  # noqa: F401
        import SpatialGlue  # noqa: F401  (真包,本机通常未装)
        # 真实流程(伪代码,服务器启用):构建两模态 AnnData + 空间图 →
        #   from SpatialGlue.preprocess import construct_neighbor_graph
        #   from SpatialGlue.SpatialGlue_pyG import Train_SpatialGlue
        #   model = Train_SpatialGlue(data, datatype='SPOTS'); emb = model.train()
        #   后再 kmeans/mclust 得 spatial domains。
        print("[spatialglue] real package detected — see README for the full GNN call")
    except Exception:
        print("[spatialglue] real SpatialGlue/torch-geometric NOT installed -> baseline only "
              "(install on server, see README)")

    import pandas as pd
    pd.DataFrame({"method": list(aris.keys()), "ARI": list(aris.values())}
                 ).to_csv(DRES / "domain_recovery_ARI.csv", index=False)
    pd.DataFrame({"x": coord[:, 0], "y": coord[:, 1], "true_domain": domain,
                  "integrated_domain": int_lab}).to_csv(DRES / "spatial_domains.csv", index=False)

    # ---- 出图(空间散点 + ARI lollipop;无条形图)------------------------
    import matplotlib.pyplot as plt
    dcol = pal(3, "okabe_ito")

    # Fig1: 空间域图 — 真值 vs 整合预测
    fig, ax = plt.subplots(1, 2, figsize=(NATURE_W2, 3.4))
    for a, (lab, ttl) in zip(ax, [(domain, "True domains"),
                                  (int_lab, f"Integrated (ARI={aris['Integrated (concat+PCA)']:.2f})")]):
        for d in np.unique(lab):
            m = lab == d
            a.scatter(coord[m, 0], coord[m, 1], s=14, c=dcol[d % 3], label=f"D{d}", linewidths=0)
        a.set_title(ttl); a.set_xticks([]); a.set_yticks([]); a.set_aspect("equal")
    ax[1].legend(title="Domain", markerscale=1.4, loc="center left", bbox_to_anchor=(1, 0.5))
    panel_labels(ax)
    fig.tight_layout(); save_fig(fig, DAST / "spatial_domains"); plt.close(fig)

    # Fig2: ARI 对照 lollipop(整合 > 单模态)
    fig, a = plt.subplots(figsize=(5.0, 2.8))
    names = list(aris.keys()); vals = list(aris.values())
    order = np.argsort(vals); names = [names[i] for i in order]; vals = [vals[i] for i in order]
    cols = [pal(3, "npg")[0] if "Integrated" in n else "#9E9E9E" for n in names]
    a.hlines(range(len(names)), 0, vals, color=cols, linewidth=2.5, alpha=0.7)
    a.scatter(vals, range(len(names)), c=cols, s=90, zorder=3)
    a.set_yticks(range(len(names))); a.set_yticklabels(names)
    a.set_xlabel("Domain-recovery ARI (vs ground truth)"); a.set_xlim(0, 1)
    a.set_title("Integration beats single modality")
    fig.tight_layout(); save_fig(fig, DAST / "ari_lollipop"); plt.close(fig)

    print("[fig] assets/: spatial_domains, ari_lollipop (.pdf+.png)")


if __name__ == "__main__":
    main()
