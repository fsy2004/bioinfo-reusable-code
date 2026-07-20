# -*- coding: utf-8 -*-
# =============================================================================
# 544 · PASTE 多空间切片对齐 / 3D 堆叠 (optimal-transport slice alignment)
# -----------------------------------------------------------------------------
# 分类: 16_spatial_communication   语言: Python   工具: PASTE (paste-bio) + POT
# 用途: 用最优传输(FGW)把相邻空间转录组切片配准对齐 —— pairwise_align 求切片间
#       spot↔spot 概率耦合 pi, stack_slices_pairwise 由 pi 解出刚体变换(旋转+平移)
#       把两张切片叠进同一坐标系, 支持部分重叠与 3D 堆叠。
#
# ★诚实基线(必须内置, 不可只报好看图):
#   对照「naive 坐标直接叠加」(不做任何对齐, 把两切片原始 spatial 坐标摞在一起)。
#   合成数据里切片 B 是切片 A 旋转 +30°、平移后的副本; 真实最优配准应把 B 转回
#   与 A 重合。本模块同时报 naive 与 PASTE 两种叠加的【配准残差】(配对 spot 间的
#   平均距离), 让 OT 对齐相对 naive 叠加修正了多少切片间偏移可量化、可证伪。
#
# Turnkey: python 544_paste2_slice_alignment.py        (合成数据 → results/+assets/)
#          换数据: --sliceA a.h5ad --sliceB b.h5ad      (两切片须各有 .obsm['spatial'])
# 真实 API(已最小试跑确认, POT 0.9.4):
#   paste.pairwise_align(A,B) -> pi (nA×nB, 总和=1 的耦合矩阵)
#   paste.stack_slices_pairwise([A,B],[pi], output_params=True)
#       -> (aligned_slices, thetas(旋转弧度), translations) ; 写回 .obsm['spatial']
#   AnnData 需 .obsm['spatial'](N×2) 且两切片共享 var_names。
#   ★坑: paste-bio 1.4.0 与 POT>=0.9.5 不兼容(line_search 回调签名变 6 参), 须 POT 0.9.4
#         (NumPy2 兼容 + 旧 5 参签名); POT 0.9.3 是 NumPy1 编译会崩。
# 复用 _framework/pubstyle.py; 图中文字英文, 注释中文。小规模 CPU 即可跑通。
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
    from pubstyle import set_pub_style, save_fig, pal, NATURE_W2
except Exception:  # 框架缺失时最小降级, 不影响分析
    def set_pub_style(*a, **k): pass
    def save_fig(fig, f, dpi=300):
        from pathlib import Path as _P
        _P(f).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(str(f) + ".pdf"); fig.savefig(str(f) + ".png", dpi=dpi)
    def pal(n=None, name="npg"):
        import matplotlib
        c = list(matplotlib.cm.tab10.colors)
        return c[:n] if n else c
    NATURE_W2 = 7.0

SEED = 42
np.random.seed(SEED)


# ---- 合成数据: 同组织相邻 2 切片, 切片 B = 切片 A 旋转/平移 + 抖动后的副本 -----
def make_synthetic(n_spots=200, n_genes=60, rot_deg=30.0, shift=(6.0, -4.0),
                   jitter=0.4, drop_frac=0.12, seed=SEED):
    """生成两张 *同组织相邻* 空间切片的合成 AnnData (synthetic, for demo only)。

    设计意图(为诚实基线服务, 且不把问题做得过于理想):
      - 切片 A: 一片带 3 个空间域 (domain) 的 spot, 每域一套 marker 基因程序;
      - 切片 B: 同一组织的相邻切片 = 把 A 的坐标整体【旋转 rot_deg + 平移 shift】,
        再叠加【spot 级随机抖动 jitter】(相邻切片 spot 网格不会像素级重合), 并随机
        丢弃 drop_frac 比例的 spot(部分重叠), 表达另抽计数 + 独立噪声(技术差异)。
      因此「真实配准」= 把 B 转回与 A 大致重合, 但不可能 0 残差;naive 直接叠加会留下
      整片旋转/平移偏移。PASTE 的 OT 对齐应显著缩小残差但非归零 —— 既证明收益又不造假。
    返回 truth: dict(rot, shift, pairs) ; pairs 给出 B-spot -> 其源 A-spot 的真值索引。
    """
    import anndata as ad
    rng = np.random.default_rng(seed)

    # spot 坐标 (切片 A): 均匀铺面
    coords = rng.uniform(0, 20, size=(n_spots, 2))
    # 3 个空间域: 按 x 坐标分 3 段, 各域一套 marker 基因升高
    domain = np.digitize(coords[:, 0], np.quantile(coords[:, 0], [1/3, 2/3]))
    dom_names = np.array(["DomainL", "DomainM", "DomainR"])[domain]

    base = rng.normal(0.2, 0.3, size=(n_spots, n_genes))
    blk = n_genes // 3
    for d in range(3):
        m = domain == d
        mk = slice(d * blk, d * blk + blk)
        base[np.ix_(np.where(m)[0], np.arange(mk.start, mk.stop))] += \
            rng.normal(1.8, 0.3, size=(m.sum(), blk))
    # 计数: 负二项 (gamma-poisson 过离散), 模拟真实 ST 计数
    mu_A = np.exp(base) * 1.0
    XA = rng.poisson(rng.gamma(2.0, mu_A / 2.0)).astype("float32")

    # 切片 B: 随机保留 (1-drop_frac) 个 spot (部分重叠), 表达另抽计数 + 噪声
    keep = np.sort(rng.choice(n_spots, size=int(round((1 - drop_frac) * n_spots)),
                              replace=False))
    mu_B = np.exp(base[keep] + rng.normal(0, 0.2, size=(len(keep), n_genes))) * 1.1
    XB = rng.poisson(rng.gamma(2.0, mu_B / 2.0)).astype("float32")

    # 切片 B 坐标 = (A 子集) 旋转 + 平移 + spot 级抖动 (这就是要被 PASTE 校正的偏移)
    th = np.deg2rad(rot_deg)
    R = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
    coords_B = coords[keep] @ R.T + np.array(shift) + rng.normal(0, jitter, size=(len(keep), 2))

    genes = [f"g{i:03d}" for i in range(n_genes)]
    A = ad.AnnData(XA); A.var_names = genes
    B = ad.AnnData(XB); B.var_names = genes
    A.obsm["spatial"] = coords.copy()
    B.obsm["spatial"] = coords_B
    A.obs["domain"] = dom_names; B.obs["domain"] = dom_names[keep]
    A.obs_names = [f"A_{i}" for i in range(n_spots)]
    B.obs_names = [f"B_{i}" for i in range(len(keep))]
    # truth: B 的第 i 个 spot 真值对应 A 的第 keep[i] 个 spot
    return A, B, {"rot": rot_deg, "shift": shift, "pairs": keep}


def pairwise_residual(coordsA, coordsB_paired):
    """配对 spot (真值一一对应) 在同一坐标系下的平均欧氏残差。越小=配准越好。"""
    return float(np.mean(np.linalg.norm(coordsA - coordsB_paired, axis=1)))


def main():
    ap = argparse.ArgumentParser(description="PASTE optimal-transport slice alignment")
    ap.add_argument("--sliceA", default=None, help="切片A h5ad (须含 .obsm['spatial']); 留空用合成数据")
    ap.add_argument("--sliceB", default=None, help="切片B h5ad")
    ap.add_argument("--alpha", type=float, default=0.1,
                    help="PASTE FGW 权衡: 0=只看表达, 1=只看空间几何 (默认 0.1)")
    ap.add_argument("--outdir", default=str(HERE / "results"))
    args = ap.parse_args()

    set_pub_style(base_size=11)
    DRES = Path(args.outdir); DAST = HERE / "assets"; DDAT = HERE / "example_data"
    for d in (DRES, DAST, DDAT):
        d.mkdir(parents=True, exist_ok=True)

    import anndata as ad
    import paste
    import ot

    # ---- 1. 数据 ------------------------------------------------------------
    if args.sliceA and args.sliceB:
        A = ad.read_h5ad(args.sliceA); B = ad.read_h5ad(args.sliceB)
        paste.filter_for_common_genes([A, B])
        truth = None
        print(f"[data] real slices: A={A.shape}, B={B.shape} (common genes filtered)")
    else:
        A, B, truth = make_synthetic()
        A.write(DDAT / "sliceA_synth.h5ad"); B.write(DDAT / "sliceB_synth.h5ad")
        print(f"[gen] synthetic 2 slices: A={A.n_obs} / B={B.n_obs} spots x {A.n_vars} genes, "
              f"3 domains; B = A rotated {truth['rot']:.0f}deg + shift {truth['shift']} "
              f"+ jitter + partial overlap (for demo only)")

    assert "spatial" in A.obsm and "spatial" in B.obsm, "两切片都需 .obsm['spatial']"
    cA0 = np.asarray(A.obsm["spatial"], float)
    cB0 = np.asarray(B.obsm["spatial"], float)
    pairs = truth["pairs"] if truth is not None else None   # B-spot -> 源 A-spot 索引

    # ---- 2. 诚实基线: naive 直接叠加 (不对齐) ------------------------------
    #   把两切片原始坐标摞在一起。合成数据有真值配对, 残差=切片间整体偏移(原始坐标系)。
    res_naive = pairwise_residual(cB0, cA0[pairs]) if pairs is not None else float("nan")

    # ---- 3. PASTE pairwise_align: 求 spot↔spot 最优传输耦合 pi ------------
    print(f"[paste] pairwise_align (alpha={args.alpha}) ... POT {ot.__version__}")
    pi = paste.pairwise_align(A, B, alpha=args.alpha, verbose=False)
    # pi: nA×nB, 行/列边际=均匀分布, 总和=1。每个 B-spot 的"软指派"质心给出匹配关系。
    pi_rowsum = pi.sum(1); pi_colsum = pi.sum(0)

    # ---- 4. stack_slices_pairwise: 由 pi 解刚体变换, 把两切片叠进同一坐标系 -
    aligned, thetas, trans = paste.stack_slices_pairwise([A, B], [pi], output_params=True)
    cA1 = np.asarray(aligned[0].obsm["spatial"], float)
    cB1 = np.asarray(aligned[1].obsm["spatial"], float)
    theta_recovered = float(thetas[-1])
    inj = -truth["rot"] if truth is not None else float("nan")
    print(f"[paste] recovered rotation = {np.rad2deg(theta_recovered):+.1f} deg "
          f"(needs {inj:+.1f} deg to undo injected {truth['rot'] if truth else float('nan')}deg)")

    # 对齐后残差(真值配对): B 的每个 spot 与其源 A-spot 在对齐坐标系下的距离。
    res_ot_truth = pairwise_residual(cB1, cA1[pairs]) if pairs is not None else float("nan")
    # pi 的 argmax 匹配: 每个 B-spot 最可能对应的 A-spot (真实无真值场景的通用关系)
    match_BtoA = pi.argmax(0)                 # len = nB, 值域 A 索引

    # 映射质量: 每个 B-spot 指派的"集中度"(max 概率 / 列和), 越高=匹配越确定
    map_conf = (pi.max(0) / np.clip(pi_colsum, 1e-12, None))

    # ---- 5. 落盘关键统计 ----------------------------------------------------
    import pandas as pd
    summary = pd.DataFrame({
        "method": ["naive overlay (no align)", "PASTE OT align"],
        "registration_residual": [res_naive, res_ot_truth],
    })
    summary.to_csv(DRES / "alignment_summary.csv", index=False)
    pd.DataFrame({
        "B_spot": B.obs_names,
        "matched_A_spot": [A.obs_names[i] for i in match_BtoA],
        "mapping_confidence": map_conf,
    }).to_csv(DRES / "spot_mapping.csv", index=False)
    np.save(DRES / "pi_coupling.npy", pi)

    print(f"[metric] registration residual (paired spots): "
          f"naive={res_naive:.3f}  ->  PASTE={res_ot_truth:.3f}  (lower=better aligned)")
    if truth is not None and res_naive > 0:
        print(f"[metric] PASTE reduces slice offset by "
              f"{100*(1-res_ot_truth/res_naive):.1f}% vs naive overlay")
    print(f"[metric] mean mapping confidence = {map_conf.mean():.3f} "
          f"(argmax-prob over uniform marginal)")

    # ---- 6. 出图 (Nature 主题; 矢量 PDF + 300dpi PNG; 禁条形图) -------------
    import matplotlib.pyplot as plt
    cols = pal(2, "npg")            # 切片 A / B 两色
    cA_col, cB_col = cols[0], cols[1]
    doms = sorted(np.unique(A.obs["domain"])) if "domain" in A.obs else None

    # Fig 1: 对齐前 vs 对齐后空间叠加散点 (一眼看出 OT 把 B 转回与 A 重合)
    fig, ax = plt.subplots(1, 2, figsize=(NATURE_W2, 3.5))
    for a, (ca, cb, ttl) in zip(
            ax,
            [(cA0, cB0, "Before — naive overlay"),
             (cA1, cB1, "After — PASTE OT align")]):
        a.scatter(ca[:, 0], ca[:, 1], s=14, c=cA_col, alpha=0.8, linewidths=0, label="Slice A")
        a.scatter(cb[:, 0], cb[:, 1], s=14, c=cB_col, alpha=0.55,
                  marker="^", linewidths=0, label="Slice B")
        a.set_title(ttl); a.set_xlabel("x"); a.set_ylabel("y")
        a.set_aspect("equal"); a.set_xticks([]); a.set_yticks([])
    ax[0].legend(loc="upper right", markerscale=1.4, frameon=False)
    if truth is not None:
        ax[0].text(0.02, 0.02, f"residual={res_naive:.2f}", transform=ax[0].transAxes,
                   fontsize=9, va="bottom", style="italic", color="grey")
        ax[1].text(0.02, 0.02, f"residual={res_ot_truth:.2f}", transform=ax[1].transAxes,
                   fontsize=9, va="bottom", style="italic", color="grey")
    fig.tight_layout(); save_fig(fig, DAST / "overlay_before_after"); plt.close(fig)

    # Fig 2: 配准叠加图 — 对齐后按空间域 (domain) 同色, 连匹配线证实生物结构对上
    fig, a = plt.subplots(figsize=(4.6, 4.2))
    if doms is not None:
        dpal = {d: pal(len(doms), "okabe_ito")[i] for i, d in enumerate(doms)}
        for d in doms:
            mA = A.obs["domain"].values == d
            a.scatter(cA1[mA, 0], cA1[mA, 1], s=22, c=dpal[d], alpha=0.9,
                      linewidths=0, label=f"A·{d}")
            mB = B.obs["domain"].values == d
            a.scatter(cB1[mB, 0], cB1[mB, 1], s=22, c=dpal[d], alpha=0.45,
                      marker="^", linewidths=0)
    # 抽样画 30 条匹配连线 (每个 B-spot 的 argmax A-spot 配对), 展示跨切片映射
    rng = np.random.default_rng(SEED)
    samp = rng.choice(len(cB1), size=min(30, len(cB1)), replace=False)
    for j in samp:
        i = match_BtoA[j]
        a.plot([cB1[j, 0], cA1[i, 0]], [cB1[j, 1], cA1[i, 1]],
                c="grey", lw=0.5, alpha=0.5, zorder=0)
    a.set_title("Registered overlay (A circles / B triangles)")
    a.set_xlabel("aligned x"); a.set_ylabel("aligned y")
    a.set_aspect("equal"); a.set_xticks([]); a.set_yticks([])
    if doms is not None:
        a.legend(loc="upper left", fontsize=7, frameon=False, ncol=1)
    fig.tight_layout(); save_fig(fig, DAST / "registration_overlay"); plt.close(fig)

    # Fig 3: 映射质量 — pi 耦合矩阵热图 (块对角=匹配清晰) + 置信度小提琴
    fig, ax = plt.subplots(1, 2, figsize=(NATURE_W2, 3.4),
                           gridspec_kw={"width_ratios": [1.15, 0.85]})
    # 3a: pi 热图 (按 A-spot 的 x 坐标排序, 真匹配应聚成对角带)
    oa = np.argsort(cA0[:, 0]); ob = np.argsort(cB0[:, 0])
    im = ax[0].imshow(pi[np.ix_(oa, ob)], cmap="viridis", aspect="auto")
    ax[0].set_title("OT coupling $\\pi$ (spot$\\times$spot)")
    ax[0].set_xlabel("Slice B spots (sorted)"); ax[0].set_ylabel("Slice A spots (sorted)")
    fig.colorbar(im, ax=ax[0], fraction=0.046, label="transport mass")
    # 3b: 映射置信度小提琴 + 抖点 (替代条形图)
    parts = ax[1].violinplot([map_conf], showmeans=True, showextrema=False)
    for b in parts["bodies"]:
        b.set_facecolor(cA_col); b.set_alpha(0.5); b.set_edgecolor("black")
    parts["cmeans"].set_color("black")
    jx = 1 + (rng.random(len(map_conf)) - 0.5) * 0.12
    ax[1].scatter(jx, map_conf, s=6, c="black", alpha=0.35, linewidths=0)
    ax[1].set_xticks([1]); ax[1].set_xticklabels(["A→B"])
    ax[1].set_ylabel("Mapping confidence (argmax / marginal)")
    ax[1].set_title("Per-spot mapping quality")
    fig.tight_layout(); save_fig(fig, DAST / "mapping_quality"); plt.close(fig)

    # Fig 4: 诚实基线对照 — 残差哑铃图 (dumbbell, 非条形): naive → PASTE 的下降
    if truth is not None:
        fig, a = plt.subplots(figsize=(4.6, 2.6))
        y = 0
        a.plot([res_naive, res_ot_truth], [y, y], c="grey", lw=2.4, zorder=1)
        a.scatter([res_naive], [y], s=190, c=pal(2, "okabe_ito")[0],
                  edgecolor="black", zorder=3, label="naive overlay")
        a.scatter([res_ot_truth], [y], s=190, c=pal(2, "okabe_ito")[1],
                  edgecolor="black", zorder=3, label="PASTE OT align")
        a.annotate(f"{res_naive:.2f}", (res_naive, y), textcoords="offset points",
                   xytext=(0, 12), ha="center", fontsize=9)
        a.annotate(f"{res_ot_truth:.2f}", (res_ot_truth, y), textcoords="offset points",
                   xytext=(0, 12), ha="center", fontsize=9)
        drop = 100 * (1 - res_ot_truth / res_naive)
        a.annotate(f"-{drop:.0f}% residual", ((res_naive + res_ot_truth) / 2, y),
                   textcoords="offset points", xytext=(0, -22), ha="center",
                   fontsize=10, color="firebrick", fontweight="bold")
        a.set_yticks([]); a.set_ylim(-1.0, 0.7)
        a.set_xlabel("Registration residual  (paired-spot mean distance, lower=better)")
        a.set_title("Honest baseline: OT alignment corrects slice offset", pad=10)
        a.legend(loc="lower center", ncol=2, frameon=False, fontsize=8,
                 bbox_to_anchor=(0.5, -0.02))
        fig.tight_layout(); save_fig(fig, DAST / "baseline_residual_dumbbell"); plt.close(fig)

    figs = "overlay_before_after, registration_overlay, mapping_quality"
    if truth is not None:
        figs += ", baseline_residual_dumbbell"
    print(f"[fig] assets/: {figs} (.pdf+.png)")

    # ---- 依赖快照 / sessionInfo (铁律6: 锁定版本以复现) ---------------------
    with open(DRES / "versions.txt", "w", encoding="utf-8") as f:
        import scipy
        f.write(f"python={sys.version.split()[0]}\n")
        f.write(f"paste-bio (import paste); POT={ot.__version__}\n")
        f.write(f"anndata={ad.__version__}\nnumpy={np.__version__}\nscipy={scipy.__version__}\n")
    try:                                  # session_info 若装则写完整快照
        import session_info
        session_info.show(write_req_file=False)
    except Exception:
        pass
    print(f"[env] POT={ot.__version__}, anndata={ad.__version__}, numpy={np.__version__}")
    print(f"[done] results -> {DRES}")


if __name__ == "__main__":
    main()
