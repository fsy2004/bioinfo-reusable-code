# -*- coding: utf-8 -*-
# =============================================================================
# 547 · ProLIF 蛋白-配体相互作用指纹 (Protein-Ligand Interaction Fingerprint)
# -----------------------------------------------------------------------------
# 做什么 : 把对接/MD 的结合 pose(单帧或多帧轨迹)解码成「残基 × 相互作用类型」
#          的相互作用指纹(IFP)。逐帧检测氢键/疏水/盐桥/π-堆叠/范德华等接触,
#          统计每个口袋残基的接触占据率(occupancy %),画 barcode、heatmap、
#          每残基接触频率 lollipop,定量回答「配体到底抓住了哪些关键残基」。
# ★诚实基线 : 不主观挑选"关键残基"。指纹直接逐帧客观检测,多帧报告每个
#          (残基, 相互作用) 的出现频率 = 占据率%(0~100)。占据率高 = 稳定接触,
#          低 = 偶发/噪声;由数据说话而非作者钦点。单帧时占据率即 0/100。
# Turnkey : python 547_prolif_interaction_fingerprint.py
#           (默认用 prolif 自带真实复合物 TOP/TRAJ:一个 FXa 类蛋白-配体 MD 轨迹)
# 换数据  : python 547_prolif_interaction_fingerprint.py --top my.pdb --traj my.xtc \
#                  --ligand "resname LIG" --protein "protein" --frames 50
#           (单帧 pose 也可:--top complex.pdb 不给 --traj,自动当 1 帧处理)
# 真实 API: prolif 2.2.0 · MDAnalysis 2.10 · rdkit 2026.03
#           fp = plf.Fingerprint(); fp.run(u.trajectory[sel], lig, prot)
#           df = fp.to_dataframe()  → 3 级列 MultiIndex (ligand, protein, interaction)
# =============================================================================
from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

# ---- 复用顶刊绘图框架 -------------------------------------------------------
import sys
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
        from pathlib import Path as _P; _P(f).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(str(f) + ".pdf"); fig.savefig(str(f) + ".png", dpi=dpi)
    def pal(n=None, name="npg"):
        import matplotlib; return list(matplotlib.cm.tab10.colors)
    def panel_labels(*a, **k): pass
    CMAP_CONT = "viridis"; CMAP_DIVERGE = "RdBu_r"; NATURE_W1 = 3.5; NATURE_W2 = 7.0

# ---- 路径(全部从脚本位置派生,绝不 hardcode / 绝不 chdir) -----------------
DIR_EX = HERE / "example_data"
DIR_RES = HERE / "results"
DIR_ASSETS = HERE / "assets"
for d in (DIR_EX, DIR_RES, DIR_ASSETS):
    d.mkdir(parents=True, exist_ok=True)

SEED = 42

# 相互作用类型 → 配色(固定映射,跨图一致)
INTERACTION_ORDER = ["Hydrophobic", "VdWContact", "HBDonor", "HBAcceptor",
                     "Cationic", "Anionic", "PiStacking", "CationPi", "PiCation",
                     "XBDonor", "XBAcceptor", "MetalDonor", "MetalAcceptor",
                     "EdgeToFace", "FaceToFace"]


def log(msg: str) -> None:
    print(f"[547] {msg}", flush=True)


# ----------------------------------------------------------------------------
# 1. 取数据:默认用 prolif 自带真实复合物;--top/--traj 覆盖
# ----------------------------------------------------------------------------
def load_universe(top, traj, lig_sel, prot_sel):
    import MDAnalysis as mda
    import prolif as plf

    if top is None:
        # prolif 官方打包的真实 MD 数据(蛋白-配体复合物,resname LIG)
        top = plf.datafiles.TOP
        traj = plf.datafiles.TRAJ
        lig_sel = "resname LIG"
        prot_sel = "protein"
        log(f"未指定 --top,使用 prolif 自带真实复合物 demo: {Path(top).name} + {Path(traj).name}")
    else:
        log(f"使用用户数据: top={top} traj={traj}")

    if traj is None:
        # 单帧 pose:仅拓扑文件,MDAnalysis 当 1 帧轨迹
        u = mda.Universe(top)
    else:
        u = mda.Universe(top, traj)

    lig = u.select_atoms(lig_sel)
    prot = u.select_atoms(prot_sel)
    if lig.n_atoms == 0:
        raise ValueError(f"配体选择 '{lig_sel}' 为空,请检查 --ligand")
    if prot.n_atoms == 0:
        raise ValueError(f"蛋白选择 '{prot_sel}' 为空,请检查 --protein")
    log(f"配体原子 {lig.n_atoms} · 蛋白原子 {prot.n_atoms} · 轨迹总帧 {len(u.trajectory)}")
    return u, lig, prot


# ----------------------------------------------------------------------------
# 2. 计算相互作用指纹
# ----------------------------------------------------------------------------
def compute_fingerprint(u, lig, prot, n_frames):
    import prolif as plf

    total = len(u.trajectory)
    n_use = min(n_frames, total)
    sel = u.trajectory[:n_use]
    log(f"在 {n_use} 帧上计算相互作用指纹(IFP)…")

    fp = plf.Fingerprint()           # 默认检测一组标准相互作用类型
    fp.run(sel, lig, prot, progress=False)
    df = fp.to_dataframe()           # 行=帧, 列=(ligand, protein, interaction) 的 3 级 MultiIndex, 值=bool
    log(f"指纹完成: {df.shape[0]} 帧 × {df.shape[1]} 个 (残基,相互作用) 列")
    return df, n_use


# ----------------------------------------------------------------------------
# 3. 诚实基线:逐帧客观占据率(不挑残基)
# ----------------------------------------------------------------------------
def occupancy_tables(df: pd.DataFrame):
    """从 IFP dataframe 派生占据率表。"""
    # 列 = (ligand, protein, interaction);把 ligand 层折叠(单配体)
    occ = df.mean(axis=0) * 100.0    # 每个 (残基,相互作用) 在所有帧里出现的频率%
    occ.index = occ.index.droplevel("ligand")
    occ_df = occ.reset_index()
    occ_df.columns = ["residue", "interaction", "occupancy_pct"]

    # 残基排序键:取残基号便于沿序列排
    def resnum(r):
        s = "".join(ch for ch in r.split(".")[0] if ch.isdigit())
        return int(s) if s else 0
    occ_df["resnum"] = occ_df["residue"].map(resnum)
    occ_df = occ_df.sort_values(["resnum", "interaction"]).reset_index(drop=True)

    # 每残基总接触频率(任意类型出现即算接触)= 诚实排序依据
    # 用 OR 折叠每残基的多种相互作用 → 该残基在多少比例帧里有接触
    # (转置后按 protein 层 groupby,避开 pandas 已弃用的 groupby(axis=1))
    res_any = df.T.groupby(level="protein").any().T.mean(axis=0) * 100.0
    per_res = res_any.reset_index()
    per_res.columns = ["residue", "contact_freq_pct"]
    per_res["resnum"] = per_res["residue"].map(resnum)
    per_res = per_res.sort_values("contact_freq_pct", ascending=False).reset_index(drop=True)
    return occ_df, per_res


# ----------------------------------------------------------------------------
# 4. 绘图(非条形)
# ----------------------------------------------------------------------------
def plot_barcode(df: pd.DataFrame, out_prefix: Path):
    """Barcode:帧(x) × (残基-相互作用)(y) 的存在/缺失栅格,直观看接触在轨迹里是否持续。"""
    mat = df.astype(int).T                       # 行=(lig,res,int), 列=帧
    res_labels = [f"{r}·{i}" for (_, r, i) in mat.index]
    # 只保留至少出现过一次的接触,并按占据率降序,信息密度最高
    occ = mat.mean(axis=1)
    keep = occ[occ > 0].sort_values(ascending=False).index
    mat = mat.loc[keep]
    res_labels = [f"{r} · {i}" for (_, r, i) in mat.index]
    n_row = len(mat)
    h = max(2.2, 0.18 * n_row)
    fig, ax = plt.subplots(figsize=(NATURE_W2, h))
    ax.imshow(mat.values, aspect="auto", cmap="Greys", interpolation="nearest",
              vmin=0, vmax=1)
    ax.set_yticks(np.arange(n_row))
    ax.set_yticklabels(res_labels, fontsize=5)
    ax.set_xlabel("Trajectory frame")
    ax.set_title("Interaction fingerprint barcode  (black = contact present)")
    ax.set_xticks(np.linspace(0, mat.shape[1] - 1, min(6, mat.shape[1])).astype(int))
    save_fig(fig, str(out_prefix))
    plt.close(fig)


def plot_heatmap(occ_df: pd.DataFrame, out_prefix: Path):
    """Heatmap:残基(y) × 相互作用类型(x),格值=占据率%。"""
    piv = occ_df.pivot_table(index="residue", columns="interaction",
                             values="occupancy_pct", aggfunc="max", fill_value=0)
    # 残基按序列号排
    def resnum(r):
        s = "".join(ch for ch in r.split(".")[0] if ch.isdigit())
        return int(s) if s else 0
    piv = piv.loc[sorted(piv.index, key=resnum)]
    # 相互作用列按固定顺序
    cols = [c for c in INTERACTION_ORDER if c in piv.columns]
    cols += [c for c in piv.columns if c not in cols]
    piv = piv[cols]

    fig, ax = plt.subplots(figsize=(max(NATURE_W1, 0.6 * len(cols) + 1.5),
                                    max(2.5, 0.22 * len(piv) + 1)))
    im = ax.imshow(piv.values, aspect="auto", cmap=CMAP_CONT, vmin=0, vmax=100)
    ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(len(piv))); ax.set_yticklabels(piv.index, fontsize=6)
    ax.set_title("Per-residue interaction occupancy (%)")
    # 在格内标数值(>0)
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            v = piv.values[i, j]
            if v > 0:
                ax.text(j, i, f"{v:.0f}", ha="center", va="center",
                        fontsize=5, color="white" if v > 55 else "black")
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("Occupancy (%)", fontsize=7)
    save_fig(fig, str(out_prefix))
    plt.close(fig)


def plot_lollipop(per_res: pd.DataFrame, out_prefix: Path):
    """Lollipop:每残基整体接触频率%(诚实基线排序,top 接触残基一目了然)。"""
    d = per_res.sort_values("contact_freq_pct", ascending=True)
    n = len(d)
    fig, ax = plt.subplots(figsize=(NATURE_W1 + 1.5, max(2.5, 0.26 * n + 0.8)))
    colors = pal(1, "npg")[0]
    y = np.arange(n)
    ax.hlines(y, 0, d["contact_freq_pct"].values, color="#B0B0B0", lw=1.4, zorder=1)
    ax.scatter(d["contact_freq_pct"].values, y, s=46, color=colors,
               edgecolor="black", lw=0.6, zorder=2)
    ax.set_yticks(y); ax.set_yticklabels(d["residue"].values, fontsize=7)
    ax.set_xlabel("Contact frequency across frames (%)")
    ax.set_xlim(0, max(105, d["contact_freq_pct"].max() * 1.1))
    ax.set_title("Per-residue contact occupancy (honest ranking)")
    for yi, v in zip(y, d["contact_freq_pct"].values):
        ax.text(v + 1.5, yi, f"{v:.0f}%", va="center", fontsize=6)
    ax.spines[["top", "right"]].set_visible(False)
    save_fig(fig, str(out_prefix))
    plt.close(fig)


# ----------------------------------------------------------------------------
# 5. 主流程
# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="547 ProLIF 蛋白-配体相互作用指纹 (turnkey)")
    ap.add_argument("--top", default=None, help="拓扑文件(.pdb/.gro/.psf 等);省略=用 prolif 自带真实 demo")
    ap.add_argument("--traj", default=None, help="轨迹文件(.xtc/.dcd 等);省略+给 top=单帧 pose")
    ap.add_argument("--ligand", default="resname LIG", help="配体 MDAnalysis 选择语句")
    ap.add_argument("--protein", default="protein", help="蛋白 MDAnalysis 选择语句")
    ap.add_argument("--frames", type=int, default=50, help="最多分析前 N 帧(demo 加速)")
    args = ap.parse_args()

    set_pub_style(base_size=11, palette="npg")
    np.random.seed(SEED)

    log("=== ProLIF interaction fingerprint ===")
    u, lig, prot = load_universe(args.top, args.traj, args.ligand, args.protein)
    df, n_use = compute_fingerprint(u, lig, prot, args.frames)

    # 把原始指纹存盘(便于复用)
    df_save = df.copy()
    df_save.columns = ["|".join(map(str, c)) for c in df_save.columns]
    df_save.to_csv(DIR_RES / "fingerprint_per_frame.csv", index_label="frame")

    log("派生诚实占据率表…")
    occ_df, per_res = occupancy_tables(df)
    occ_df.drop(columns="resnum").to_csv(DIR_RES / "occupancy_by_residue_interaction.csv", index=False)
    per_res.drop(columns="resnum").to_csv(DIR_RES / "contact_freq_per_residue.csv", index=False)
    log(f"检测到 {per_res['residue'].nunique()} 个接触残基; "
        f"top: {per_res.iloc[0]['residue']} ({per_res.iloc[0]['contact_freq_pct']:.0f}%)")

    log("绘图(barcode / heatmap / lollipop)…")
    plot_barcode(df, DIR_ASSETS / "547_barcode")
    plot_heatmap(occ_df, DIR_ASSETS / "547_heatmap_occupancy")
    plot_lollipop(per_res, DIR_ASSETS / "547_lollipop_contact_freq")

    # 版本写盘
    import prolif as plf, MDAnalysis as mda, rdkit, matplotlib as mpl
    (DIR_RES / "versions.txt").write_text(
        f"prolif {plf.__version__}\nMDAnalysis {mda.__version__}\n"
        f"rdkit {rdkit.__version__}\nnumpy {np.__version__}\n"
        f"pandas {pd.__version__}\nmatplotlib {mpl.__version__}\n"
        f"n_frames_analyzed {n_use}\nseed {SEED}\n", encoding="utf-8")

    log("完成。结果见 results/,图见 assets/")


if __name__ == "__main__":
    main()
