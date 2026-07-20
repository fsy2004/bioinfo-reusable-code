# -*- coding: utf-8 -*-
# ============================================================================
# 556 · PoseBusters 对接 pose 物理有效性面板 (PoseBusters validity panel)
# ----------------------------------------------------------------------------
# 做什么 : 用 PoseBusters 对一批对接/AI 生成的小分子 3D pose 做物理有效性体检
#          (键长 / 键角 / 内部位阻冲突 / 芳环平面性 / 双键平面性 / 内能 等 ~12 项),
#          逐 pose 给出每项 PASS/FAIL,并汇总 PB-valid 整体通过率。
# ★诚实基线: 不只报"好看的通过率"。本模块内置 好构象 vs 坏构象 对照——
#          故意把好分子的原子重叠(位阻冲突)、拉长键(异常键长)造成坏构象,
#          直接报每项检查的通过率%,并证明 PoseBusters 能把坏构象抓出来。
#          这呼应"DL pose 生成器(DiffDock 等)+ PoseBusters 守门"的领域铁律:
#          再高的打分,过不了物理体检就不能信。
# Turnkey 用法 : python 556_posebusters_validity_panel.py
#                (零改动即跑:首次自动用 RDKit 合成 好/坏 构象到 example_data/)
# 换数据用法   : python 556_posebusters_validity_panel.py --input my_poses.sdf
#                (--input 给一个多构象 .sdf;每个构象当作一个 pose 体检)
# ----------------------------------------------------------------------------
# 真实 API (posebusters 0.6.5, 已实测):
#   from posebusters import PoseBusters
#   bust = PoseBusters(config="mol")              # "mol" 单分子体检,无需蛋白/参考
#   df = bust.bust([mol_pred], mol_true, mol_cond, full_report=False)
#       → 返回 DataFrame,每行一个 pose,各列为布尔检查项(True=PASS)。
# 绘图 : 复用 _framework/pubstyle.py;图中文字英文;不用平凡条形图。
# ============================================================================
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---- 复用顶刊绘图框架(定位 _framework 并 import pubstyle) -------------------
import sys

HERE = Path(__file__).resolve().parent
for up in [HERE, *HERE.parents]:
    if (up / "_framework" / "pubstyle.py").exists():
        sys.path.insert(0, str(up / "_framework"))
        break
try:
    from pubstyle import (
        set_pub_style, save_fig, pal, panel_labels,
        CMAP_CONT, CMAP_DIVERGE, NATURE_W1, NATURE_W2,
    )
except Exception:
    def set_pub_style(*a, **k):
        pass

    def save_fig(fig, f, dpi=300):
        from pathlib import Path as _P
        _P(f).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(str(f) + ".pdf")
        fig.savefig(str(f) + ".png", dpi=dpi)

    def pal(n=None, name="npg"):
        import matplotlib
        return list(matplotlib.cm.tab10.colors)

    def panel_labels(*a, **k):
        pass

    CMAP_CONT = "viridis"
    CMAP_DIVERGE = "RdBu_r"
    NATURE_W1 = 3.5
    NATURE_W2 = 7.0

# ---- 路径(全部从脚本位置派生,绝不 hardcode 绝对路径 / 绝不 chdir) ----------
DIR_EXAMPLE = HERE / "example_data"
DIR_RESULTS = HERE / "results"
DIR_ASSETS = HERE / "assets"

SEED = 42

# 一批演示分子(SMILES, 名称):覆盖常见药物骨架
DEMO_MOLS = [
    ("aspirin",      "CC(=O)Oc1ccccc1C(=O)O"),
    ("ibuprofen",    "CC(C)Cc1ccc(cc1)C(C)C(=O)O"),
    ("paracetamol",  "CC(=O)Nc1ccc(O)cc1"),
    ("caffeine",     "Cn1cnc2c1c(=O)n(C)c(=O)n2C"),
    ("naproxen",     "COc1ccc2cc(ccc2c1)C(C)C(=O)O"),
    ("benzamide",    "NC(=O)c1ccccc1"),
]


# ============================================================================
# 1) 合成示例数据:每个分子生成一个"好"构象 + 一个"坏"构象
# ============================================================================
def build_demo_sdf(out_sdf: Path) -> None:
    """用 RDKit 嵌入 3D 构象;偶数 pose = 好构象,奇数 pose = 故意破坏的坏构象。

    坏构象的破坏方式(物理上不合理,应被 PoseBusters 抓出):
      - 把某两个相邻原子坐标设为重合  → internal_steric_clash / 键长异常
      - 把另一个原子沿 x 平移 5 Å      → bond_lengths / bond_angles 异常
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Geometry import Point3D

    np.random.seed(SEED)               # 全局种子(可复现)
    rng = np.random.default_rng(SEED)  # 本地 RNG,显式传种子
    writer = Chem.SDWriter(str(out_sdf))
    n_written = 0
    for name, smi in DEMO_MOLS:
        base = Chem.MolFromSmiles(smi)
        if base is None:
            continue
        # ---- 好构象 ----
        good = Chem.AddHs(base)
        # 每个分子一个稳定但可复现的随机种子(显式传入 → 可复现)
        seed_g = int(rng.integers(0, 1_000_000))
        if AllChem.EmbedMolecule(good, randomSeed=seed_g) != 0:
            continue
        AllChem.MMFFOptimizeMolecule(good)
        good.SetProp("_Name", f"{name}__good")
        good.SetProp("pose_id", f"{name}__good")
        good.SetProp("expected", "valid")
        writer.write(good)
        n_written += 1

        # ---- 坏构象(从好构象复制后破坏) ----
        bad = Chem.Mol(good)
        conf = bad.GetConformer()
        heavy = [a.GetIdx() for a in bad.GetAtoms() if a.GetAtomicNum() > 1]
        if len(heavy) >= 3:
            a0, a1, a2 = heavy[0], heavy[1], heavy[2]
            p0 = conf.GetAtomPosition(a0)
            # 原子重叠 → 位阻冲突 + 键长塌缩
            conf.SetAtomPosition(a1, Point3D(p0.x, p0.y, p0.z))
            # 远离 → 异常长键
            p2 = conf.GetAtomPosition(a2)
            conf.SetAtomPosition(a2, Point3D(p2.x + 5.0, p2.y, p2.z))
        bad.SetProp("_Name", f"{name}__bad")
        bad.SetProp("pose_id", f"{name}__bad")
        bad.SetProp("expected", "invalid")
        writer.write(bad)
        n_written += 1
    writer.close()
    print(f"[data] 合成 {n_written} 个 pose(好/坏成对)→ {out_sdf.name}")


# ============================================================================
# 2) 读 pose 并跑 PoseBusters 体检
# ============================================================================
def load_poses(sdf_path: Path):
    """读 SDF;返回 (mols, pose_ids, expected_labels)。"""
    from rdkit import Chem

    mols, ids, exp = [], [], []
    supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False, sanitize=True)
    for i, m in enumerate(supplier):
        if m is None:
            print(f"[warn] SDF 第 {i} 个分子无法读取,跳过")
            continue
        pid = m.GetProp("pose_id") if m.HasProp("pose_id") else (
            m.GetProp("_Name") if m.HasProp("_Name") else f"pose_{i}")
        ex = m.GetProp("expected") if m.HasProp("expected") else "unknown"
        mols.append(m)
        ids.append(pid)
        exp.append(ex)
    return mols, ids, exp


def run_posebusters(mols, ids) -> pd.DataFrame:
    """对每个 pose 跑 PoseBusters(config='mol');返回 pose × check 的布尔表。"""
    from posebusters import PoseBusters

    # max_workers=0 → 不并行(Windows 下稳定、可复现)
    bust = PoseBusters(config="mol", max_workers=0)
    rows = []
    for pid, m in zip(ids, mols):
        df = bust.bust([m], None, None, full_report=False)
        rec = df.iloc[0].to_dict()
        rec = {k: rec[k] for k in df.columns}  # 仅保留检查列
        rec["pose_id"] = pid
        rows.append(rec)
        n_pass = sum(bool(v) for k, v in rec.items()
                     if k != "pose_id" and v is True)
        print(f"[busters] {pid:<22s}  {n_pass}/{len(df.columns)} checks PASS")
    out = pd.DataFrame(rows).set_index("pose_id")
    # 统一成数值布尔(NaN/None → False:无法评估视为未通过,守门从严)。
    # PoseBusters 在严重坏构象上 internal_energy 可能返回 NaN,显式按 False 处理。
    out = out.apply(lambda col: col.map(lambda v: bool(v) if pd.notna(v) else False))
    out = out.astype(float)
    return out


# ============================================================================
# 3) 绘图(全部非条形:tick-heatmap / lollipop / dumbbell)
# ============================================================================
def _short(name: str) -> str:
    """检查项名缩短为图内英文标签。"""
    return name.replace("_", " ")


def plot_tick_heatmap(checks: pd.DataFrame, out: Path):
    """图1: pose × check 通过/失败 tick heatmap(绿=PASS,红=FAIL)。"""
    M = checks.astype(float).values  # 1=PASS, 0=FAIL
    n_pose, n_chk = M.shape
    from matplotlib.colors import ListedColormap

    cmap = ListedColormap(["#D6604D", "#4DAF67"])  # red=fail, green=pass
    fig, ax = plt.subplots(figsize=(NATURE_W2, 0.42 * n_pose + 1.6))
    ax.imshow(M, cmap=cmap, vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(n_chk))
    ax.set_xticklabels([_short(c) for c in checks.columns],
                       rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(n_pose))
    ax.set_yticklabels(list(checks.index), fontsize=7)
    # 网格线
    ax.set_xticks(np.arange(-0.5, n_chk, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_pose, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.2)
    ax.tick_params(which="minor", length=0)
    # 在失败格画 ✗
    for i in range(n_pose):
        for j in range(n_chk):
            if M[i, j] < 0.5:
                ax.text(j, i, "x", ha="center", va="center",
                        color="white", fontsize=8, fontweight="bold")
    ax.set_title("PoseBusters per-pose physical validity (green = PASS, red = FAIL)")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color="#4DAF67", label="PASS"),
                       Patch(color="#D6604D", label="FAIL")],
              loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=8)
    save_fig(fig, str(out))
    plt.close(fig)


def plot_passrate_lollipop(checks: pd.DataFrame, out: Path):
    """图2: 每项检查的通过率 lollipop(按通过率排序;100% 一条参考线)。"""
    rate = checks.mean(axis=0).sort_values() * 100.0
    fig, ax = plt.subplots(figsize=(NATURE_W1 + 1.2, 0.36 * len(rate) + 1.2))
    y = np.arange(len(rate))
    colors = ["#D6604D" if r < 100 else "#3C5488" for r in rate.values]
    ax.hlines(y, 0, rate.values, color="#B8B8B8", linewidth=1.4, zorder=1)
    ax.scatter(rate.values, y, s=70, color=colors, zorder=2, edgecolor="white")
    for yi, r in zip(y, rate.values):
        ax.text(r + 1.5, yi, f"{r:.0f}%", va="center", fontsize=7)
    ax.axvline(100, color="#999999", ls="--", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([_short(c) for c in rate.index], fontsize=7)
    ax.set_xlim(0, 112)
    ax.set_xlabel("Pass rate across all poses (%)")
    ax.set_title("Per-check pass rate")
    save_fig(fig, str(out))
    plt.close(fig)


def plot_good_bad_dumbbell(checks: pd.DataFrame, expected: dict, out: Path):
    """图3: 好-vs-坏 构象 dumbbell。

    对每项检查,比较 expected=valid 组 与 expected=invalid 组的通过率,
    用哑铃连线展示二者落差——坏构象应在键长/键角/位阻等项明显掉下来,
    即"诚实基线":证明 PoseBusters 真的把坏 pose 抓出来了。
    """
    is_valid = pd.Series(expected).reindex(checks.index)
    good = checks[is_valid == "valid"]
    bad = checks[is_valid == "invalid"]
    if good.empty or bad.empty:
        # 无成对标签时退化为整体通过率单点图
        rate = checks.mean(axis=0) * 100
        fig, ax = plt.subplots(figsize=(NATURE_W1 + 1, 0.36 * len(rate) + 1.2))
        ax.scatter(rate.values, range(len(rate)), s=60, color="#3C5488")
        ax.set_yticks(range(len(rate)))
        ax.set_yticklabels([_short(c) for c in rate.index], fontsize=7)
        ax.set_xlabel("Pass rate (%)")
        ax.set_title("Per-check pass rate (no good/bad labels)")
        save_fig(fig, str(out))
        plt.close(fig)
        return

    gr = good.mean(axis=0) * 100
    br = bad.mean(axis=0) * 100
    # 按落差(good - bad)排序,落差大的排上面
    order = (gr - br).sort_values(ascending=True).index
    gr, br = gr[order], br[order]
    y = np.arange(len(order))
    fig, ax = plt.subplots(figsize=(NATURE_W1 + 1.6, 0.36 * len(order) + 1.3))
    ax.hlines(y, br.values, gr.values, color="#C7C7C7", linewidth=2.0, zorder=1)
    ax.scatter(br.values, y, s=70, color="#D6604D", zorder=3,
               edgecolor="white", label="bad conformers")
    ax.scatter(gr.values, y, s=70, color="#4DAF67", zorder=3,
               edgecolor="white", label="good conformers")
    ax.set_yticks(y)
    ax.set_yticklabels([_short(c) for c in order], fontsize=7)
    ax.set_xlim(-5, 112)
    ax.set_xlabel("Pass rate within group (%)")
    ax.set_title("Good vs. bad conformers — PoseBusters gatekeeping")
    ax.legend(loc="lower left", fontsize=8)
    save_fig(fig, str(out))
    plt.close(fig)


# ============================================================================
# 4) 主流程
# ============================================================================
def main():
    ap = argparse.ArgumentParser(
        description="PoseBusters 对接 pose 物理有效性面板 (556)")
    ap.add_argument("--input", type=str, default=None,
                    help="多构象 .sdf;每个构象当作一个 pose 体检。缺省用合成示例。")
    args = ap.parse_args()

    DIR_EXAMPLE.mkdir(parents=True, exist_ok=True)
    DIR_RESULTS.mkdir(parents=True, exist_ok=True)
    DIR_ASSETS.mkdir(parents=True, exist_ok=True)

    set_pub_style(base_size=11, palette="npg")

    # --- 选数据 ---
    if args.input:
        sdf = Path(args.input).resolve()
        print(f"[input] 使用用户 SDF: {sdf}")
    else:
        sdf = DIR_EXAMPLE / "demo_poses.sdf"
        if not sdf.exists():
            print("[data] 未找到示例 SDF,首次运行自动合成...")
            build_demo_sdf(sdf)
        else:
            print(f"[data] 复用已存在的示例 SDF: {sdf.name}")

    # --- 读 pose ---
    mols, ids, exp = load_poses(sdf)
    print(f"[load] 读入 {len(mols)} 个 pose")
    if not mols:
        raise SystemExit("没有可用 pose,退出。")

    # --- PoseBusters 体检 ---
    print("[step] 跑 PoseBusters 物理有效性检查...")
    checks = run_posebusters(mols, ids)
    expected = dict(zip(ids, exp))

    # --- 汇总 PB-valid 通过率(诚实基线:直接报百分比) ---
    per_pose_pass = checks.all(axis=1)
    pb_valid_rate = 100.0 * per_pose_pass.mean()
    per_check_rate = (checks.mean(axis=0) * 100.0).round(1)

    # 好/坏对照
    is_valid = pd.Series(expected)
    good_ids = [i for i in ids if expected.get(i) == "valid"]
    bad_ids = [i for i in ids if expected.get(i) == "invalid"]
    good_rate = 100.0 * checks.loc[good_ids].all(axis=1).mean() if good_ids else float("nan")
    bad_rate = 100.0 * checks.loc[bad_ids].all(axis=1).mean() if bad_ids else float("nan")

    print(f"[result] 整体 PB-valid 通过率(全部检查都过): {pb_valid_rate:.1f}%")
    if good_ids and bad_ids:
        print(f"[result] 诚实对照  好构象 PB-valid={good_rate:.0f}%  "
              f"坏构象 PB-valid={bad_rate:.0f}%  → 坏构象被守门拦下")

    # --- 写结果表 ---
    checks_out = checks.copy()
    checks_out.insert(0, "expected", [expected.get(i, "unknown") for i in checks.index])
    checks_out["pb_valid_all"] = per_pose_pass.values
    checks_out.to_csv(DIR_RESULTS / "posebusters_checks.csv", encoding="utf-8")
    per_check_rate.rename("pass_rate_pct").to_csv(
        DIR_RESULTS / "per_check_pass_rate.csv", encoding="utf-8")
    with open(DIR_RESULTS / "summary.txt", "w", encoding="utf-8") as fh:
        fh.write(f"n_poses\t{len(checks)}\n")
        fh.write(f"n_checks\t{checks.shape[1]}\n")
        fh.write(f"PB_valid_rate_pct\t{pb_valid_rate:.1f}\n")
        fh.write(f"good_conformer_PB_valid_pct\t{good_rate:.1f}\n")
        fh.write(f"bad_conformer_PB_valid_pct\t{bad_rate:.1f}\n")

    # --- 绘图 ---
    print("[step] 绘图(tick-heatmap / lollipop / dumbbell)...")
    plot_tick_heatmap(checks, DIR_ASSETS / "pose_check_heatmap")
    plot_passrate_lollipop(checks, DIR_ASSETS / "per_check_passrate_lollipop")
    plot_good_bad_dumbbell(checks, expected, DIR_ASSETS / "good_vs_bad_dumbbell")

    # --- 环境快照(依赖版本) ---
    import posebusters
    import rdkit
    import matplotlib as mpl
    with open(DIR_RESULTS / "session_info.txt", "w", encoding="utf-8") as fh:
        fh.write(f"python\t{sys.version.split()[0]}\n")
        fh.write(f"posebusters\t{posebusters.__version__}\n")
        fh.write(f"rdkit\t{rdkit.__version__}\n")
        fh.write(f"numpy\t{np.__version__}\n")
        fh.write(f"pandas\t{pd.__version__}\n")
        fh.write(f"matplotlib\t{mpl.__version__}\n")
        fh.write(f"seed\t{SEED}\n")

    print(f"[done] 结果 → {DIR_RESULTS}  图 → {DIR_ASSETS}")


if __name__ == "__main__":
    main()
