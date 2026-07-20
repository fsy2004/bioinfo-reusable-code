"""573 · Proseg — probabilistic cell segmentation for imaging-based spatial transcriptomics.

Proseg (Jones et al., Nat Methods 2025) treats segmentation as a *cell simulation*: it runs a
Bayesian model over a voxel field to reassign individual transcripts to cells, recovering
cell boundaries from the transcript point cloud rather than from a stained image.
It is a **Rust CLI**, not a Python package — this module therefore does two things:

  1) BASELINE (always runs, local deps only): the naive comparator that every vendor
     pipeline ships — nearest-nucleus expansion. Each transcript is assigned to the nearest
     nucleus centroid within a fixed radius, everything else is dropped. We sweep the radius,
     score against ground truth, and emit counts / cell metadata in the same spirit as
     proseg's own outputs. This is the floor any probabilistic method must beat.

  2) PROSEG PATH (--run-proseg): a guarded wrapper that locates the `proseg` executable and
     builds a command line from flags read off the official README. If the binary is absent
     it prints the real install command and exits cleanly. No fake Python API is invented.

Upstream repo : https://github.com/dcjones/proseg  (GPLv3, LICENSE.md)
API 核对来源  : proseg **3.2.0** 的 Rust 源码 src/main.rs 的 clap `struct Args` +
                src/output.rs 的 writer(不是只读 README);逐 flag 行号见
                build_proseg_cmd 的 docstring。
Paper         : Jones DC, Elz AE, Hadadianpour A, Ryu H, Glass DR, Newell EW.
                "Cell simulation as cell segmentation." Nat Methods 2025;22(6):1331-1342.
                doi:10.1038/s41592-025-02697-0 · PMID 40404994  (both verified)
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")
ASSETS = os.path.join(HERE, "assets")
EXAMPLE = os.path.join(HERE, "example_data", "transcripts.csv")

# 框架出图样式(Arial/期刊配色/矢量导出)
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..", "..", "_framework")))
from pubstyle import set_pub_style, pal, save_fig, CMAP_CONT, NATURE_W1, NATURE_W2  # noqa: E402

SEED = 573


# ---------------------------------------------------------------------------
# 合成数据:模拟 Xenium 式转录本点云(细胞核 + 胞质外扩 + 环境噪声)
# ---------------------------------------------------------------------------
def simulate_transcripts(n_cells=120, n_genes=24, n_types=4, seed=SEED) -> pd.DataFrame:
    """生成带 ground-truth 的转录本表。

    模拟真实成像空间组的三个痛点:①胞质转录本远离核、②细胞挨得近导致边界混淆、
    ③环境/扩散(diffusion)噪声转录本不属于任何细胞。Proseg 的卖点正是处理 ①③。
    """
    rng = np.random.default_rng(seed)
    # 细胞核中心:抖动网格,保证有的细胞挨得很近(边界难题)
    side = int(np.ceil(np.sqrt(n_cells)))
    gx, gy = np.meshgrid(np.arange(side), np.arange(side))
    centers = np.c_[gx.ravel(), gy.ravel()][:n_cells].astype(float) * 30.0
    centers += rng.normal(0, 6.0, centers.shape)

    cell_type = rng.integers(0, n_types, n_cells)
    # 每种细胞类型一套基因表达谱(狄利克雷,保证类型可分)
    type_profile = rng.dirichlet(np.full(n_genes, 0.35), n_types)
    genes = [f"GENE{i:02d}" for i in range(n_genes)]

    nucleus_r = rng.uniform(4.0, 6.5, n_cells)      # 核半径 (µm 量级)
    cell_sigma = rng.uniform(7.0, 12.0, n_cells)    # 胞质转录本扩散尺度
    n_tx = rng.integers(60, 160, n_cells)           # 每细胞转录本数

    rows = []
    for c in range(n_cells):
        k = int(n_tx[c])
        ang = rng.uniform(0, 2 * np.pi, k)
        # 半径:一部分落在核内,其余按半正态扩到胞质(长尾 → 远离核)
        in_nuc = rng.random(k) < 0.35
        rad = np.where(in_nuc,
                       nucleus_r[c] * np.sqrt(rng.random(k)),
                       np.abs(rng.normal(0, cell_sigma[c], k)) + nucleus_r[c])
        x = centers[c, 0] + rad * np.cos(ang)
        y = centers[c, 1] + rad * np.sin(ang)
        g = rng.choice(n_genes, k, p=type_profile[cell_type[c]])
        rows.append(pd.DataFrame({
            "x_location": x, "y_location": y,
            "feature_name": [genes[i] for i in g],
            "true_cell_id": c,
            "overlaps_nucleus": (rad <= nucleus_r[c]).astype(int),
        }))

    tx = pd.concat(rows, ignore_index=True)

    # 环境噪声转录本:不属于任何细胞 (true_cell_id = -1)
    n_amb = int(0.12 * len(tx))
    amb = pd.DataFrame({
        "x_location": rng.uniform(tx.x_location.min(), tx.x_location.max(), n_amb),
        "y_location": rng.uniform(tx.y_location.min(), tx.y_location.max(), n_amb),
        "feature_name": rng.choice(genes, n_amb),
        "true_cell_id": -1,
        "overlaps_nucleus": 0,
    })
    tx = pd.concat([tx, amb], ignore_index=True)
    tx = tx.sample(frac=1.0, random_state=seed).reset_index(drop=True)

    # 先验核分配:只有落在核内的转录本带 cell_id,其余为 -1(proseg 的输入要求)
    tx["cell_id"] = np.where(tx.overlaps_nucleus == 1, tx.true_cell_id, -1)
    tx.insert(0, "transcript_id", np.arange(len(tx)))
    tx["z_location"] = 0.0
    tx["true_cell_type"] = np.where(tx.true_cell_id >= 0,
                                    cell_type[tx.true_cell_id.clip(0)], -1)
    return tx[["transcript_id", "x_location", "y_location", "z_location",
               "feature_name", "overlaps_nucleus", "cell_id",
               "true_cell_id", "true_cell_type"]]


# ---------------------------------------------------------------------------
# 基线:最近核外扩 (nearest-nucleus expansion) —— 厂商默认做法
# ---------------------------------------------------------------------------
def nucleus_centroids(tx: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """由先验核分配 (cell_id >= 0) 求核质心。基线只能看到这些,看不到 ground truth。"""
    nuc = tx[tx.cell_id >= 0]
    cen = nuc.groupby("cell_id")[["x_location", "y_location"]].mean()
    return cen.index.to_numpy(), cen.to_numpy()


def expand_assign(tx: pd.DataFrame, ids: np.ndarray, cen: np.ndarray,
                  radius: float) -> np.ndarray:
    """把每个转录本分给半径内最近的核;超出半径 → 未分配 (-1)。"""
    from scipy.spatial import cKDTree
    tree = cKDTree(cen)
    d, j = tree.query(np.c_[tx.x_location, tx.y_location], k=1)
    return np.where(d <= radius, ids[j], -1)


def score(truth: np.ndarray, pred: np.ndarray) -> dict:
    """转录本级评分。核心指标只在「真属于某细胞」的转录本上算,避免噪声灌水。"""
    from sklearn.metrics import adjusted_rand_score
    real = truth >= 0
    correct = (pred == truth) & real
    return {
        "assigned_frac": float((pred >= 0).mean()),
        # 召回:真实细胞内转录本被正确归位的比例
        "recall_on_real": float(correct.sum() / max(real.sum(), 1)),
        # 精确:被分配出去的转录本里分对的比例(噪声被分配即算错)
        "precision": float(correct.sum() / max((pred >= 0).sum(), 1)),
        # 噪声污染:环境转录本被错误塞进细胞的比例
        "ambient_leak": float(((pred >= 0) & ~real).sum() / max((~real).sum(), 1)),
        "ari": float(adjusted_rand_score(truth, pred)),
    }


def run_baseline(tx: pd.DataFrame, outdir: str, radii) -> tuple[pd.DataFrame, float, np.ndarray]:
    ids, cen = nucleus_centroids(tx)
    truth = tx.true_cell_id.to_numpy()
    recs = []
    for r in radii:
        pred = expand_assign(tx, ids, cen, r)
        recs.append({"radius": r, **score(truth, pred)})
    sweep = pd.DataFrame(recs)
    # F1 选最优半径:这就是朴素法的软肋——结果强依赖一个手调的全局半径
    f1 = 2 * sweep.precision * sweep.recall_on_real / (sweep.precision + sweep.recall_on_real).clip(1e-9)
    sweep["f1"] = f1
    best_r = float(sweep.loc[f1.idxmax(), "radius"])
    best_pred = expand_assign(tx, ids, cen, best_r)

    sweep.to_csv(os.path.join(outdir, "baseline_radius_sweep.csv"), index=False)

    # 输出对齐 proseg 的产物语义:cell-by-gene counts + cell metadata
    a = tx.assign(pred_cell=best_pred).query("pred_cell >= 0")
    counts = pd.crosstab(a.pred_cell, a.feature_name)
    counts.to_csv(os.path.join(outdir, "baseline_cell_by_gene_counts.csv"))
    meta = a.groupby("pred_cell").agg(
        centroid_x=("x_location", "mean"), centroid_y=("y_location", "mean"),
        n_transcripts=("transcript_id", "size"),
        spread=("x_location", "std")).reset_index()
    meta.to_csv(os.path.join(outdir, "baseline_cell_metadata.csv"), index=False)
    return sweep, best_r, best_pred


# ---------------------------------------------------------------------------
# Proseg 守卫式封装:真实 CLI,不发明 Python API
# ---------------------------------------------------------------------------
INSTALL_HINT = (
    "proseg 是 Rust 二进制,不是 Python 包。安装:\n"
    "    cargo install proseg\n"
    "  或 bioconda(上游 README 的 badge,recipe 名 rust-proseg): conda install -c bioconda rust-proseg\n"
    "  或从源码: git clone https://github.com/dcjones/proseg && cd proseg && cargo build --release"
)


def build_proseg_cmd(exe: str, transcripts: str, outdir: str, preset: str = "xenium",
                     nthreads: int | None = None) -> list[str]:
    """拼 proseg 命令行。每个 flag 都对照过 proseg 3.2.0 的 Rust 源码 src/main.rs。

    用法形如 `proseg [arguments...] transcripts.csv.gz`(`transcript_csv` 是唯一位置参数,
    src/main.rs:54)。逐条核对(proseg 3.2.0,src/main.rs 的 clap `struct Args`):

      --xenium / --cosmx / --cosmx-micron / --merscope / --visiumhd  : L58 / L63 / L67 / L71 / L103
      --nthreads (-t, Option<usize>, 默认用满所有核)                  : L264-265
      --output-counts            : L371  → write_sparse_mtx(src/output.rs:140),永远 gzip
                                   matrix-market,不看扩展名,故 .mtx.gz 合法
      --output-cell-metadata     : L389  → 格式由扩展名推断(infer_format_from_filename,
      --output-transcript-metadata : L396  src/output.rs:95),**只认 .csv.gz / .csv / .parquet**,
                                   其他扩展名会 panic,故这里用 .csv.gz
      --output-cell-polygons     : L425  → write_cell_multipolygons(src/output.rs:854),
                                   gzip 后的 GeoJSON,扩展名任意
      --output-spatialdata       : L359,**默认值就是 "proseg-output.zarr"**

    最后一条是个坑:proseg 3.x 默认必写 spatialdata zarr。不显式指定的话它会落在当前工作
    目录而不是 outdir,所以这里显式指到 outdir 下。该 zarr 已存在时 proseg 会拒绝覆盖,
    需要 --overwrite(L366);本模块不代加,以免默默删掉别人的结果。

    模型参数(--ncomponents L232 默认 10 / --voxel-size L296 / --diffusion-probability L312 /
    --cell-compactness L269 默认 0.04 / --samples L252)确实存在,但本模块不代为固定,
    请以 `proseg --help` 为准。
    """
    cmd = [exe, f"--{preset}"]
    if nthreads:
        cmd += ["--nthreads", str(nthreads)]
    cmd += [
        "--output-counts", os.path.join(outdir, "counts.mtx.gz"),
        "--output-cell-metadata", os.path.join(outdir, "cell-metadata.csv.gz"),
        "--output-transcript-metadata", os.path.join(outdir, "transcript-metadata.csv.gz"),
        "--output-cell-polygons", os.path.join(outdir, "cell-polygons.geojson.gz"),
        # 显式给出,否则 3.x 的默认值会把 zarr 写进当前工作目录
        "--output-spatialdata", os.path.join(outdir, "proseg-output.zarr"),
        transcripts,
    ]
    return cmd


def run_proseg(transcripts: str, outdir: str, preset: str, nthreads, execute: bool) -> dict:
    exe = shutil.which("proseg")
    if exe is None:
        return {"status": "skipped", "reason": "`proseg` 不在 PATH", "install": INSTALL_HINT}
    cmd = build_proseg_cmd(exe, transcripts, outdir, preset, nthreads)
    if not execute:
        return {"status": "ready", "exe": exe, "command": " ".join(cmd)}
    p = subprocess.run(cmd, capture_output=True, text=True)
    return {"status": "ok" if p.returncode == 0 else "failed",
            "returncode": p.returncode, "command": " ".join(cmd),
            "stderr_tail": p.stderr[-800:]}


# ---------------------------------------------------------------------------
# 出图(全部英文标注;无条形图)
# ---------------------------------------------------------------------------
def figures(tx: pd.DataFrame, sweep: pd.DataFrame, best_r: float, pred: np.ndarray,
            outdir: str):
    import matplotlib.pyplot as plt
    set_pub_style(base_size=9)
    cols = pal(6, "npg")
    figs = {}

    # Fig 1 — 空间散点:ground truth vs 基线分割
    truth = tx.true_cell_id.to_numpy()
    rng = np.random.default_rng(SEED)
    ncell = int(truth.max()) + 1
    lut = rng.permutation(ncell)
    fig, axes = plt.subplots(1, 2, figsize=(NATURE_W2, NATURE_W2 / 2.15))
    for ax, lab, v in ((axes[0], "Ground truth", truth), (axes[1],
                       f"Nucleus expansion (r = {best_r:g})", pred)):
        m = v >= 0
        ax.scatter(tx.x_location[~m], tx.y_location[~m], s=1.2, c="#D9D9D9",
                   linewidths=0, rasterized=True, label="unassigned / ambient")
        ax.scatter(tx.x_location[m], tx.y_location[m], s=1.6,
                   c=lut[v[m]] % 20, cmap="tab20", linewidths=0, rasterized=True)
        ax.set_title(lab, fontsize=9)
        ax.set_xlabel("x (µm)"); ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([])
    axes[0].set_ylabel("y (µm)")
    axes[1].legend(loc="upper right", markerscale=6, fontsize=6)
    save_fig(fig, os.path.join(ASSETS, "fig1_segmentation_map"))
    plt.close(fig); figs["fig1_segmentation_map.png"] = "spatial scatter, truth vs baseline"

    # Fig 2 — 半径扫描:朴素法对全局半径的敏感性(线+点,非条形)
    fig, ax = plt.subplots(figsize=(NATURE_W1 * 1.25, NATURE_W1))
    for k, lab, c in (("recall_on_real", "Recall (true-cell tx)", cols[0]),
                      ("precision", "Precision", cols[1]),
                      ("ari", "Adjusted Rand index", cols[2]),
                      ("ambient_leak", "Ambient leak-in", cols[3])):
        ax.plot(sweep.radius, sweep[k], "-o", ms=3.5, lw=1.2, color=c, label=lab)
    ax.axvline(best_r, ls="--", lw=0.9, color="#666666")
    # 标注放在轴内左侧竖排,避免与标题重叠
    ax.text(best_r - 0.6, 0.02, f"best F1  r = {best_r:g}", rotation=90,
            ha="right", va="bottom", fontsize=6, color="#444444")
    ax.set_xlabel("Nucleus expansion radius (µm)"); ax.set_ylabel("Score")
    ax.set_ylim(-0.02, 1.05); ax.legend(fontsize=6, loc="center right")
    ax.set_title("Naive baseline depends on one hand-tuned radius", fontsize=8)
    save_fig(fig, os.path.join(ASSETS, "fig2_radius_sweep"))
    plt.close(fig); figs["fig2_radius_sweep.png"] = "radius sweep, line+dot"

    # Fig 3 — 每细胞转录本数:violin + 抖动点(raincloud 风格)
    t_cnt = pd.Series(truth[truth >= 0]).value_counts().sort_index()
    p_cnt = pd.Series(pred[pred >= 0]).value_counts().reindex(t_cnt.index).fillna(0)
    fig, ax = plt.subplots(figsize=(NATURE_W1, NATURE_W1))
    parts = ax.violinplot([t_cnt.values, p_cnt.values], positions=[0, 1],
                          showextrema=False, widths=0.75)
    for b, c in zip(parts["bodies"], (cols[4], cols[1])):
        b.set_facecolor(c); b.set_alpha(0.45); b.set_edgecolor("black"); b.set_linewidth(0.7)
    jr = np.random.default_rng(SEED)
    for i, (v, c) in enumerate(((t_cnt.values, cols[4]), (p_cnt.values, cols[1]))):
        ax.scatter(i + jr.normal(0, 0.045, len(v)), v, s=5, color=c,
                   alpha=0.7, linewidths=0)
        ax.hlines(np.median(v), i - 0.22, i + 0.22, color="black", lw=1.3)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Ground truth", "Baseline"])
    ax.set_ylabel("Transcripts per cell")
    # 标题据实描述:在 F1 最优半径下基线是「偏多」还是「偏少」由数据决定,不预设结论
    dm = float(np.median(p_cnt.values) - np.median(t_cnt.values))
    ax.set_title(f"Per-cell yield at r = {best_r:g} (median {dm:+.0f} vs truth)", fontsize=8)
    save_fig(fig, os.path.join(ASSETS, "fig3_counts_per_cell"))
    plt.close(fig); figs["fig3_counts_per_cell.png"] = "violin + jitter"

    # Fig 4 — 每细胞计数散点 truth vs baseline(y=x 参考线)
    fig, ax = plt.subplots(figsize=(NATURE_W1, NATURE_W1))
    lim = max(t_cnt.max(), p_cnt.max()) * 1.05
    ax.plot([0, lim], [0, lim], ls="--", lw=0.9, color="#888888", zorder=0)
    sc = ax.scatter(t_cnt.values, p_cnt.values, c=p_cnt.values / t_cnt.values.clip(1),
                    s=16, cmap=CMAP_CONT, linewidths=0.3, edgecolors="white")
    plt.colorbar(sc, ax=ax, label="Recovered fraction", shrink=0.8)
    from scipy.stats import pearsonr
    r = pearsonr(t_cnt.values, p_cnt.values)[0]
    ax.set_xlabel("True transcripts per cell"); ax.set_ylabel("Baseline transcripts per cell")
    ax.set_title(f"Per-cell count recovery (Pearson r = {r:.2f})", fontsize=8)
    ax.set_xlim(0, lim); ax.set_ylim(0, lim)
    save_fig(fig, os.path.join(ASSETS, "fig4_count_recovery"))
    plt.close(fig); figs["fig4_count_recovery.png"] = "scatter vs identity line"

    # Fig 5 — 细胞类型 × 基因 表达热图(基线分割后的下游可用性)
    a = tx.assign(pred_cell=pred).query("pred_cell >= 0")
    ct = tx[tx.true_cell_id >= 0].groupby("true_cell_id").true_cell_type.first()
    cbg = pd.crosstab(a.pred_cell, a.feature_name)
    cbg = cbg.div(cbg.sum(1), axis=0)
    grp = cbg.groupby(ct.reindex(cbg.index).values).mean()
    z = (grp - grp.mean()) / grp.std().replace(0, 1)
    fig, ax = plt.subplots(figsize=(NATURE_W2 * 0.8, NATURE_W1 * 0.75))
    im = ax.imshow(z.values, cmap="RdBu_r", aspect="auto", vmin=-2, vmax=2)
    ax.set_xticks(range(z.shape[1])); ax.set_xticklabels(z.columns, rotation=90, fontsize=5)
    ax.set_yticks(range(z.shape[0]))
    ax.set_yticklabels([f"Type {int(i)}" for i in z.index], fontsize=7)
    plt.colorbar(im, ax=ax, label="Row-scaled mean fraction", shrink=0.85)
    ax.set_title("Cell-type expression profiles after baseline segmentation", fontsize=8)
    save_fig(fig, os.path.join(ASSETS, "fig5_celltype_heatmap"))
    plt.close(fig); figs["fig5_celltype_heatmap.png"] = "heatmap"

    return figs


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", default=EXAMPLE, help="转录本表 CSV(默认 example_data/)")
    ap.add_argument("--outdir", default=RESULTS)
    ap.add_argument("--radii", default="6,9,12,15,18,21,25,30",
                    help="基线外扩半径扫描列表 (µm)")
    ap.add_argument("--run-proseg", action="store_true", help="尝试调用真实 proseg 二进制")
    ap.add_argument("--proseg-execute", action="store_true",
                    help="不仅打印命令,还真的执行(需 proseg 在 PATH 且输入为其支持格式)")
    ap.add_argument("--preset", default="xenium",
                    choices=["xenium", "cosmx", "cosmx-micron", "merscope", "visiumhd"])
    ap.add_argument("--nthreads", type=int, default=None)
    a = ap.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    os.makedirs(ASSETS, exist_ok=True)
    np.random.seed(SEED)

    print("[573] Step 1 · 读取转录本点云")
    if not os.path.exists(a.input):
        print(f"       {a.input} 不存在,生成合成示例数据")
        os.makedirs(os.path.dirname(a.input), exist_ok=True)
        simulate_transcripts().to_csv(a.input, index=False)
    tx = pd.read_csv(a.input)
    need = {"x_location", "y_location", "feature_name", "cell_id"}
    miss = need - set(tx.columns)
    if miss:
        sys.exit(f"输入缺列 {sorted(miss)};proseg 要求转录本表带位置 + 先验核/细胞分配")
    has_truth = "true_cell_id" in tx.columns
    print(f"       {len(tx):,} transcripts · {tx.feature_name.nunique()} genes · "
          f"{(tx.cell_id >= 0).sum():,} 带先验核分配")
    if not has_truth:
        sys.exit("本模块的基线评分需要 true_cell_id 列(合成数据自带);真实数据请只用 --run-proseg 路径")

    print("[573] Step 2 · 基线:最近核外扩 + 半径扫描")
    radii = [float(x) for x in a.radii.split(",")]
    sweep, best_r, pred = run_baseline(tx, a.outdir, radii)
    best = sweep.loc[sweep.f1.idxmax()]
    for k in ("radius", "assigned_frac", "recall_on_real", "precision", "ambient_leak", "ari", "f1"):
        print(f"       {k:>16}: {best[k]:.3f}")

    print("[573] Step 3 · Proseg 路径")
    pinfo = ({"status": "not requested (--run-proseg)"} if not a.run_proseg
             else run_proseg(a.input, a.outdir, a.preset, a.nthreads, a.proseg_execute))
    for k, v in pinfo.items():
        print(f"       {k}: {v}")

    print("[573] Step 4 · 出图")
    figs = figures(tx, sweep, best_r, pred, a.outdir)
    for f, d in figs.items():
        print(f"       {f}  ({d})")

    with open(os.path.join(a.outdir, "573_summary.json"), "w", encoding="utf-8") as fh:
        json.dump({"n_transcripts": int(len(tx)), "best_radius": best_r,
                   "baseline_best": {k: float(best[k]) for k in sweep.columns},
                   "proseg": pinfo, "figures": list(figs)}, fh, indent=1,
                  ensure_ascii=False, default=str)
    print(f"[573] done · {a.outdir}")


if __name__ == "__main__":
    main()
