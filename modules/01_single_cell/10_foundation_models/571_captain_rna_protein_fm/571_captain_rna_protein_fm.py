"""571 · CAPTAIN — 配对 RNA+表面蛋白(CITE-seq)多模态基础模型 · 蛋白填补基准台。

本模块做的事:在配对 CITE-seq 数据上评测「从 RNA 预测表面蛋白丰度」这一 CAPTAIN 的
旗舰下游任务,并**始终先跑一条本机可跑的朴素基线**:
  B0 同源基因基线 (matched-gene) : 用蛋白对应基因的归一化表达直接当预测值。
  B1 岭回归基线   (PCA + Ridge)  : RNA → PCA(仅在训练集上拟合) → 每蛋白岭回归。
评测口径:按细胞划分 train/test,per-protein 在**留出细胞**上算 Pearson r / RMSE。
预处理(HVG、scaler、PCA、Ridge)一律只在训练集拟合再 transform 测试集,防数据泄漏。

CAPTAIN 本体路径是**守卫式封装**:CAPTAIN 没有 PyPI 包,也没有可 import 的稳定 API,
上游是「克隆仓库 + 下载权重 + 跑 downstream_tasks/<任务>/*.py 脚本」的用法。
本脚本因此只做环境探测并打印**从上游仓库实际读到的**命令,绝不伪造 `import captain`。

上游仓库 : https://github.com/iamjiboya/CAPTAIN
API 核对来源:本地克隆的上游源码(2026-07-21 逐行 grep 确认,非网页转述)
  - downstream_tasks/cell_surface_protein_prediction/zeroshot_genate.py
      L35 argparse.ArgumentParser · L38-49 全部 flag · L545/548 两个输出 pickle
  - 同目录 finetune.py / genate.py / Tutorial_Protein_Prediction[_Zero_shot].ipynb
  - token_dict/{vocab.json,csp_token_dict.pickle,csp_align_dict.pickle,
    human_mouse_align.pickle} · prior_knowledge/final_{human,mouse}_prior_knwo.npy
  - README.md L27(docker)L45/56(conda+pip)L71-76(权重 Drive 链接)· LICENSE(MIT)
论文 : Ji B, Hu T, Wang J, et al. CAPTAIN: a multimodal foundation model pretrained on
       co-assayed single-cell RNA and protein. Nat Commun 2026 May 7;17(1):6161.
       doi:10.1038/s41467-026-72882-y · PMID 42098152 · PMCID PMC13365403
       (PMID/DOI/卷期页/作者名单均经 NCBI E-utilities efetch 核实)
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).parent
sys.path.append(str(HERE.parents[2] / "_framework"))   # modules/_framework/pubstyle.py
from pubstyle import (CMAP_DIVERGE, NATURE_W1, NATURE_W2, pal,  # noqa: E402
                      save_fig, set_pub_style)

SEED = 20260720


# ============================================================================
# 数据读取与归一化
# ============================================================================
def load_data(rna_csv: Path, adt_csv: Path, map_csv: Path):
    """读入配对 CITE-seq 计数矩阵(行=细胞,列=基因/蛋白)与蛋白→同源基因映射。"""
    rna = pd.read_csv(rna_csv, index_col=0, comment="#")
    adt = pd.read_csv(adt_csv, index_col=0, comment="#")
    pmap = pd.read_csv(map_csv, comment="#")
    cells = rna.index.intersection(adt.index)
    if len(cells) == 0:
        sys.exit("RNA 与 ADT 矩阵没有共同细胞条码 —— 两个文件的行名必须是同一批细胞")
    return rna.loc[cells], adt.loc[cells], pmap


def lognorm(counts: pd.DataFrame, target: float = 1e4) -> pd.DataFrame:
    """CPM-like 文库大小归一 + log1p(scanpy 标准做法)。"""
    lib = counts.sum(axis=1).replace(0, np.nan)
    return np.log1p(counts.div(lib, axis=0).fillna(0) * target)


def clr(counts: pd.DataFrame) -> pd.DataFrame:
    """ADT 的标准 CLR(centered log-ratio)归一,按细胞跨蛋白做几何均值中心化。

    ★方向对应 Seurat 的 **margin = 2**(不是 1)。Seurat 文档把 margin 写成
    "normalize across features (1) or cells (2)",措辞有歧义;以实现为准:
    `Seurat:::CustomNormalize` 做 `apply(data, MARGIN=margin, FUN=...)`,而 Seurat
    的矩阵是 **features × cells**,故 MARGIN=1 = 逐**基因/蛋白**跨细胞,
    MARGIN=2 = 逐**细胞**跨蛋白。本函数是后者。
    (实测:本实现与 Seurat margin=2 的逐元素 Pearson r = 0.94,与 margin=1 只有 0.64。)

    ★公式也不与 Seurat 逐位相同。Seurat 的 CLR 是
    `log1p(x / exp(sum(log1p(x[x>0])) / length(x)))`(几何均值只用正值求和、却除以
    全长,外层是 log1p 比值);本函数用的是教科书 CLR:log1p 后按细胞减去均值,
    即 log((1+x_i) / geomean_j(1+x_j)),中心化后每个细胞行均值严格为 0。
    两者单调相关但不等价,不要在论文里写成"用 Seurat CLR"。
    """
    x = np.log1p(counts.astype(float))
    return x.sub(x.mean(axis=1), axis=0)


# ============================================================================
# 基线:B0 同源基因 / B1 PCA+Ridge
# ============================================================================
def run_baselines(rna: pd.DataFrame, adt: pd.DataFrame, pmap: pd.DataFrame,
                  n_pcs: int = 30, test_frac: float = 0.3, alpha: float = 10.0):
    """训练/测试按细胞切分;所有拟合(scaler/PCA/Ridge)只见训练集,防泄漏。"""
    from sklearn.decomposition import PCA
    from sklearn.linear_model import Ridge
    from sklearn.preprocessing import StandardScaler

    rng = np.random.default_rng(SEED)
    n = rna.shape[0]
    perm = rng.permutation(n)
    n_test = max(30, int(round(n * test_frac)))
    te, tr = perm[:n_test], perm[n_test:]

    X = lognorm(rna)
    Y = clr(adt)

    # ---- B1: 训练集上拟合 scaler → PCA → 每蛋白 Ridge ----------------------
    sc = StandardScaler().fit(X.values[tr])
    n_pcs = int(min(n_pcs, len(tr) - 1, X.shape[1]))
    pca = PCA(n_components=n_pcs, random_state=SEED).fit(sc.transform(X.values[tr]))
    Ztr, Zte = pca.transform(sc.transform(X.values[tr])), pca.transform(sc.transform(X.values[te]))
    ridge = Ridge(alpha=alpha, random_state=SEED).fit(Ztr, Y.values[tr])
    pred_ridge = ridge.predict(Zte)

    # ---- B0: 同源基因(z-score 到蛋白的尺度,只用训练集的均值/方差)---------
    cog = dict(zip(pmap["protein"], pmap["cognate_gene"]))
    pred_cog = np.full_like(pred_ridge, np.nan)
    for j, p in enumerate(Y.columns):
        g = cog.get(p)
        if g is None or g not in X.columns:
            continue
        gx = X[g].values
        mu, sd = gx[tr].mean(), gx[tr].std() or 1.0
        ymu, ysd = Y.values[tr, j].mean(), Y.values[tr, j].std() or 1.0
        pred_cog[:, j] = (gx[te] - mu) / sd * ysd + ymu

    # ---- 评测:per-protein 留出细胞 Pearson r / RMSE -------------------------
    def _score(pred):
        r, rmse = [], []
        for j in range(Y.shape[1]):
            yt, yp = Y.values[te, j], pred[:, j]
            if np.all(np.isnan(yp)) or np.std(yp) < 1e-12 or np.std(yt) < 1e-12:
                r.append(np.nan); rmse.append(np.nan); continue
            r.append(float(np.corrcoef(yt, yp)[0, 1]))
            rmse.append(float(np.sqrt(np.mean((yt - yp) ** 2))))
        return np.array(r), np.array(rmse)

    r_cog, rmse_cog = _score(pred_cog)
    r_rid, rmse_rid = _score(pred_ridge)

    metrics = pd.DataFrame({
        "protein": Y.columns,
        "cognate_gene": [cog.get(p, "") for p in Y.columns],
        "r_matched_gene": r_cog, "rmse_matched_gene": rmse_cog,
        "r_pca_ridge": r_rid, "rmse_pca_ridge": rmse_rid,
    })
    metrics["delta_r"] = metrics["r_pca_ridge"] - metrics["r_matched_gene"]

    obs = pd.DataFrame(Y.values[te], index=Y.index[te], columns=Y.columns)
    return {
        "metrics": metrics,
        "observed": obs,
        "pred_ridge": pd.DataFrame(pred_ridge, index=obs.index, columns=Y.columns),
        "pred_cognate": pd.DataFrame(pred_cog, index=obs.index, columns=Y.columns),
        "n_train": int(len(tr)), "n_test": int(len(te)), "n_pcs": n_pcs,
    }


# ============================================================================
# CAPTAIN 守卫式封装(不 import 不存在的包,只探测 + 打印真实命令)
# ============================================================================
CAPTAIN_REPO = "https://github.com/iamjiboya/CAPTAIN"
# 以下路径/文件名逐条核对过本地克隆的上游源码(2026-07-21),不是推测。
# 代码文件:存在即可
CAPTAIN_CODE_FILES = [
    "downstream_tasks/cell_surface_protein_prediction/zeroshot_genate.py",
    "downstream_tasks/cell_surface_protein_prediction/finetune.py",
    "downstream_tasks/cell_surface_protein_prediction/genate.py",
    "token_dict/vocab.json",
    "token_dict/csp_token_dict.pickle",
    "token_dict/csp_align_dict.pickle",
    "token_dict/human_mouse_align.pickle",
]
# 数据资产:★上游仓库里这些文件是「装着 Google Drive 链接的几十字节纯文本占位符」,
# 不是真文件(实测 prior_knowledge/*.npy = 86 B、各任务目录下 CAPTAIN_Base.pt = 88 B,
# 内容就是一行 https://drive.google.com/...)。因此不能只判存在,必须判是不是占位符。
CAPTAIN_DATA_ASSETS = ["prior_knowledge/final_human_prior_knwo.npy"]
CAPTAIN_WEIGHTS = "CAPTAIN_Base.pt"
_PLACEHOLDER_MAX_BYTES = 4096


def _is_placeholder(path: Path) -> bool:
    """上游用「内容为一行 Drive 链接的小文本文件」占位真正的权重/先验矩阵。"""
    try:
        if path.stat().st_size > _PLACEHOLDER_MAX_BYTES:
            return False
        head = path.read_bytes()[:200].lstrip()
        return head.startswith(b"http")
    except OSError:
        return False


def probe_captain(repo_dir: str | None):
    """探测本机是否具备跑 CAPTAIN 的条件,并返回一份诚实的状态报告。

    CAPTAIN **没有 PyPI 包**(上游仓库无 setup.py / pyproject.toml,`captain/` 目录下
    也没有 `__init__.py`),上游用法是克隆仓库后直接跑 downstream_tasks 下的脚本
    (README 的安装方式为 docker 镜像或 `conda create -n captain python==3.10.0` +
    `pip install -r requirements.txt && pip install scgpt`)。因此这里既不 `import
    captain`,也不假装有函数式 API。
    """
    rep: dict = {"repo": CAPTAIN_REPO, "pypi_package": False}
    try:
        import torch  # noqa: F401
        rep["torch"] = torch.__version__
        rep["cuda_available"] = bool(torch.cuda.is_available())
    except ImportError:
        rep["torch"] = None
        rep["cuda_available"] = False

    if not repo_dir:
        rep["status"] = "skipped"
        rep["reason"] = "未提供 --captain-repo;CAPTAIN 需先克隆仓库并下载权重"
        return rep
    root = Path(repo_dir)
    if not root.is_dir():
        rep["status"] = "skipped"
        rep["reason"] = f"--captain-repo 路径不存在: {root}"
        return rep
    missing = [f for f in CAPTAIN_CODE_FILES if not (root / f).exists()]
    # 数据资产:缺失 或 仍是 Drive 链接占位符,都算「没有」
    stub_assets = [f for f in CAPTAIN_DATA_ASSETS
                   if not (root / f).exists() or _is_placeholder(root / f)]
    found = sorted(root.rglob(CAPTAIN_WEIGHTS))
    weights = sorted(str(p.relative_to(root)) for p in found if not _is_placeholder(p))
    stub_weights = sorted(str(p.relative_to(root)) for p in found if _is_placeholder(p))
    rep["missing_files"] = missing
    rep["weights_found"] = weights
    rep["placeholder_weights"] = stub_weights     # 仓库自带的 Drive 链接占位符
    rep["placeholder_or_missing_assets"] = stub_assets
    if missing or stub_assets or not weights or not rep["cuda_available"]:
        rep["status"] = "skipped"
        rep["reason"] = (
            "仓库代码文件缺失" if missing else
            f"先验矩阵仍是 Drive 链接占位符或缺失: {stub_assets}" if stub_assets else
            ("仓库里的 CAPTAIN_Base.pt 只是 Google Drive 链接占位符(约 88 B),"
             "需按上游 README 从 Drive 下载真权重覆盖" if stub_weights else
             "未找到 CAPTAIN_Base.pt 权重") if not weights else
            "无可用 CUDA GPU;CAPTAIN 推理/微调在 CPU 上不现实")
    else:
        rep["status"] = "ready"
    # 下列 flag 全部逐条核对过 zeroshot_genate.py 的 argparse(第 38-49 行),
    # 输出文件名来自同文件第 545 / 548 行。
    rep["how_to_run"] = (
        "cd <repo>/downstream_tasks/cell_surface_protein_prediction && "
        "python zeroshot_genate.py --data_rna_path <rna_test.h5ad> "
        "--data_protein_path <adt_test.h5ad> --model_filename CAPTAIN_Base.pt "
        "--species human   # 零样本;微调走同目录 finetune.py。"
        "权重默认从 --load_model_dir(默认=脚本所在目录)加载;"
        "输出写到 --save_dir(默认=同目录 results/)下的 "
        "true_adt_data_scale.pickle / predicted_adt_scale.pickle"
    )
    return rep


# ============================================================================
# 出图(全部非条形图:dumbbell / 散点 / heatmap)
# ============================================================================
def plot_dumbbell(metrics: pd.DataFrame, outstem: Path):
    """Dumbbell:每个蛋白两种基线的留出集 Pearson r,连线看谁赢多少。"""
    import matplotlib.pyplot as plt

    d = metrics.dropna(subset=["r_pca_ridge"]).sort_values("r_pca_ridge")
    y = np.arange(len(d))
    c = pal(4, "npg")
    fig, ax = plt.subplots(figsize=(NATURE_W1 * 1.6, 0.30 * len(d) + 1.6))
    ax.hlines(y, d["r_matched_gene"], d["r_pca_ridge"], color="#BBBBBB", lw=1.6, zorder=1)
    ax.scatter(d["r_matched_gene"], y, s=46, color=c[1], zorder=2, label="Matched-gene baseline")
    ax.scatter(d["r_pca_ridge"], y, s=46, color=c[0], zorder=3, label="PCA + Ridge baseline")
    ax.axvline(0, color="black", lw=0.8, ls="--")
    ax.set_yticks(y); ax.set_yticklabels(d["protein"])
    ax.set_xlabel("Pearson r (held-out cells)")
    ax.set_title("RNA-to-protein imputation, per surface protein")
    ax.legend(loc="upper left")   # 左上区域在按 r 排序后必为空白,避免压点
    save_fig(fig, outstem); plt.close(fig)


def plot_scatter(res: dict, outstem: Path, top: int = 6):
    """散点小多图:留出细胞上 observed vs predicted(取 ridge 表现最好的若干蛋白)。"""
    import matplotlib.pyplot as plt

    m = res["metrics"].dropna(subset=["r_pca_ridge"]).nlargest(top, "r_pca_ridge")
    ncol = 3
    nrow = int(np.ceil(len(m) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(NATURE_W2, 2.3 * nrow), squeeze=False)
    c = pal(4, "npg")
    for ax, (_, row) in zip(axes.ravel(), m.iterrows()):
        p = row["protein"]
        x, yv = res["observed"][p].values, res["pred_ridge"][p].values
        ax.scatter(x, yv, s=9, alpha=0.55, color=c[0], edgecolors="none")
        lo, hi = np.nanmin([x, yv]), np.nanmax([x, yv])
        ax.plot([lo, hi], [lo, hi], color="#666666", lw=0.9, ls="--")
        ax.set_title(f"{p}  (r = {row['r_pca_ridge']:.2f})")
        ax.set_xlabel("Observed CLR"); ax.set_ylabel("Predicted")
    for ax in axes.ravel()[len(m):]:
        ax.axis("off")
    fig.suptitle("Held-out cells: observed vs predicted surface protein", y=1.01)
    fig.tight_layout()
    save_fig(fig, outstem); plt.close(fig)


def plot_heatmap(res: dict, outstem: Path):
    """Heatmap:预测蛋白 × 真实蛋白的相关矩阵。对角占优 = 预测有蛋白特异性,
    off-diagonal 亮 = 模型其实只学到了共享的细胞类型信号。"""
    import matplotlib.pyplot as plt

    O, P = res["observed"], res["pred_ridge"]
    keep = [p for p in O.columns if np.std(P[p].values) > 1e-12]
    M = np.array([[np.corrcoef(P[a].values, O[b].values)[0, 1] for b in keep] for a in keep])
    fig, ax = plt.subplots(figsize=(NATURE_W1 * 1.7, NATURE_W1 * 1.55))
    im = ax.imshow(M, cmap=CMAP_DIVERGE, vmin=-1, vmax=1)
    ax.set_xticks(range(len(keep))); ax.set_xticklabels(keep, rotation=90)
    ax.set_yticks(range(len(keep))); ax.set_yticklabels(keep)
    ax.set_xlabel("Observed protein"); ax.set_ylabel("Predicted protein")
    ax.set_title("Specificity check: predicted vs observed")
    fig.colorbar(im, ax=ax, shrink=0.75, label="Pearson r")
    save_fig(fig, outstem); plt.close(fig)


# ============================================================================
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rna", default=str(HERE / "example_data" / "citeseq_rna_counts.csv"))
    ap.add_argument("--adt", default=str(HERE / "example_data" / "citeseq_adt_counts.csv"))
    ap.add_argument("--protein-map", default=str(HERE / "example_data" / "protein_gene_map.csv"))
    ap.add_argument("--outdir", default=str(HERE / "results"))
    ap.add_argument("--assets", default=str(HERE / "assets"))
    ap.add_argument("--n-pcs", type=int, default=30)
    ap.add_argument("--test-frac", type=float, default=0.3)
    ap.add_argument("--alpha", type=float, default=10.0, help="Ridge 正则强度")
    ap.add_argument("--captain-repo", default=None,
                    help="本地 CAPTAIN 仓库克隆路径(给出后探测能否跑真模型)")
    a = ap.parse_args()

    out, assets = Path(a.outdir), Path(a.assets)
    out.mkdir(parents=True, exist_ok=True); assets.mkdir(parents=True, exist_ok=True)
    set_pub_style()

    print("Step 1 · 读入配对 CITE-seq 矩阵")
    rna, adt, pmap = load_data(Path(a.rna), Path(a.adt), Path(a.protein_map))
    print(f"       cells={rna.shape[0]}  genes={rna.shape[1]}  proteins={adt.shape[1]}")

    print("Step 2 · 基线:matched-gene 与 PCA+Ridge(仅训练集拟合,留出细胞评测)")
    res = run_baselines(rna, adt, pmap, n_pcs=a.n_pcs, test_frac=a.test_frac, alpha=a.alpha)
    m = res["metrics"]
    print(f"       train={res['n_train']}  test={res['n_test']}  PCs={res['n_pcs']}")
    print(f"       median r  matched-gene={np.nanmedian(m['r_matched_gene']):.3f}  "
          f"PCA+Ridge={np.nanmedian(m['r_pca_ridge']):.3f}")
    m.to_csv(out / "571_baseline_metrics.csv", index=False)
    res["pred_ridge"].to_csv(out / "571_pred_pca_ridge.csv")
    res["observed"].to_csv(out / "571_observed_clr.csv")

    print("Step 3 · 出图(dumbbell / 散点 / heatmap)")
    for stem in (out, assets):
        plot_dumbbell(m, stem / "571_protein_r_dumbbell")
        plot_scatter(res, stem / "571_obs_vs_pred_scatter")
        plot_heatmap(res, stem / "571_specificity_heatmap")

    print("Step 4 · CAPTAIN 本体探测(守卫式,不伪造 API)")
    cap = probe_captain(a.captain_repo)
    for k, v in cap.items():
        print(f"       {k}: {v}")

    summary = {
        "n_cells": int(rna.shape[0]), "n_genes": int(rna.shape[1]),
        "n_proteins": int(adt.shape[1]),
        "n_train": res["n_train"], "n_test": res["n_test"], "n_pcs": res["n_pcs"],
        "median_r_matched_gene": float(np.nanmedian(m["r_matched_gene"])),
        "median_r_pca_ridge": float(np.nanmedian(m["r_pca_ridge"])),
        "n_proteins_ridge_better": int((m["delta_r"] > 0).sum()),
        "captain": cap, "seed": SEED,
    }
    with open(out / "571_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=1, ensure_ascii=False, default=str)
    print(f"Done · results -> {out}")


if __name__ == "__main__":
    main()
