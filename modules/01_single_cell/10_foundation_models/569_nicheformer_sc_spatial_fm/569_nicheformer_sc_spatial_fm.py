# -*- coding: utf-8 -*-
# =============================================================================
# 569 · Nicheformer —— 单细胞 + 空间联合预训练基础模型(空间感知嵌入 / 标签·生态位转移)
# -----------------------------------------------------------------------------
# Nicheformer (Tejada-Lapuerta, Schaar et al., Nat Methods 2025) 在 SpatialCorpus-110M
# 上联合预训练 —— 摘要原文:"over 57 million dissociated and 53 million spatially
# resolved cells across 73 tissues"(即 5700 万解离 + 5300 万空间,共约 1.1 亿细胞;
# 5700 万只是解离那一半,别写成总量)。产出「带空间上下文」的细胞嵌入,用于
# 生态位(niche)标签预测与解离↔空间之间的标签迁移。
#
# ★本模块的定位(诚实):完整 Nicheformer 需要下载官方 checkpoint(Mendeley Data)
#   + tokenization 流程 + GPU,本机不装包,因此走【守卫式引用封装】(--run-nicheformer)。
#   模块主体是一条**本机零改动可跑的诚实基线**:
#     intrinsic  = 细胞自身表达 PCA(无空间上下文)
#     niche-aware= 自身表达 PCA ⊕ 空间 kNN 邻域均值表达 PCA(线性"邻域 token"对照)
#   两者在同一套交叉验证下比 niche / cell-type 标签的 macro-F1,再做解离参考→空间
#   query 的跨模态标签迁移。任何"基础模型更好"的说法都必须先打赢这条线性基线
#   (Genome Biol 2025 zero-shot 打不过 PCA;Nat Methods 2025 扰动 DL 打不过线性基线)。
#
# Turnkey: python 569_nicheformer_sc_spatial_fm.py       (默认读 example_data/ → results/+assets/)
#          换数据: --sp_expr a.csv --sp_meta b.csv --ref_expr c.csv --ref_meta d.csv
# 复用 _framework/pubstyle.py;图中文字英文,注释中文;不用条形图。
#
# 论文  : Tejada-Lapuerta A, Schaar AC, Gutgesell R, ... Theis FJ.
#         Nicheformer: a foundation model for single-cell and spatial omics.
#         Nat Methods. 2025. doi:10.1038/s41592-025-02814-z · PMID 41168487  (已核实)
# 仓库  : https://github.com/theislab/nicheformer
# =============================================================================
from __future__ import annotations

import argparse
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

# ---- 定位脚本目录 + 载入顶刊绘图框架(向上搜 _framework)--------------------
HERE = Path(__file__).resolve().parent
for up in [HERE, *HERE.parents]:
    if (up / "_framework" / "pubstyle.py").exists():
        sys.path.insert(0, str(up / "_framework"))
        break
try:
    from pubstyle import (set_pub_style, save_fig, pal, panel_labels,
                          NATURE_W1, NATURE_W2, CMAP_CONT)
except Exception:                       # 框架缺失时最小降级,不影响分析
    def set_pub_style(*a, **k): pass
    def save_fig(fig, f, dpi=300):
        Path(f).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(str(f) + ".pdf"); fig.savefig(str(f) + ".png", dpi=dpi)
    def pal(n=None, name="npg"):
        base = ["#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F"]
        return base[:n] if n else base
    def panel_labels(*a, **k): pass
    NATURE_W1, NATURE_W2, CMAP_CONT = 3.5, 7.0, "viridis"

SEED = 42
EX = HERE / "example_data"


# =============================================================================
# 1. 数据读取与预处理
# =============================================================================
def load_inputs(sp_expr, sp_meta, ref_expr, ref_meta):
    """读表达矩阵(cells × genes)与元数据;返回对齐后的 4 个对象。"""
    Xs = pd.read_csv(sp_expr, index_col=0)
    Ms = pd.read_csv(sp_meta).set_index("cell_id").loc[Xs.index]
    Xr = pd.read_csv(ref_expr, index_col=0)
    Mr = pd.read_csv(ref_meta).set_index("cell_id").loc[Xr.index]
    shared = [g for g in Xs.columns if g in set(Xr.columns)]   # 只用共有基因
    if len(shared) < 5:
        sys.exit("空间与参考数据共有基因 < 5,无法做跨模态迁移")
    return Xs[shared], Ms, Xr[shared], Mr


def norm_log(X: pd.DataFrame, target: float = 1e4) -> np.ndarray:
    """CPM-like 归一化 + log1p —— 与 scanpy normalize_total/log1p 等价的最小实现。"""
    A = X.to_numpy(dtype=float)
    tot = A.sum(1, keepdims=True)
    tot[tot == 0] = 1.0
    return np.log1p(A / tot * target)


# =============================================================================
# 2. 两种表征:intrinsic(无空间)vs niche-aware(线性邻域 token)
# =============================================================================
def spatial_knn_mean(L: np.ndarray, coords: np.ndarray, k: int) -> np.ndarray:
    """空间 kNN 邻域均值表达(排除自身)。

    这是 Nicheformer「niche token」(X_niche_0..4:每个细胞的 k 个空间近邻)的
    线性对照版 —— 论文用 transformer 编码邻域,这里用邻域均值,故意做得朴素。
    """
    from sklearn.neighbors import NearestNeighbors
    nn = NearestNeighbors(n_neighbors=k + 1).fit(coords)
    idx = nn.kneighbors(coords, return_distance=False)[:, 1:]   # 去掉自身
    return L[idx].mean(axis=1)


def build_reps(L_sp: np.ndarray, coords: np.ndarray, n_pcs: int, k: int, seed: int):
    """返回 {'intrinsic': ..., 'niche-aware': ...} 两套嵌入。"""
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    sc_i = StandardScaler().fit(L_sp)
    p_i = PCA(n_components=n_pcs, random_state=seed).fit(sc_i.transform(L_sp))
    E_int = p_i.transform(sc_i.transform(L_sp))

    Lnb = spatial_knn_mean(L_sp, coords, k)
    sc_n = StandardScaler().fit(Lnb)
    p_n = PCA(n_components=n_pcs, random_state=seed).fit(sc_n.transform(Lnb))
    E_nb = p_n.transform(sc_n.transform(Lnb))

    # 两块各自 z-score 后拼接,避免某一块方差主导
    z = lambda M: (M - M.mean(0)) / (M.std(0) + 1e-8)
    return {"intrinsic": E_int, "niche-aware": np.hstack([z(E_int), z(E_nb)])}


# =============================================================================
# 3. 评估:同折交叉验证下比较两种表征
# =============================================================================
def cv_scores(E: np.ndarray, y: np.ndarray, folds: int, seed: int, k: int = 15):
    """分层 K 折 kNN 分类的 macro-F1(每折一个值)+ 汇总的混淆矩阵。

    防泄漏:标签只在 fit 折内使用;表征本身是无监督构建的(PCA/邻域均值不看标签)。
    """
    from sklearn.model_selection import StratifiedKFold
    from sklearn.neighbors import KNeighborsClassifier
    from sklearn.metrics import f1_score, confusion_matrix

    labs = np.unique(y)
    skf = StratifiedKFold(n_splits=folds, shuffle=True, random_state=seed)
    f1s, pred_all = [], np.empty_like(y)
    for tr, te in skf.split(E, y):
        clf = KNeighborsClassifier(n_neighbors=k).fit(E[tr], y[tr])
        p = clf.predict(E[te])
        pred_all[te] = p
        f1s.append(f1_score(y[te], p, average="macro", labels=labs))
    cm = confusion_matrix(y, pred_all, labels=labs)
    return np.array(f1s), cm, labs


def cross_modality_transfer(L_ref, y_ref, L_sp, y_sp, n_pcs, seed, k=15):
    """解离参考 → 空间 query 的细胞类型标签迁移(Nicheformer 的核心应用场景之一)。

    朴素做法:在参考集上拟合 scaler+PCA,把 query 投到同一空间,kNN 传标签。
    这条线故意不做批次校正 —— 它是"不校正也能到多少"的地板。
    """
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    from sklearn.neighbors import KNeighborsClassifier
    from sklearn.metrics import f1_score

    sc = StandardScaler().fit(L_ref)
    p = PCA(n_components=n_pcs, random_state=seed).fit(sc.transform(L_ref))
    clf = KNeighborsClassifier(n_neighbors=k).fit(p.transform(sc.transform(L_ref)), y_ref)
    pred = clf.predict(p.transform(sc.transform(L_sp)))
    return {"accuracy": float((pred == y_sp).mean()),
            "macro_f1": float(f1_score(y_sp, pred, average="macro"))}


# =============================================================================
# 4. 出图(lollipop / slopegraph / 散点 / heatmap —— 不用条形图)
# =============================================================================
def fig_tissue_and_embedding(coords, niche, ct, reps, outdir):
    """Fig1:组织切片(niche / cell type)+ 两种嵌入的 PC 散点(按 niche 着色)。"""
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(NATURE_W2, 6.2))
    nl, cl = np.unique(niche), np.unique(ct)
    cmap_n = dict(zip(nl, pal(len(nl), "npg")))
    cmap_c = dict(zip(cl, pal(len(cl), "jama")))

    for lab in nl:
        m = niche == lab
        axes[0, 0].scatter(coords[m, 0], coords[m, 1], s=7, c=cmap_n[lab], label=lab, lw=0)
    axes[0, 0].set(title="Tissue section: spatial niche", xlabel="x (um)", ylabel="y (um)")
    # 图例放到画布外,避免压住散点(顶刊制图习惯)
    axes[0, 0].legend(markerscale=2, fontsize=6, loc="upper left",
                      bbox_to_anchor=(1.0, 1.0), borderaxespad=0.2)

    for lab in cl:
        m = ct == lab
        axes[0, 1].scatter(coords[m, 0], coords[m, 1], s=7, c=cmap_c[lab], label=lab, lw=0)
    axes[0, 1].set(title="Tissue section: cell type", xlabel="x (um)", ylabel="y (um)")
    axes[0, 1].legend(markerscale=2, fontsize=6, loc="upper left",
                      bbox_to_anchor=(1.0, 1.0), borderaxespad=0.2)

    for ax, name in zip(axes[1], ["intrinsic", "niche-aware"]):
        E = reps[name]
        for lab in nl:
            m = niche == lab
            ax.scatter(E[m, 0], E[m, 1], s=7, c=cmap_n[lab], label=lab, lw=0)
        ax.set(title=f"{name} embedding (PC1/PC2)", xlabel="PC1", ylabel="PC2")
    panel_labels(axes.ravel())
    fig.tight_layout()
    save_fig(fig, outdir / "fig1_tissue_and_embeddings")
    plt.close(fig)


def fig_slopegraph(res, outdir):
    """Fig2:每折 macro-F1 的 slopegraph(intrinsic → niche-aware),两个任务分色。"""
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(NATURE_W1 + 0.6, 3.6))
    cols = pal(len(res), "npg")
    for (task, d), c in zip(res.items(), cols):
        a, b = d["intrinsic"], d["niche-aware"]
        for i in range(len(a)):
            ax.plot([0, 1], [a[i], b[i]], color=c, alpha=0.35, lw=1.0, zorder=1)
            ax.scatter([0, 1], [a[i], b[i]], s=16, color=c, zorder=2, lw=0)
        ax.plot([0, 1], [a.mean(), b.mean()], color=c, lw=2.6, zorder=3, label=task)
        ax.scatter([0, 1], [a.mean(), b.mean()], s=60, color=c,
                   edgecolor="white", lw=1.2, zorder=4)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["intrinsic\n(expression PCA)",
                                               "niche-aware\n(+ spatial kNN context)"])
    ax.set_xlim(-0.28, 1.28); ax.set_ylabel("macro-F1 (5-fold CV)")
    ax.set_title("Does spatial context help?\nlinear baseline, per-fold")
    ax.legend(loc="center left", fontsize=7)
    fig.tight_layout()
    save_fig(fig, outdir / "fig2_representation_slopegraph")
    plt.close(fig)


def fig_confusion(cms, outdir):
    """Fig3:niche 标签的行归一化混淆矩阵热图(两种表征并排)。"""
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, len(cms), figsize=(NATURE_W2, 3.2))
    axes = np.atleast_1d(axes)
    for pi, (ax, (name, (cm, labs))) in enumerate(zip(axes, cms.items())):
        M = cm / cm.sum(1, keepdims=True)
        im = ax.imshow(M, cmap=CMAP_CONT, vmin=0, vmax=1)
        ax.set_xticks(range(len(labs))); ax.set_xticklabels(labs, rotation=45, ha="right")
        ax.set_yticks(range(len(labs))); ax.set_yticklabels(labs)
        # 只有最左侧 panel 保留 y 轴标题,否则会压到左邻 panel 的刻度标签
        ax.set(title=f"niche recall · {name}", xlabel="predicted",
               ylabel="true" if pi == 0 else "")
        for i in range(len(labs)):
            for j in range(len(labs)):
                ax.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center", fontsize=7,
                        color="white" if M[i, j] < 0.6 else "black")
    fig.colorbar(im, ax=axes.tolist(), fraction=0.025, label="row-normalised fraction")
    save_fig(fig, outdir / "fig3_niche_confusion_heatmap")
    plt.close(fig)


# =============================================================================
# 5. 守卫式引用封装:真正的 Nicheformer 路径
# =============================================================================
def run_nicheformer(checkpoint: str | None):
    """Nicheformer 官方路径。需要 pip install -e (源码安装) + 官方 checkpoint + GPU。

    以下 API 名称逐个核对过上游源码(theislab/nicheformer @ main),括号内为源码位置:

    公共 API(经 nicheformer.models 导出,可直接 import):
      nicheformer.models.Nicheformer            models/__init__.py:1 → models/_nicheformer.py:15
        pl.LightningModule;__init__ 见 _nicheformer.py:16-36
      nicheformer.models.NicheformerFineTune    models/__init__.py:2 → models/_nicheformer_fine_tune.py:18
      Nicheformer.get_embeddings(batch, layer=-1, with_context=False)
        _nicheformer.py:272;layer<0 表示倒数第 layer 层;with_context=False 时丢掉前 3 个
        上下文 token(embeddings[:, 3:, :]),再对 token 维求均值;返回 (batch, dim_model) 张量

    ★非公共 API(源码在,但**不是可调用的库接口**,原样调用必失败):
      src/nicheformer/_embeddings.py::get_embeddings_model(config)   _embeddings.py:10
        config 键 'checkpoint_path'/'fine_tuned_checkpoint_path'/'organ'  _embeddings.py:17,26
          (取值样例见 config_files/_config_embeddings.py:2-4)
        merlin 列名 'X'/'X_niche_0'..'X_niche_4'/'idx'/'assay'/'specie'/'modality'
          _embeddings.py:32,48;这些列随后成为 batch 键 _embeddings.py:109-122
        输出 embeddings.npy(512 维,_embeddings.py:100,138)+ metadata_embeddings.csv:139
        ⚠ 三点实测问题:① 顶部是 `from models._nicheformer import ...`(非包内相对导入),
          装成包后 import 不到;② path_organ 全是作者集群的 /lustre/... 硬编码路径;
          ③ 同目录入口 get_embeddings.py:3 导入的是 `get_embeddings_organ`,而 _embeddings.py
          里并不存在这个名字。故它只能当**参考实现**读,不能当 API 用。
    数据必须先经官方 tokenization notebook 处理:notebooks/tokenization/ 下实际只有 3 个 —
    cosmx_human_lung / xenium_human_lung / scRNAseq_mouse_brain_haviv(没有 merfish 版)。
    """
    try:
        import nicheformer                                     # noqa: F401
        from nicheformer.models import Nicheformer             # noqa: F401
    except Exception as e:
        return {"status": "skipped",
                "reason": f"nicheformer 未安装 ({type(e).__name__}: {e})",
                "install": ("git clone https://github.com/theislab/nicheformer && "
                            "cd nicheformer && pip install -e .   (Python >=3.9)"),
                "weights": ("官方 checkpoint 托管于 Mendeley Data,见仓库 README "
                            "https://github.com/theislab/nicheformer  "
                            "(注意:不是 HuggingFace)")}
    try:
        import torch
        gpu = torch.cuda.is_available()
    except Exception:
        gpu = False
    if not checkpoint:
        return {"status": "skipped", "reason": "未提供 --checkpoint(官方预训练权重路径)",
                "gpu_available": gpu}
    if not Path(checkpoint).exists():
        return {"status": "skipped", "reason": f"checkpoint 不存在: {checkpoint}", "gpu_available": gpu}
    return {"status": "ready",
            "gpu_available": gpu,
            "checkpoint": checkpoint,
            "next": ("按官方 notebooks/ 走:tokenization → Nicheformer.load_from_checkpoint(ckpt) "
                     "→ get_embeddings(batch, layer=-1) → 下游 niche/label 转移;"
                     "签名与 config 完整键值以官方 notebook 为准,本模块未固定")}


# =============================================================================
# 6. 主流程
# =============================================================================
def main():
    ap = argparse.ArgumentParser(description="569 · Nicheformer 空间感知嵌入 + 诚实线性基线")
    ap.add_argument("--sp_expr", default=str(EX / "spatial_query_expression.csv"))
    ap.add_argument("--sp_meta", default=str(EX / "spatial_query_meta.csv"))
    ap.add_argument("--ref_expr", default=str(EX / "dissociated_reference_expression.csv"))
    ap.add_argument("--ref_meta", default=str(EX / "dissociated_reference_meta.csv"))
    ap.add_argument("--niche_key", default="niche")
    ap.add_argument("--label_key", default="cell_type")
    ap.add_argument("--n_pcs", type=int, default=15)
    ap.add_argument("--k_spatial", type=int, default=15, help="空间邻域 token 的近邻数")
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--outdir", default=str(HERE / "results"))
    ap.add_argument("--assets", default=str(HERE / "assets"))
    ap.add_argument("--run-nicheformer", dest="run_nf", action="store_true",
                    help="尝试真正的 Nicheformer 路径(需源码安装 + 官方 checkpoint + GPU)")
    ap.add_argument("--checkpoint", default=None, help="官方预训练权重路径")
    a = ap.parse_args()

    np.random.seed(a.seed)
    set_pub_style(base_size=9)
    outdir, assets = Path(a.outdir), Path(a.assets)
    outdir.mkdir(parents=True, exist_ok=True); assets.mkdir(parents=True, exist_ok=True)

    print("Step 1 · 读取空间 query 与解离参考数据")
    Xs, Ms, Xr, Mr = load_inputs(a.sp_expr, a.sp_meta, a.ref_expr, a.ref_meta)
    L_sp, L_ref = norm_log(Xs), norm_log(Xr)
    coords = Ms[["x", "y"]].to_numpy(float)
    niche = Ms[a.niche_key].to_numpy(str)
    ct_sp = Ms[a.label_key].to_numpy(str)
    ct_ref = Mr[a.label_key].to_numpy(str)
    print(f"       spatial {L_sp.shape} · reference {L_ref.shape} · "
          f"niches {len(np.unique(niche))} · cell types {len(np.unique(ct_sp))}")

    print(f"Step 2 · 构建两种表征(n_pcs={a.n_pcs}, spatial k={a.k_spatial})")
    reps = build_reps(L_sp, coords, a.n_pcs, a.k_spatial, a.seed)
    for n, E in reps.items():
        print(f"       {n:12s} → {E.shape}")

    print(f"Step 3 · 同折交叉验证比较两种表征({a.folds}-fold, kNN)")
    tasks = {"niche": niche, "cell type": ct_sp}
    res, cms_niche, tidy = {}, {}, []
    for tname, y in tasks.items():
        res[tname] = {}
        for rname, E in reps.items():
            f1s, cm, labs = cv_scores(E, y, a.folds, a.seed)
            res[tname][rname] = f1s
            tidy += [{"task": tname, "representation": rname, "fold": i,
                      "macro_f1": float(v)} for i, v in enumerate(f1s)]
            if tname == "niche":
                cms_niche[rname] = (cm, labs)
            print(f"       {tname:10s} · {rname:12s} macro-F1 = "
                  f"{f1s.mean():.3f} ± {f1s.std():.3f}")
        d = res[tname]["niche-aware"].mean() - res[tname]["intrinsic"].mean()
        print(f"       {tname:10s} · delta(niche-aware - intrinsic) = {d:+.3f}")

    print("Step 4 · 解离参考 → 空间 query 的细胞类型标签迁移(未校正地板)")
    tr = cross_modality_transfer(L_ref, ct_ref, L_sp, ct_sp, a.n_pcs, a.seed)
    print(f"       accuracy = {tr['accuracy']:.3f} · macro-F1 = {tr['macro_f1']:.3f}")

    print("Step 5 · 出图")
    for d in (outdir, assets):
        fig_tissue_and_embedding(coords, niche, ct_sp, reps, d)
        fig_slopegraph(res, d)
        fig_confusion(cms_niche, d)

    print("Step 6 · 落盘结果")
    pd.DataFrame(tidy).to_csv(outdir / "569_cv_macro_f1_per_fold.csv", index=False)
    summary = {
        "seed": a.seed, "n_pcs": a.n_pcs, "k_spatial": a.k_spatial, "folds": a.folds,
        "baseline_cv_macro_f1": {t: {r: {"mean": float(v.mean()), "sd": float(v.std())}
                                     for r, v in d.items()} for t, d in res.items()},
        "delta_niche_aware_minus_intrinsic": {
            t: float(d["niche-aware"].mean() - d["intrinsic"].mean()) for t, d in res.items()},
        "cross_modality_label_transfer": tr,
    }
    nf = run_nicheformer(a.checkpoint) if a.run_nf else \
        {"status": "not requested", "hint": "加 --run-nicheformer --checkpoint <path> 走官方路径"}
    summary["nicheformer_path"] = nf
    # 依赖版本快照(铁律6:可复现)—— save_session 的 Python 等价物
    import platform, sklearn, matplotlib
    summary["session"] = {"python": platform.python_version(), "numpy": np.__version__,
                          "pandas": pd.__version__, "scikit-learn": sklearn.__version__,
                          "matplotlib": matplotlib.__version__}
    print(f"       nicheformer path: {nf.get('status')} · {nf.get('reason', nf.get('hint', ''))}")
    (outdir / "569_summary.json").write_text(
        json.dumps(summary, indent=1, ensure_ascii=False), encoding="utf-8")
    print(f"[569] done → {outdir}")


if __name__ == "__main__":
    main()
