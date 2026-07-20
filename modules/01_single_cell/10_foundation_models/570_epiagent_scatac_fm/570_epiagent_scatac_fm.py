"""570 · EpiAgent — foundation model for single-cell chromatin accessibility (scATAC).

本模块两条路径:
  1) 基线(永远可跑, CPU): TF-IDF + SVD 得到 LSI 嵌入 → 聚类 / 细胞类型分类 /
     低秩重构填补。这是 scATAC 领域的朴素标准做法(Signac/ArchR 的 LSI),
     用来给 EpiAgent 的三项声明能力(embedding / cell-type prediction / imputation)
     各配一个可量化的对照。基础模型若不比它好, 就不该被报道。
  2) EpiAgent 路径(--run-epiagent, 需装包 + CUDA GPU + flash-attn): 守卫式封装。
     只做「能否运行」的探测并打印上游真实 API 与官方教程指引, 不臆造调用序列。

上游 API 来源(逐文件实读, 见 README「API 溯源」表), 例如:
  epiagent.preprocessing.construct_cell_by_ccre_matrix / global_TFIDF
  epiagent.tokenization.tokenization
  epiagent.model.EpiAgent / EpiAgent_supervised / EpiAgent_BC / EpiAgent_PT
  epiagent.dataset.CellDataset / collate_fn ...
  epiagent.inference.infer_cell_embeddings / infer_cell_types / infer_reconstructed_signals

论文: Chen X, et al. EpiAgent: foundation model for single-cell epigenomics.
      Nat Methods 2025. doi:10.1038/s41592-025-02822-z · PMID 40999099 (已核实)
仓库: https://github.com/xy-chen16/EpiAgent
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import warnings

warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")
ASSETS = os.path.join(HERE, "assets")
EXAMPLE = os.path.join(HERE, "example_data")

# 框架统一出图样式
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..", "..", "_framework")))
from pubstyle import (  # noqa: E402
    CMAP_CONT, NATURE_W1, NATURE_W2, pal, save_fig, set_pub_style,
)

SEED = 0


# --------------------------------------------------------------------------- #
# 合成示例数据: 细胞 × cCRE 二值可及性矩阵                                        #
# --------------------------------------------------------------------------- #
def make_synthetic(n_cells=600, n_peaks=900, n_types=4, n_batches=2, seed=SEED):
    """生成小型合成 scATAC 计数矩阵(synthetic, for demo only)。

    结构: 每个细胞类型有一组特异开放的 peak(细胞类型信号) + 一批管家 peak(共开放)
          + 批次特异的整体捕获效率差异(技术批次效应) + 泊松/伯努利稀疏化。
    真值标签写进 obs, 供基线做有监督评估(注意: 评估时严格 train/test 划分, 防泄漏)。
    """
    import numpy as np
    import pandas as pd

    rng = np.random.default_rng(seed)
    cell_type = rng.integers(0, n_types, n_cells)
    batch = rng.integers(0, n_batches, n_cells)

    # peak 分区: 前 n_types*k 个为类型特异, 其余为管家/背景
    # 信噪比刻意压低(类型特异 0.16 vs 背景 0.04), 让基线不至于饱和到 ARI=1.0
    # ——饱和的对照没有区分力, 基础模型再好也看不出增益。
    k = 120
    prob = np.full((n_cells, n_peaks), 0.04)          # 背景可及性
    house = slice(n_types * k, n_types * k + 200)
    prob[:, house] = 0.25                              # 管家 peak 全细胞开放
    for t in range(n_types):
        idx = slice(t * k, (t + 1) * k)
        prob[cell_type == t, idx] = 0.16               # 类型特异 peak(弱信号)

    # 相邻细胞类型共享一半特异 peak → 制造易混淆的类型对(真实组织常见)
    for t in range(n_types - 1):
        shared = slice(t * k, t * k + k // 2)
        prob[cell_type == t + 1, shared] = 0.13

    # 批次效应: 整体捕获效率(测序深度)差异, 是纯技术项
    depth = np.where(batch == 0, 1.0, 0.62)[:, None]
    # 细胞层面的深度抖动(文库复杂度差异)
    depth = depth * rng.lognormal(0, 0.25, (n_cells, 1))
    prob = np.clip(prob * depth, 0, 0.95)

    X = (rng.random((n_cells, n_peaks)) < prob).astype("float32")

    obs = pd.DataFrame({
        "cell_type": [f"CT{t}" for t in cell_type],
        "batch": [f"B{b}" for b in batch],
    }, index=[f"Cell{i:04d}" for i in range(n_cells)])
    var = pd.DataFrame(index=[f"cCRE{j:05d}" for j in range(n_peaks)])
    return X, obs, var


def write_example(outdir=EXAMPLE):
    """把合成矩阵落盘为 example_data/(csv, 便于人眼检查)。"""
    import pandas as pd

    os.makedirs(outdir, exist_ok=True)
    X, obs, var = make_synthetic()
    mat = pd.DataFrame(X.astype(int), index=obs.index, columns=var.index)
    header = "# synthetic, for demo only -- NOT real scATAC data\n"
    p_mat = os.path.join(outdir, "cell_by_ccre_counts.csv")
    p_obs = os.path.join(outdir, "cell_metadata.csv")
    with open(p_mat, "w", newline="") as fh:
        fh.write(header)
        mat.to_csv(fh)
    with open(p_obs, "w", newline="") as fh:
        fh.write(header)
        obs.to_csv(fh)
    return p_mat, p_obs


def load_input(mat_csv, obs_csv):
    """读用户矩阵(行=细胞, 列=cCRE/peak)+ 元数据。"""
    import pandas as pd

    mat = pd.read_csv(mat_csv, index_col=0, comment="#")
    obs = pd.read_csv(obs_csv, index_col=0, comment="#")
    obs = obs.loc[mat.index]
    return mat.values.astype("float32"), obs, pd.DataFrame(index=mat.columns)


# --------------------------------------------------------------------------- #
# 基线: TF-IDF + SVD (LSI)                                                      #
# --------------------------------------------------------------------------- #
def tfidf(X):
    """scATAC 标准 TF-IDF。概念上对应上游 epiagent.preprocessing.global_TFIDF,
    但此处是本模块自实现的朴素版本, 不冒充上游函数。"""
    import numpy as np

    Xb = (X > 0).astype("float32")
    depth = Xb.sum(1, keepdims=True)
    depth[depth == 0] = 1.0
    tf = Xb / depth
    df = Xb.sum(0)
    idf = np.log(1.0 + Xb.shape[0] / (1.0 + df))
    return np.log1p(tf * idf * 1e4)


def lsi_embedding(X, n_comp=30, drop_first=True, seed=SEED):
    """TF-IDF → 截断 SVD。第 1 主成分通常与测序深度高度相关, 默认丢弃(领域惯例)。"""
    from sklearn.decomposition import TruncatedSVD

    Z = tfidf(X)
    svd = TruncatedSVD(n_components=n_comp, random_state=seed)
    emb = svd.fit_transform(Z)
    return emb[:, 1:] if drop_first else emb


def eval_clustering(emb, labels, seed=SEED):
    """无监督: KMeans 聚类 vs 真值细胞类型, 报 ARI / NMI。"""
    import numpy as np
    from sklearn.cluster import KMeans
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

    k = len(np.unique(labels))
    pred = KMeans(n_clusters=k, n_init=10, random_state=seed).fit_predict(emb)
    return {
        "ARI": float(adjusted_rand_score(labels, pred)),
        "NMI": float(normalized_mutual_info_score(labels, pred)),
    }, pred


def eval_celltype_clf(emb, labels, seed=SEED):
    """有监督: 嵌入 + 逻辑回归做细胞类型预测。

    防泄漏: 嵌入在全体细胞上无监督拟合(不看标签), 分类器只在 train 折拟合,
    指标只在 test 折报告。这是 EpiAgent「cell-type prediction」的朴素对照。
    """
    import numpy as np
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import accuracy_score, confusion_matrix, f1_score
    from sklearn.model_selection import train_test_split

    idx = np.arange(len(labels))
    tr, te = train_test_split(idx, test_size=0.3, random_state=seed, stratify=labels)
    clf = LogisticRegression(max_iter=2000, random_state=seed).fit(emb[tr], labels[tr])
    pred = clf.predict(emb[te])
    classes = sorted(np.unique(labels))
    cm = confusion_matrix(labels[te], pred, labels=classes)
    return {
        "accuracy": float(accuracy_score(labels[te], pred)),
        "macro_F1": float(f1_score(labels[te], pred, average="macro")),
    }, cm, classes


def eval_imputation(X, n_comp=20, mask_frac=0.15, seed=SEED):
    """填补: 随机遮蔽一部分「已观测为开放」的位点 + 等量闭合位点, 用低秩 SVD
    重构打分, 报 AUROC。这是 EpiAgent「imputation」的朴素对照。"""
    import numpy as np
    from sklearn.decomposition import TruncatedSVD
    from sklearn.metrics import roc_auc_score

    rng = np.random.default_rng(seed)
    Xb = (X > 0).astype("float32")
    pos = np.argwhere(Xb == 1)
    neg = np.argwhere(Xb == 0)
    n_mask = int(len(pos) * mask_frac)
    pos_sel = pos[rng.choice(len(pos), n_mask, replace=False)]
    neg_sel = neg[rng.choice(len(neg), n_mask, replace=False)]

    Xtrain = Xb.copy()
    Xtrain[pos_sel[:, 0], pos_sel[:, 1]] = 0          # 把开放位点打成 dropout

    Z = tfidf(Xtrain)
    svd = TruncatedSVD(n_components=n_comp, random_state=seed)
    rec = svd.inverse_transform(svd.fit_transform(Z))

    held = np.vstack([pos_sel, neg_sel])
    y = np.r_[np.ones(len(pos_sel)), np.zeros(len(neg_sel))]
    score = rec[held[:, 0], held[:, 1]]
    return {"imputation_AUROC": float(roc_auc_score(y, score)),
            "n_masked_each_class": int(n_mask)}


def eval_batch_mixing(emb, batch, seed=SEED):
    """批次混合度: kNN 中异批次邻居比例(越高越混合)。同时给一个「完美混合」期望值。"""
    import numpy as np
    from sklearn.neighbors import NearestNeighbors

    nn = NearestNeighbors(n_neighbors=21).fit(emb)
    _, ind = nn.kneighbors(emb)
    ind = ind[:, 1:]
    same = (batch[ind] == batch[:, None]).mean()
    _, cnt = np.unique(batch, return_counts=True)
    p = cnt / cnt.sum()
    return {"kNN_same_batch_frac": float(same),
            "expected_if_well_mixed": float((p ** 2).sum())}


# --------------------------------------------------------------------------- #
# EpiAgent 路径(守卫式)                                                         #
# --------------------------------------------------------------------------- #
def run_epiagent(probe_only=True):
    """守卫式封装: 只探测环境并回报上游真实 API, 不构造未经官方教程验证的调用序列。

    真实签名(逐行核对上游源码 v0.0.3, 行号见 README「API 溯源」表):
      preprocessing.construct_cell_by_ccre_matrix(intersect_file, ccre_bed_path)  # preprocessing.py:8
      preprocessing.global_TFIDF(adata, cCRE_document_frequency)                  # preprocessing.py:52
      tokenization.tokenization(adata, num_cCREs=1355445)                         # tokenization.py:6
      dataset.CellDataset(cell_sentences, max_length=8192, is_random=True)        # dataset.py:6
      dataset.collate_fn(data)                                                    # dataset.py:75
      model.EpiAgent(vocab_size=1355449, num_layers=18, embedding_dim=512,        # model.py:8
                     num_attention_heads=8, max_rank_embeddings=8192, ...,
                     use_flash_attn=True)
      model.EpiAgent_supervised(..., num_classes=10, use_flash_attn=True)         # model.py:268
      inference.infer_cell_embeddings(model, device, dataloader)                  # inference.py:78
      inference.infer_cell_types(model, device, dataloader,                       # inference.py:117
                                 need_cell_embeddings=True)
      inference.infer_reconstructed_signals(model, device, dataloader,            # inference.py:7
                                            need_cell_embeddings=True,
                                            predicted_cCRE_indices=None)
    注意: 上游 epiagent/__init__.py 是空文件, 不 re-export 也不定义 __version__,
          必须显式 import 子模块。权重加载方式与端到端串联顺序以官方 demo 为准, 此处未固定。
    """
    info = {}
    try:
        import epiagent  # noqa: F401
    except ImportError as e:
        return {"status": "skipped",
                "reason": f"epiagent 未安装 ({getattr(e, 'name', e)})",
                "install": ("conda create -n EpiAgent python=3.11 && conda activate EpiAgent && "
                            "pip install torch==2.0.1 torchvision==0.15.2 torchaudio==2.0.2 && "
                            "pip uninstall -y ninja && pip install ninja && "
                            "pip install flash-attn==2.5.8 --no-build-isolation && "
                            "pip install epiagent"),
                "note": "默认 use_flash_attn=True → 需要 NVIDIA GPU + flash-attn; CPU 不可行"}
    try:
        import torch
        info["cuda"] = bool(torch.cuda.is_available())
    except ImportError:
        info["cuda"] = False
    try:
        import flash_attn  # noqa: F401
        info["flash_attn"] = True
    except ImportError:
        info["flash_attn"] = False

    if not info["cuda"]:
        info.update(status="skipped", reason="无 CUDA GPU; EpiAgent 为 18 层 Transformer, CPU 推理不现实")
        return info

    # __init__.py 为空 → 顶层 import 成功不代表子模块可用(如 flash_attn 缺失会在
    # import epiagent.model 时才炸)。逐个显式探测, 如实回报。
    submods = {}
    for name in ("preprocessing", "tokenization", "dataset", "model", "inference"):
        try:
            __import__(f"epiagent.{name}")
            submods[name] = "ok"
        except Exception as e:                      # noqa: BLE001
            submods[name] = f"{type(e).__name__}: {e}"
    info["submodules"] = submods

    info.update(
        status="ready" if all(v == "ok" for v in submods.values()) else "partial",
        version=getattr(__import__("epiagent"), "__version__",
                        "未定义(上游 __init__.py 为空)"),
        next=("按官方 demo notebook 串联: preprocessing.global_TFIDF → tokenization.tokenization "
              "→ dataset.CellDataset(+collate_fn) → 载入预训练权重到 model.EpiAgent / "
              "EpiAgent_supervised → inference.infer_cell_embeddings / infer_cell_types / "
              "infer_reconstructed_signals"),
        demos="https://github.com/xy-chen16/EpiAgent/tree/main/demo/",
        weights=("预训练权重(EpiAgent / EpiAgent-B / EpiAgent-NT)+ 示例文件 + cCRE_frequency.npy "
                 "由上游放在 Google Drive(见上游 README 链接), 非 GitHub release"),
        repo="https://github.com/xy-chen16/EpiAgent",
        probe_only=probe_only,
    )
    return info


# --------------------------------------------------------------------------- #
# 出图(全部非条形图)                                                            #
# --------------------------------------------------------------------------- #
def fig_embedding(emb, obs, outstem):
    """图1: LSI 嵌入散点, 左=按细胞类型着色, 右=按批次着色。"""
    import matplotlib.pyplot as plt
    import numpy as np

    fig, axes = plt.subplots(1, 2, figsize=(NATURE_W2, 3.1))
    for ax, key, title in zip(axes, ["cell_type", "batch"],
                              ["Colored by cell type", "Colored by batch"]):
        lev = sorted(obs[key].unique())
        cols = pal(len(lev), "npg" if key == "cell_type" else "okabe_ito")
        for lv, c in zip(lev, cols):
            m = (obs[key] == lv).values
            ax.scatter(emb[m, 0], emb[m, 1], s=6, c=c, alpha=0.75,
                       linewidths=0, label=lv)
        ax.set_xlabel("LSI-1"); ax.set_ylabel("LSI-2"); ax.set_title(title)
        ax.legend(markerscale=2.2, loc="best", handletextpad=0.2)
    fig.suptitle("Baseline scATAC embedding (TF-IDF + SVD)", y=1.03, fontweight="bold")
    fig.tight_layout()
    save_fig(fig, outstem)
    plt.close(fig)


def fig_metrics_lollipop(metrics, outstem):
    """图2: 基线各任务指标 lollipop(不用条形图), 附随机/朴素参照线。"""
    import matplotlib.pyplot as plt
    import numpy as np

    items = [(k, v, ref) for k, v, ref in metrics]
    items = items[::-1]
    names = [i[0] for i in items]
    vals = np.array([i[1] for i in items])
    refs = np.array([i[2] for i in items])
    y = np.arange(len(names))
    cols = pal(len(names), "npg")

    fig, ax = plt.subplots(figsize=(NATURE_W1 + 1.6, 0.46 * len(names) + 1.5))
    ax.hlines(y, refs, vals, color="#BBBBBB", lw=2, zorder=1)
    ax.scatter(refs, y, s=52, facecolors="white", edgecolors="#777777",
               lw=1.3, zorder=3, label="Naive reference")
    ax.scatter(vals, y, s=74, c=cols, zorder=4, label="LSI baseline")
    for yi, v in zip(y, vals):
        ax.text(v + 0.025, yi, f"{v:.3f}", va="center", fontsize=9)
    ax.set_yticks(y); ax.set_yticklabels(names)
    ax.set_xlim(0, 1.16); ax.set_xlabel("Score")
    ax.set_ylim(-0.6, len(names) - 0.4)
    # 图例放到坐标区上方横排, 避免压住最下面一行的数值标签
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=2,
              columnspacing=1.4, handletextpad=0.3)
    ax.set_title("Baseline performance vs naive reference", pad=34)
    fig.tight_layout()
    save_fig(fig, outstem)
    plt.close(fig)


def fig_confusion(cm, classes, outstem):
    """图3: 细胞类型预测混淆矩阵热图(行归一化)。"""
    import matplotlib.pyplot as plt
    import numpy as np

    cmn = cm / np.clip(cm.sum(1, keepdims=True), 1, None)
    fig, ax = plt.subplots(figsize=(NATURE_W1 + 0.6, NATURE_W1 + 0.1))
    im = ax.imshow(cmn, cmap=CMAP_CONT, vmin=0, vmax=1)
    ax.set_xticks(range(len(classes)), classes)
    ax.set_yticks(range(len(classes)), classes)
    ax.set_xlabel("Predicted"); ax.set_ylabel("True")
    ax.set_title("Cell-type prediction (held-out 30%)")
    for i in range(len(classes)):
        for j in range(len(classes)):
            ax.text(j, i, f"{cmn[i, j]:.2f}", ha="center", va="center",
                    fontsize=8, color="white" if cmn[i, j] < 0.6 else "black")
    fig.colorbar(im, ax=ax, shrink=0.8, label="Row-normalized fraction")
    fig.tight_layout()
    save_fig(fig, outstem)
    plt.close(fig)


def fig_depth_violin(X, obs, outstem):
    """图4: 每批次的 per-cell 开放位点数分布(violin + 抖动散点), 暴露技术批次效应。"""
    import matplotlib.pyplot as plt
    import numpy as np

    rng = np.random.default_rng(SEED)
    depth = (X > 0).sum(1)
    lev = sorted(obs["batch"].unique())
    data = [depth[(obs["batch"] == lv).values] for lv in lev]
    cols = pal(len(lev), "okabe_ito")

    fig, ax = plt.subplots(figsize=(NATURE_W1, 3.0))
    vp = ax.violinplot(data, showextrema=False, widths=0.8)
    for body, c in zip(vp["bodies"], cols):
        body.set_facecolor(c); body.set_alpha(0.35); body.set_edgecolor(c)
    for i, (d, c) in enumerate(zip(data, cols), start=1):
        ax.scatter(i + rng.normal(0, 0.045, len(d)), d, s=4, c=c, alpha=0.5, linewidths=0)
        ax.hlines(np.median(d), i - 0.28, i + 0.28, color="black", lw=1.6, zorder=5)
    ax.set_xticks(range(1, len(lev) + 1), lev)
    ax.set_xlabel("Batch"); ax.set_ylabel("Open cCREs per cell")
    ax.set_title("Per-cell coverage by batch")
    fig.tight_layout()
    save_fig(fig, outstem)
    plt.close(fig)


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description="570 EpiAgent scATAC 基础模型 + LSI 基线")
    ap.add_argument("--matrix", default=os.path.join(EXAMPLE, "cell_by_ccre_counts.csv"),
                    help="细胞×cCRE 计数矩阵 csv(行=细胞)")
    ap.add_argument("--meta", default=os.path.join(EXAMPLE, "cell_metadata.csv"),
                    help="细胞元数据 csv, 需含 cell_type / batch 列")
    ap.add_argument("--outdir", default=RESULTS)
    ap.add_argument("--n-comp", type=int, default=30, help="SVD 成分数")
    ap.add_argument("--run-epiagent", action="store_true",
                    help="尝试 EpiAgent 路径(需 pip install epiagent + CUDA + flash-attn)")
    ap.add_argument("--regen-example", action="store_true", help="重新生成 example_data")
    a = ap.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    os.makedirs(ASSETS, exist_ok=True)
    set_pub_style()

    if a.regen_example or not os.path.exists(a.matrix):
        print("[570] Step 0 · 生成合成示例数据 example_data/")
        write_example()

    print("[570] Step 1 · 读入 cell × cCRE 矩阵")
    import numpy as np
    X, obs, var = load_input(a.matrix, a.meta)
    print(f"       {X.shape[0]} cells × {X.shape[1]} cCREs, "
          f"sparsity={1 - (X > 0).mean():.3f}")

    print("[570] Step 2 · 基线嵌入 TF-IDF + SVD (LSI)")
    emb = lsi_embedding(X, n_comp=a.n_comp)

    labels = obs["cell_type"].values
    batch = obs["batch"].values

    print("[570] Step 3 · 基线评估(聚类 / 细胞类型预测 / 填补 / 批次混合)")
    clu, clu_pred = eval_clustering(emb, labels)
    clf, cm, classes = eval_celltype_clf(emb, labels)
    imp = eval_imputation(X)
    bmx = eval_batch_mixing(emb, batch)
    for d in (clu, clf, imp, bmx):
        for k, v in d.items():
            print(f"       {k}: {v}")

    # 朴素参照: ARI/NMI 随机划分≈0; 分类按最大类占比(多数类猜测); AUROC 随机=0.5
    maj = float(max((labels == c).mean() for c in classes))
    metrics = [
        ("Clustering ARI",        clu["ARI"],               0.0),
        ("Clustering NMI",        clu["NMI"],               0.0),
        ("Cell-type accuracy",    clf["accuracy"],          maj),
        ("Cell-type macro-F1",    clf["macro_F1"],          1.0 / len(classes)),
        ("Imputation AUROC",      imp["imputation_AUROC"],  0.5),
    ]

    print("[570] Step 4 · 出图")
    fig_embedding(emb, obs, os.path.join(a.outdir, "570_embedding_scatter"))
    fig_metrics_lollipop(metrics, os.path.join(a.outdir, "570_baseline_metrics_lollipop"))
    fig_confusion(cm, classes, os.path.join(a.outdir, "570_celltype_confusion"))
    fig_depth_violin(X, obs, os.path.join(a.outdir, "570_coverage_violin"))
    # 复制展示图到 assets/(README 引用)
    import shutil
    for stem in ("570_embedding_scatter", "570_baseline_metrics_lollipop",
                 "570_celltype_confusion", "570_coverage_violin"):
        shutil.copyfile(os.path.join(a.outdir, stem + ".png"),
                        os.path.join(ASSETS, stem + ".png"))

    print("[570] Step 5 · EpiAgent 路径")
    ep = run_epiagent() if a.run_epiagent else {
        "status": "not requested", "hint": "加 --run-epiagent 探测环境"}
    for k, v in ep.items():
        print(f"       {k}: {v}")

    import pandas as pd
    pd.DataFrame({"cell": obs.index, "cell_type": labels, "batch": batch,
                  "kmeans_cluster": clu_pred,
                  **{f"LSI{i+1}": emb[:, i] for i in range(min(5, emb.shape[1]))}}
                 ).to_csv(os.path.join(a.outdir, "570_baseline_embedding.csv"), index=False)
    # 依赖版本快照(铁律6 可复现: 关键数字须能追溯到具体环境)
    import platform
    import sklearn
    import matplotlib
    session = {"python": platform.python_version(), "numpy": np.__version__,
               "pandas": pd.__version__, "scikit-learn": sklearn.__version__,
               "matplotlib": matplotlib.__version__}
    summary = {"input": {"n_cells": int(X.shape[0]), "n_ccres": int(X.shape[1])},
               "baseline": {**clu, **clf, **imp, **bmx, "majority_class_frac": maj},
               "epiagent": ep, "seed": SEED, "session_info": session}
    with open(os.path.join(a.outdir, "570_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1, default=str)
    print(f"[570] done · results → {a.outdir}")


if __name__ == "__main__":
    main()
