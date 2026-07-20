# -*- coding: utf-8 -*-
"""587 · RegFormer — GRN 重建评测台(Mamba 单细胞基础模型 vs 朴素基因嵌入基线)

RegFormer 的 GRN 重建下游任务本质上是:**基因嵌入 → 余弦相似度 → 只以 TF 为源的 top-k 有向边
→ 谱聚类分模块 → 富集**。本模块把这条「嵌入 → 图 → 评测」的下游链路完整地本地跑通,
并用**朴素基因嵌入(共表达谱 / PCA)**作为基线;RegFormer 的基础模型嵌入只是把第一步换掉。

因此本模块的定位是**评测台**,而不是 RegFormer 的封装:
  · 无需装 RegFormer 即可跑完(基线路径,CPU 分钟级);
  · 装了 RegFormer 并跑完官方 GRN 任务后,把它输出的 `gene_embedding.npy` + `gene_names.csv`
    用 `--embedding/--gene-names` 喂进来,即可在**完全相同的图构建与评测口径**下与基线对比。
  · `--run-regformer` 只做守卫式探测(仓库/依赖/GPU),不模拟、不伪造模型输出。

图构建规则逐条对齐上游源码 downstream_task/regformer_grn.py 的 `GrnTaskMamba.construct_grn`
(余弦相似度、自相似置 -1、阈值 >0.2、每个 TF 取 top-k=20、只有 TF 可作源节点),
谱聚类对齐 `evaluate_grn`(SpectralClustering, affinity='precomputed', random_state=42,
assign_labels='kmeans', k=5..30 步长 5)。

论文: Hu L, Qin H, Zhang Y, et al. RegFormer: a single-cell foundation model powered by
gene regulatory hierarchies. Nat Commun. 2026;17(1). doi:10.1038/s41467-026-72198-x · PMID 42086551
仓库: https://github.com/BGIResearch/RegFormer (API 读自 master 分支源码,见 README)
"""
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
for _up in [HERE, *HERE.parents]:
    if (_up / "_framework" / "pubstyle.py").exists():
        sys.path.insert(0, str(_up / "_framework"))
        break
try:
    from pubstyle import (set_pub_style, save_fig, pal,          # noqa: E402
                          NATURE_W1, NATURE_W2, CMAP_CONT, CMAP_DIVERGE)
except Exception:                       # 框架缺失时最小降级,不影响分析
    def set_pub_style(*a, **k): pass
    def save_fig(fig, f, dpi=300):
        Path(f).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(str(f) + ".pdf"); fig.savefig(str(f) + ".png", dpi=dpi)
    def pal(n=None, name="npg"):
        base = ["#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F"]
        return base[:n] if n else base
    NATURE_W1, NATURE_W2, CMAP_CONT, CMAP_DIVERGE = 3.5, 7.0, "viridis", "RdBu_r"

SEED = 42
EX = HERE / "example_data"
RESULTS = HERE / "results"
ASSETS = HERE / "assets"


# =============================================================================
# 1. 数据读取
# =============================================================================
def load_inputs(counts_csv: Path, tf_file: Path, truth_edges: Path | None,
                truth_modules: Path | None):
    """读取 cells×genes 计数表、TF 名单、(可选)合成金标准。"""
    X = pd.read_csv(counts_csv, comment="#", index_col=0)
    # 上游 tf_file 为无表头 tsv,取第 1 列(见 construct_grn: read_csv(header=None, sep='\t')[0])
    tfs = pd.read_csv(tf_file, header=None, sep="\t", comment="#")[0].astype(str).tolist()
    tfs = [t for t in tfs if t in X.columns]
    truth = pd.read_csv(truth_edges, comment="#") if truth_edges and Path(truth_edges).exists() else None
    mods = pd.read_csv(truth_modules, comment="#") if truth_modules and Path(truth_modules).exists() else None
    return X, tfs, truth, mods


def normalize(X: pd.DataFrame) -> np.ndarray:
    """标准 CPM-like 归一化 + log1p(cells×genes)。"""
    lib = X.values.sum(axis=1, keepdims=True)
    lib[lib == 0] = 1.0
    return np.log1p(X.values / lib * 1e4)


# =============================================================================
# 2. 基因嵌入:朴素基线(这就是 RegFormer 必须打败的地板)
# =============================================================================
def embed_coexpression(Xn: np.ndarray) -> np.ndarray:
    """最朴素的地板:基因嵌入 = 它自己的跨细胞表达谱(中心化)。
    在这条谱上做余弦相似度 ≈ 皮尔逊共表达,几十年来的 GRN 默认做法。"""
    Z = Xn - Xn.mean(axis=0, keepdims=True)
    return Z.T                                   # genes × cells


def embed_pca(Xn: np.ndarray, n_comp: int = 20) -> np.ndarray:
    """第二个地板:对 基因×细胞 矩阵做 PCA 得到低维基因嵌入。
    对应领域内的红线——基础模型的 zero-shot 嵌入常常打不过 PCA
    (Kedzierska et al., Genome Biol 2025 的结论),所以这里必须报出来。"""
    from sklearn.decomposition import PCA
    G = (Xn - Xn.mean(axis=0, keepdims=True)).T
    k = int(min(n_comp, G.shape[0] - 1, G.shape[1] - 1))
    return PCA(n_components=k, random_state=SEED).fit_transform(G)


def load_external_embedding(npy: Path, names_csv: Path, genes: list[str]):
    """载入外部基因嵌入(可直接喂 RegFormer GRN 任务输出的
    `gene_embedding.npy` + `gene_names.csv`,文件名读自上游 regformer_grn.py)。"""
    E = np.load(npy)
    nm = pd.read_csv(names_csv)
    nm = nm.iloc[:, 0].astype(str).tolist()
    if len(nm) != E.shape[0]:
        raise ValueError(f"嵌入行数 {E.shape[0]} 与基因名数 {len(nm)} 不一致")
    idx = {g: i for i, g in enumerate(nm)}
    keep = [g for g in genes if g in idx]
    if len(keep) < 10:
        raise ValueError("外部嵌入与表达矩阵基因名重叠不足 10 个,无法评测")
    return E[[idx[g] for g in keep]], keep


# =============================================================================
# 3. 图构建 —— 逐条对齐上游 GrnTaskMamba.construct_grn
# =============================================================================
def construct_grn(emb: np.ndarray, genes: list[str], tfs: list[str],
                  top_k: int = 20, thr: float = 0.2):
    """余弦相似度 → 只以 TF 为源 → 自相似置 -1 → 相似度 > thr → 每 TF 取 top_k。

    规则来自上游 downstream_task/regformer_grn.py::GrnTaskMamba.construct_grn(L202-269)
    (阈值 0.2 见 L249、top_k=20 见 L202 签名默认值,此处保持一致以确保口径可比)。

    偏差说明(诚实标注):上游 L237 把 `G.add_nodes_from(genes)` 注释掉了,孤立基因不进图;
    本模块保留全部基因作节点,使 n_nodes 恒等于基因数、谱聚类覆盖全基因(否则模块 ARI
    只在有边基因上定义,不同嵌入之间节点集不同、无法横向比)。边集规则与上游完全一致。
    """
    import networkx as nx
    from sklearn.metrics.pairwise import cosine_similarity

    S = cosine_similarity(emb)
    tfset = set(tfs)
    G = nx.DiGraph()
    G.add_nodes_from(genes)
    for i, g in enumerate(genes):
        if g not in tfset:
            continue
        s = S[i].copy()
        s[i] = -1
        valid = np.where(s > thr)[0]
        if len(valid) > top_k:
            valid = valid[np.argsort(s[valid])[-top_k:]]
        for j in valid:
            if not G.has_edge(g, genes[j]):
                G.add_edge(g, genes[j], weight=float(S[i, j]))
    return G, S


def spectral_modules(G, k_list=range(5, 31, 5)):
    """谱聚类分模块,对齐上游 evaluate_grn 的参数。

    偏差说明(诚实标注):上游把有向邻接直接传给 affinity='precomputed',
    而 sklearn 的谱聚类要求对称亲和矩阵,这里对称化为 (A+Aᵀ)/2 后再聚类。
    """
    import networkx as nx
    from sklearn.cluster import SpectralClustering
    nodes = list(G.nodes())
    A = nx.to_numpy_array(G, nodelist=nodes, weight="weight")
    A = np.abs((A + A.T) / 2.0)
    out = {}
    for k in k_list:
        if k >= len(nodes):
            continue
        try:
            lab = SpectralClustering(n_clusters=k, affinity="precomputed",
                                     random_state=SEED, assign_labels="kmeans").fit_predict(A)
            out[k] = pd.Series(lab, index=nodes)
        except Exception as e:
            print(f"       spectral k={k} failed: {type(e).__name__}")
    return out


# =============================================================================
# 4. 评测(仅在提供金标准时;真实数据上自动跳过)
# =============================================================================
def score_edges(S: np.ndarray, genes: list[str], tfs: list[str], G, truth: pd.DataFrame):
    """对 TF→gene 全部候选对按相似度排序算 AUPRC,并对实际出的边算 P/R/F1。"""
    from sklearn.metrics import average_precision_score, precision_recall_curve
    gi = {g: i for i, g in enumerate(genes)}
    true_set = {(r.source, r.target) for r in truth.itertuples()
                if r.source in gi and r.target in gi}
    pairs, y, sc = [], [], []
    for t in tfs:
        i = gi[t]
        for g in genes:
            if g == t:
                continue
            pairs.append((t, g)); y.append(int((t, g) in true_set)); sc.append(S[i, gi[g]])
    y = np.asarray(y); sc = np.asarray(sc)
    pred = set(G.edges())
    tp = len(pred & true_set)
    prec = tp / len(pred) if pred else 0.0
    rec = tp / len(true_set) if true_set else 0.0
    f1 = 2 * prec * rec / (prec + rec) if prec + rec else 0.0
    p_curve, r_curve, _ = precision_recall_curve(y, sc)
    return {
        "auprc": float(average_precision_score(y, sc)),
        "prevalence": float(y.mean()),
        "precision": float(prec), "recall": float(rec), "f1": float(f1),
        "n_pred_edges": len(pred), "n_true_edges": len(true_set), "n_tp": tp,
    }, (p_curve, r_curve)


def score_modules(memb: pd.Series, mods: pd.DataFrame) -> float:
    """谱聚类模块 vs 真实调控程序的 ARI(只在有金标准基因上算)。"""
    from sklearn.metrics import adjusted_rand_score
    m = mods.set_index("gene")["module"]
    common = [g for g in memb.index if g in m.index]
    return float(adjusted_rand_score(m.loc[common].values, memb.loc[common].values))


def per_tf_recall(G, truth: pd.DataFrame) -> pd.DataFrame:
    """每个 TF 召回了多少条真实靶基因边(用于 lollipop 图)。"""
    pred = {}
    for s, t in G.edges():
        pred.setdefault(s, set()).add(t)
    rows = []
    for tf, grp in truth.groupby("source"):
        tset = set(grp["target"])
        hit = len(pred.get(tf, set()) & tset)
        rows.append({"tf": tf, "n_true": len(tset), "n_hit": hit,
                     "recall": hit / len(tset) if tset else np.nan})
    return pd.DataFrame(rows).sort_values("recall", ascending=False)


# =============================================================================
# 5. 出图(全部非条形图:PR 曲线 / slopegraph / heatmap / lollipop)
# =============================================================================
def fig_pr(curves: dict, prevalence: float, out: Path):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(NATURE_W1, NATURE_W1 * 0.85))
    cols = pal(max(len(curves), 3), "npg")
    for (name, (p, r)), c in zip(curves.items(), cols):
        ax.plot(r, p, lw=1.6, color=c, label=name)
    ax.axhline(prevalence, ls="--", lw=1.0, color="#999999")
    ax.text(0.98, prevalence, "random", ha="right", va="bottom", fontsize=7, color="#666666")
    ax.set_xlabel("Recall"); ax.set_ylabel("Precision")
    ax.set_title("TF-target edge recovery")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.legend(loc="upper right", fontsize=7)
    save_fig(fig, out); plt.close(fig)


def fig_slope(metrics: dict, out: Path):
    """slopegraph:每个指标一条竖轴,嵌入方法之间连线对比(替代条形图)。"""
    import matplotlib.pyplot as plt
    keys = ["auprc", "precision", "recall", "f1", "module_ari"]
    labels = ["AUPRC", "Precision", "Recall", "F1", "Module ARI"]
    names = list(metrics.keys())
    cols = pal(max(len(names), 3), "npg")
    fig, ax = plt.subplots(figsize=(NATURE_W2 * 0.62, NATURE_W1 * 0.95))
    x = np.arange(len(keys))
    ymat = np.array([[metrics[nm].get(k, np.nan) for k in keys] for nm in names], dtype=float)
    for si, (nm, c) in enumerate(zip(names, cols)):
        y = ymat[si]
        ax.plot(x, y, "-o", color=c, lw=1.5, ms=5, label=nm)
        for xi, yi in zip(x, y):
            if np.isnan(yi):
                continue
            # 标签避让:同一指标上,取值高的标在上、低的标在下,避免数字互相压字
            col = ymat[:, xi]
            above = yi >= np.nanmax(col) - 1e-12
            ax.annotate(f"{yi:.2f}", (xi, yi), textcoords="offset points",
                        xytext=(0, 8 if above else -13), ha="center",
                        fontsize=6.5, color=c)
    for xi in x:
        ax.axvline(xi, color="#EDEDED", lw=0.8, zorder=0)
    ax.set_xticks(x); ax.set_xticklabels(labels)
    ax.set_ylabel("Score"); ax.set_ylim(-0.03, 1.05)
    ax.set_title("Gene-embedding comparison on the identical GRN harness")
    ax.legend(loc="lower left", fontsize=7, ncol=len(names))
    save_fig(fig, out); plt.close(fig)


def fig_heatmap(S: np.ndarray, genes: list[str], mods: pd.DataFrame, out: Path,
                n_show: int = 60):
    """基因-基因余弦相似度热图,按真实调控程序排序(看嵌入有没有抓到块结构)。"""
    import matplotlib.pyplot as plt
    m = mods.set_index("gene")["module"]
    # 按真实调控程序排序;未连基因(module=-1)排到最后,便于对比「有结构 vs 无结构」
    order = sorted([g for g in genes if g in m.index],
                   key=lambda g: (m[g] if m[g] >= 0 else 10**6, g))[:n_show]
    idx = [genes.index(g) for g in order]
    M = S[np.ix_(idx, idx)]
    fig, ax = plt.subplots(figsize=(NATURE_W1 * 1.15, NATURE_W1 * 1.05))
    im = ax.imshow(M, cmap=CMAP_DIVERGE, vmin=-1, vmax=1, interpolation="nearest")
    lab = m.loc[order].values
    bounds = np.where(np.diff(lab) != 0)[0] + 0.5
    for b in bounds:
        ax.axhline(b, color="black", lw=0.6); ax.axvline(b, color="black", lw=0.6)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("Genes (ordered by true program)"); ax.set_ylabel("Genes")
    ax.set_title("Cosine similarity of gene embeddings")
    fig.colorbar(im, ax=ax, shrink=0.75, label="cosine similarity")
    save_fig(fig, out); plt.close(fig)


def fig_lollipop(df: pd.DataFrame, out: Path):
    import matplotlib.pyplot as plt
    d = df.sort_values("recall")
    y = np.arange(len(d))
    fig, ax = plt.subplots(figsize=(NATURE_W1, max(2.0, 0.22 * len(d) + 1.0)))
    ax.hlines(y, 0, d["recall"], color="#C9C9C9", lw=1.2)
    ax.scatter(d["recall"], y, s=34, color=pal(1, "npg")[0], zorder=3)
    ax.set_yticks(y); ax.set_yticklabels(d["tf"])
    ax.set_xlabel("Recall of true targets"); ax.set_xlim(-0.02, 1.02)
    ax.set_title("Per-TF target recovery (baseline)")
    save_fig(fig, out); plt.close(fig)


# =============================================================================
# 6. RegFormer 守卫路径 —— 只探测,不伪造
# =============================================================================
def probe_regformer(repo_dir: str | None):
    """守卫式探测:仓库/依赖/GPU 是否齐备。任何一项缺失就诚实退出并打印真实命令。

    这里给出的入口点全部读自上游 master 分支源码,未做任何推测:
      · CLI    : python downstream_task/regformer_grn.py --config_file grn.toml
      · 类     : downstream_task/regformer_grn.py::GrnTaskMamba(config_file)
      · 方法   : .get_gene_expression_embedding() / .construct_grn(embeddings, top_k=20)
                 / .evaluate_grn(g_nx, gene_names) / .run_grn_analysis()
      · 产物   : <save_dir>/<run_name>/gene_embedding.npy, gene_names.csv, edges.csv,
                 grn_enrichment.csv
    """
    info = {
        "cli": "python downstream_task/regformer_grn.py --config_file grn.toml",
        "entry_class": "downstream_task/regformer_grn.py::GrnTaskMamba",
        "entry_methods": ["get_gene_expression_embedding()", "construct_grn(embeddings, top_k=20)",
                          "evaluate_grn(g_nx, gene_names)", "run_grn_analysis()"],
        "artifacts": ["gene_embedding.npy", "gene_names.csv", "edges.csv", "grn_enrichment.csv"],
        "config_keys_required": ["data_path", "vocab_file", "graph_path", "tf_file", "gene_sets",
                                 "load_model", "save_dir", "run_name", "model_name",
                                 "graph_sort", "max_seq_len", "n_bins", "device"],
        # 下面这段照抄上游 README「Installation」原文;注意上游写的 `cd RegFormer-Official`
        # 与实际 clone 出来的目录名 RegFormer 不一致(上游笔误),按实际目录名 cd。
        "install": ("git clone https://github.com/BGIResearch/RegFormer && cd RegFormer && "
                    "conda create -n regformer python=3.9 && pip install -r requirements.txt"),
    }
    if repo_dir:
        p = Path(repo_dir)
        if not (p / "regformer").is_dir():
            return {**info, "status": "skipped", "reason": f"{p} 下没有 regformer/ 包目录"}
        sys.path.insert(0, str(p))
    try:
        import regformer                                    # noqa: F401
    except Exception as e:
        # 区分两种失败:根本没拿到仓库 vs 拿到了但 CUDA 扩展没装
        if isinstance(e, ModuleNotFoundError) and getattr(e, "name", "") == "regformer":
            tail = "RegFormer 不在 PyPI(pypi.org 查无此包),须 clone 仓库后用 --repo-dir 指到仓库根目录"
        else:
            tail = ("仓库已找到,但其依赖未就绪(如 causal_conv1d_cuda / mamba_ssm 的 CUDA 扩展);"
                    "按上游 requirements.txt 在带 CUDA 的环境里安装后重试")
        return {**info, "status": "skipped",
                "reason": f"import regformer 失败({type(e).__name__}: {e});{tail}"}
    try:
        import mamba_ssm                                    # noqa: F401
        mamba = True
    except Exception:
        mamba = False
    try:
        import torch
        gpu = bool(torch.cuda.is_available())
    except Exception:
        gpu = False
    if not (mamba and gpu):
        return {**info, "status": "skipped",
                "reason": f"mamba_ssm={mamba}, cuda={gpu};上游 requirements.txt 钉的是 CUDA 构建"
                          f"(torch==2.0.0+cuda11.6 / dgl==1.1.2+cu118 / causal_conv1d / mamba_ssm),"
                          f"且 Docs/configs/grn_10k.toml 里 device='cuda',故按需 GPU 处理"}
    return {**info, "status": "ready",
            "next": "按 Docs/configs/grn_10k.toml 填好 data_path/load_model/vocab_file/graph_path/tf_file 后跑 CLI,"
                    "再把 gene_embedding.npy + gene_names.csv 用 --embedding/--gene-names 喂回本模块对比"}


# =============================================================================
# 7. 主流程
# =============================================================================
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--counts", default=str(EX / "expression_counts.csv"), help="cells×genes 计数 CSV")
    ap.add_argument("--tf-file", default=str(EX / "tf_list.txt"), help="TF 名单(无表头,第1列)")
    ap.add_argument("--truth-edges", default=str(EX / "ground_truth_edges.csv"), help="金标准边(可选)")
    ap.add_argument("--truth-modules", default=str(EX / "ground_truth_modules.csv"), help="金标准模块(可选)")
    ap.add_argument("--embedding", help="外部基因嵌入 .npy(如 RegFormer 的 gene_embedding.npy)")
    ap.add_argument("--gene-names", help="与 --embedding 配套的 gene_names.csv")
    ap.add_argument("--top-k", type=int, default=20, help="每个 TF 保留的边数(上游默认 20)")
    ap.add_argument("--sim-thr", type=float, default=0.2, help="相似度阈值(上游默认 0.2)")
    ap.add_argument("--n-pcs", type=int, default=20, help="PCA 基线的主成分数")
    ap.add_argument("--run-regformer", action="store_true", help="探测 RegFormer 环境(不伪造输出)")
    ap.add_argument("--repo-dir", help="RegFormer 仓库根目录(clone 后的路径)")
    ap.add_argument("--outdir", default=str(RESULTS))
    a = ap.parse_args()

    np.random.seed(SEED)
    set_pub_style()
    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    ASSETS.mkdir(parents=True, exist_ok=True)

    print("[587] Step 1 载入数据")
    X, tfs, truth, mods = load_inputs(Path(a.counts), Path(a.tf_file),
                                      a.truth_edges, a.truth_modules)
    genes = list(X.columns)
    print(f"       {X.shape[0]} cells × {X.shape[1]} genes, {len(tfs)} TFs, "
          f"truth={'yes' if truth is not None else 'no'}")
    Xn = normalize(X)

    print("[587] Step 2 基因嵌入(基线 + 可选外部嵌入)")
    embs: dict[str, tuple[np.ndarray, list[str]]] = {
        "Co-expression (baseline)": (embed_coexpression(Xn), genes),
        f"PCA-{a.n_pcs} (baseline)": (embed_pca(Xn, a.n_pcs), genes),
    }
    if a.embedding and a.gene_names:
        E, keep = load_external_embedding(Path(a.embedding), Path(a.gene_names), genes)
        embs["External embedding"] = (E, keep)
        print(f"       外部嵌入载入: {E.shape[0]} genes × {E.shape[1]} dims")
    elif a.embedding or a.gene_names:
        sys.exit("--embedding 与 --gene-names 必须同时给出")

    print("[587] Step 3 构图 + 谱聚类 + 评测(口径对齐上游 construct_grn/evaluate_grn)")
    metrics, curves, graphs = {}, {}, {}
    for name, (E, gset) in embs.items():
        G, S = construct_grn(E, gset, [t for t in tfs if t in gset], a.top_k, a.sim_thr)
        graphs[name] = (G, S, gset)
        rec = {"n_nodes": G.number_of_nodes(), "n_edges": G.number_of_edges()}
        if truth is not None:
            sc, pr = score_edges(S, gset, [t for t in tfs if t in gset], G, truth)
            rec.update(sc); curves[name] = pr
        if mods is not None:
            memb = spectral_modules(G)
            best = max(((k, score_modules(v, mods)) for k, v in memb.items()),
                       key=lambda kv: kv[1], default=(None, np.nan))
            rec["module_ari"], rec["module_best_k"] = best[1], best[0]
            if best[0] is not None:
                memb[best[0]].rename("module").to_frame().to_csv(
                    out / f"modules_{name.split()[0].lower()}.csv")
        metrics[name] = rec
        print(f"       {name}: " + ", ".join(
            f"{k}={v:.3f}" if isinstance(v, float) else f"{k}={v}" for k, v in rec.items()))

    # 边列表按上游 edges.csv 格式(nx.to_pandas_edgelist)落盘
    import networkx as nx
    for name, (G, _, _) in graphs.items():
        nx.to_pandas_edgelist(G).to_csv(
            out / f"edges_{name.split()[0].lower()}.csv", index=False)

    print("[587] Step 4 出图")
    made = []
    if curves:
        prev = next(iter(metrics.values())).get("prevalence", 0.0)
        for dst in (ASSETS, out):
            fig_pr(curves, prev, dst / "fig1_pr_curve")
        made.append("fig1_pr_curve.png")
    if truth is not None or mods is not None:
        for dst in (ASSETS, out):
            fig_slope(metrics, dst / "fig2_metric_slopegraph")
        made.append("fig2_metric_slopegraph.png")
    if mods is not None:
        base_name = "Co-expression (baseline)"
        G, S, gset = graphs[base_name]
        for dst in (ASSETS, out):
            fig_heatmap(S, gset, mods, dst / "fig3_similarity_heatmap")
        made.append("fig3_similarity_heatmap.png")
    if truth is not None:
        ptf = per_tf_recall(graphs["Co-expression (baseline)"][0], truth)
        ptf.to_csv(out / "per_tf_recall.csv", index=False)
        for dst in (ASSETS, out):
            fig_lollipop(ptf, dst / "fig4_tf_recall_lollipop")
        made.append("fig4_tf_recall_lollipop.png")
    print("       " + ", ".join(made))

    print("[587] Step 5 RegFormer 路径")
    if a.run_regformer:
        rf = probe_regformer(a.repo_dir)
        print(f"       status: {rf['status']}")
        print(f"       {rf.get('reason', rf.get('next', ''))}")
    else:
        rf = {"status": "not requested",
              "hint": "加 --run-regformer [--repo-dir <RegFormer clone>] 探测环境"}
        print("       未请求(--run-regformer);仅跑基线对比")

    # 依赖版本快照(铁律6 可复现):关键数字与包版本一起落盘
    import platform, sklearn, networkx as _nx, matplotlib as _mpl
    session = {"python": platform.python_version(), "numpy": np.__version__,
               "pandas": pd.__version__, "scikit-learn": sklearn.__version__,
               "networkx": _nx.__version__, "matplotlib": _mpl.__version__}

    summary = {"seed": SEED, "n_cells": int(X.shape[0]), "n_genes": int(X.shape[1]),
               "n_tf": len(tfs), "top_k": a.top_k, "sim_thr": a.sim_thr,
               "metrics": metrics, "regformer": rf, "figures": made,
               "session_info": session}
    (out / "587_summary.json").write_text(
        json.dumps(summary, indent=1, default=str), encoding="utf-8")
    pd.DataFrame(metrics).T.to_csv(out / "metrics_table.csv")
    print(f"[587] done → {out / '587_summary.json'}")


if __name__ == "__main__":
    main()
