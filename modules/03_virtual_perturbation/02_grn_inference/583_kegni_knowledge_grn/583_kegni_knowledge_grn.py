"""583 · KEGNI —— 知识图增强的基因调控网络推断(可跑基线 + 上游守卫式封装)。

上游方法
    KEGNI (Knowledge graph-Enhanced Gene regulatory Network Inference)
    Li P, Li L, Nan J, Chen J, Sun J, Cao Y. KEGNI: knowledge graph enhanced framework
    for gene regulatory network inference. Genome Biology 2025;26:294.
    PMID 40983951 · doi:10.1186/s13059-025-03780-7   (两者本地已用 E-utilities / doi.org 核实)
    仓库 https://github.com/Lipxiao/KEGNI  (无 pip 包,是命令行研究代码)

本模块做两件事
    (A) 基线套件 —— 只用本机已有依赖(numpy/scipy/sklearn/networkx/matplotlib)就能跑完:
        相关性 / 秩相关 / PCA 基因嵌入点积 / 纯知识先验 / 知识-表达秩融合,
        全部按 BEELINE 口径用 EPR + AUPRC + AUROC 对同一套 ground truth 打分。
        这是"知识嵌入到底加了多少分"的诚实地板 —— 任何 KEGNI 结果都必须跟它比。
    (B) KEGNI 真身 —— 守卫式封装:检查 torch/dgl 与本地仓库,凑齐就打印/执行经核实的
        真实命令行;缺依赖就优雅退出并打印真实安装与运行命令,绝不伪造 API。

KEGNI 的接口是命令行,不是 python 函数。以下参数逐一读自本地克隆的 utils/args.py 与 run_mESC.sh:
    python train.py --input <expr.csv> --data_path <kg.tsv> --genes <int>
                    --n_neighbors <int> --num_hidden <int> --num_heads <int>
                    --num_layers <int> --norm -1 --max_steps <int> --lambda_kge <float>
    python eval.py -p <pred.csv> -t <ground_truth.csv>
`--eval` 只对 BEELINE 命名的数据集有效(ground truth 路径在 trainer.py 里写死),
故默认不加;详见 kegni_guard 的文档字符串。
逐条 API 溯源(文件:行号)见 README「API 溯源」一节。
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
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..", "..", "_framework")))
from pubstyle import (CMAP_CONT, NATURE_W1, NATURE_W2, pal, save_fig,  # noqa: E402
                      set_pub_style)

import matplotlib.pyplot as plt  # noqa: E402
import networkx as nx  # noqa: E402
from scipy.stats import rankdata, spearmanr  # noqa: E402
from sklearn.decomposition import PCA  # noqa: E402
from sklearn.metrics import (average_precision_score, precision_recall_curve,  # noqa: E402
                             roc_auc_score)

SEED = 2024

# 上游经核实的真实命令行默认值(来自 configs.yml 的 mESC 档 + run_mESC.sh)
KEGNI_DEFAULTS = dict(n_neighbors=30, num_hidden=512, num_heads=4, num_layers=4,
                      norm=-1, max_steps=1000, lambda_kge=1.0)


# --------------------------------------------------------------------------- #
# 一、数据读取(格式严格对齐上游)
# --------------------------------------------------------------------------- #
def load_inputs(expr_file: str, kg_file: str, gt_file: str):
    """读表达矩阵 / 知识图三元组 / ground truth。基因名一律大写(上游同样这么做)。"""
    expr = pd.read_csv(expr_file, index_col=0, header=0)
    expr.index = expr.index.str.upper()

    kg = pd.read_csv(kg_file, sep="\t", header=None)      # head / relation / tail
    kg.columns = ["head", "relation", "tail"]
    kg["head"] = kg["head"].str.upper()
    kg["tail"] = kg["tail"].str.upper()

    gt = pd.read_csv(gt_file, header=0)
    gt.columns = ["Gene1", "Gene2"] + list(gt.columns[2:])
    gt["Gene1"] = gt["Gene1"].str.upper()
    gt["Gene2"] = gt["Gene2"].str.upper()
    gt = gt[gt["Gene1"] != gt["Gene2"]].drop_duplicates(subset=["Gene1", "Gene2"])
    return expr, kg, gt


def candidate_pairs(expr: pd.DataFrame, gt: pd.DataFrame):
    """候选边 = (ground truth 中出现过的调控子) × (表达矩阵中其余全部基因)。

    与上游口径的差别(已逐行比对 utils/utility.py,照实写):
      · 上游 computeScores 的候选集是 product(unique(GT.Gene1 ∪ GT.Gene2), repeat=2),
        **两端都限定在 ground truth 里出现过的基因**,且 selfEdges=True 默认含自环;
      · 本模块把尾端放宽到表达矩阵的全部基因(自环已排除)。
    因此本模块的 AUPRC 基线率更低、数值不可与上游论文表格直接对比;
    它只用于**同一张榜单内部**各方法之间的横向比较。
    """
    genes = list(expr.index)
    regulators = [g for g in dict.fromkeys(gt["Gene1"]) if g in set(genes)]
    pairs = [(a, b) for a in regulators for b in genes if a != b]
    truth = set(map(tuple, gt[["Gene1", "Gene2"]].values))
    y = np.array([1 if p in truth else 0 for p in pairs], dtype=int)
    return regulators, pairs, y


# --------------------------------------------------------------------------- #
# 二、基线打分器(全部本机可跑)
# --------------------------------------------------------------------------- #
def score_correlation(expr, pairs, method="pearson"):
    """|相关系数| —— BEELINE 里最经典的朴素对照。"""
    X = expr.values.astype(float)
    if method == "spearman":
        X = np.apply_along_axis(rankdata, 1, X)
    C = np.corrcoef(X)
    C = np.nan_to_num(C)
    idx = {g: i for i, g in enumerate(expr.index)}
    return np.array([abs(C[idx[a], idx[b]]) for a, b in pairs])


def score_pca_dot(expr, pairs, n_comp=10, seed=SEED):
    """基因嵌入点积。

    这是对 KEGNI 读出方式的"浅层同构对照":上游 --norm -1 即用基因嵌入的点积构造
    预测边权,只不过它的嵌入来自 GAT 掩码自编码器 + 知识图嵌入联合训练,这里换成 PCA。
    深度嵌入若打不过这个线性嵌入,就说明增益不来自"深度"。
    """
    X = expr.values.astype(float)
    X = (X - X.mean(1, keepdims=True)) / (X.std(1, keepdims=True) + 1e-9)
    Z = PCA(n_components=min(n_comp, min(X.shape) - 1), random_state=seed).fit_transform(X)
    Z = Z / (np.linalg.norm(Z, axis=1, keepdims=True) + 1e-9)
    idx = {g: i for i, g in enumerate(expr.index)}
    return np.array([abs(float(Z[idx[a]] @ Z[idx[b]])) for a, b in pairs])


def score_kg_prior(kg, pairs):
    """纯知识先验:知识图上的最短路邻近度 1/(1+d),不看任何表达。

    对应 KEGNI 的"知识"一侧被单独拿出来当对照 —— 如果知识先验单独就能达到融合的分数,
    那么表达那一路其实没贡献。这条基线就是用来暴露这件事的。
    """
    G = nx.Graph()
    for h, r, t in kg[["head", "relation", "tail"]].values:
        G.add_edge(h, t)
    out = np.zeros(len(pairs))
    cache: dict[str, dict] = {}
    for i, (a, b) in enumerate(pairs):
        if a not in G or b not in G:
            continue
        if a not in cache:
            cache[a] = nx.single_source_shortest_path_length(G, a, cutoff=4)
        d = cache[a].get(b)
        out[i] = 0.0 if d is None else 1.0 / (1.0 + d)
    return out


def score_fusion(s_expr, s_kg, w=0.5):
    """知识-表达秩融合:两路各自转百分位秩后加权平均。

    透明、无需训练,是"知识嵌入应该带来增益"这一主张的最低成本实现;
    它 **不是** KEGNI 的算法,只是同一动机下的诚实下界。
    """
    r1 = rankdata(s_expr) / len(s_expr)
    r2 = rankdata(s_kg) / len(s_kg)
    return (1 - w) * r1 + w * r2


# --------------------------------------------------------------------------- #
# 三、评估(BEELINE 口径;与上游 eval.py 的差别已逐条注明,勿当作同一数字)
#
#   上游 utils/utility.py 的实现(已读源码):
#     · EarlyPrec()  返回 **裸的 early precision**(len(命中)/len(取出的边)),不除随机期望;
#     · computeScores() 的 AUPRC = auc(recall, prec) 梯形积分,不是 average_precision_score;
#   本模块用的是 BEELINE 论文更常引用的 EPR **比值** 与 sklearn 的 AP,
#   两者单调性一致但数值不同。跨模块/跨论文比较时请以各自定义为准。
# --------------------------------------------------------------------------- #
def early_precision_ratio(y, s):
    """EPR:取前 k 条(k = 真实边数)预测的精确率 / 随机期望精确率。"""
    k = int(y.sum())
    if k == 0:
        return float("nan")
    top = np.argsort(-s, kind="stable")[:k]
    ep = y[top].sum() / k
    return float(ep / (k / len(y)))


def evaluate(y, s):
    return dict(EPR=early_precision_ratio(y, s),
                AUPRC=float(average_precision_score(y, s)),
                AUROC=float(roc_auc_score(y, s)))


# --------------------------------------------------------------------------- #
# 四、KEGNI 守卫式封装(不臆造 API;缺依赖就诚实退出)
# --------------------------------------------------------------------------- #
def kegni_guard(repo: str | None, expr_file: str, kg_file: str, gt_file: str,
                outdir: str, run: bool, genes: int, max_steps: int,
                use_eval: bool = False):
    """检查真实依赖与仓库,凑齐则按核实过的命令行调用,否则打印真实安装指引。

    ★ 关于 --eval(读 train/trainer.py:128-165 得到的硬事实):
      加了 --eval,上游会把 ground truth 路径**写死**成
      ``<--dir>data/GroundTruth/TFs500/<name>/<name>-{STRING,NonSpe,ChIP}-network.csv``
      (mESC 另多一个 ``-lofgof-``;``--genes 1000`` 时整体切到 ``TFs1000/``,见 :148-160),
      其中 ``name`` 取自日志文件名 basename 的 ``split('_')[0]``,而日志名 =
      ``<input 去扩展名>_<时间戳>``(train.py:21,29 + trainer.py:89),等价于
      ``basename(--input).split('_')[0]``。也就是说 --eval 只对 BEELINE 那
      7 个数据集(命名形如 ``mESC_exp.csv``)有效,喂任意数据会直接 FileNotFoundError。
      上游自己的非 BEELINE 例子 run_pbmc_naiveCD4T.sh 也**没有**传 --eval。
      故本模块默认不加 --eval,要加请显式 --kegni-eval 并自行摆好 BEELINE 目录结构。
    """
    report: dict = {"attempted": bool(run), "ran": False}
    missing = []
    # transformers 是硬依赖:train/trainer.py:6 有 `from transformers import Trainer`
    for mod in ("torch", "dgl", "transformers"):
        try:
            __import__(mod)
        except Exception:
            missing.append(mod)
    report["missing_python_deps"] = missing

    repo_ok = bool(repo) and os.path.isfile(os.path.join(repo or "", "train.py"))
    report["repo_path"] = repo
    report["repo_found"] = repo_ok

    cmd = ["python", "train.py"] + (["--eval"] if use_eval else []) + [
           "--input", os.path.abspath(expr_file),
           "--data_path", os.path.abspath(kg_file),
           "--genes", str(genes),
           "--n_neighbors", str(KEGNI_DEFAULTS["n_neighbors"]),
           "--num_hidden", str(KEGNI_DEFAULTS["num_hidden"]),
           "--num_heads", str(KEGNI_DEFAULTS["num_heads"]),
           "--num_layers", str(KEGNI_DEFAULTS["num_layers"]),
           "--norm", str(KEGNI_DEFAULTS["norm"]),
           "--max_steps", str(max_steps),
           "--lambda_kge", str(KEGNI_DEFAULTS["lambda_kge"])]
    report["command"] = " ".join(cmd)

    if missing or not repo_ok:
        print("    [KEGNI] 未运行上游模型。原因:", end=" ")
        if missing:
            print(f"缺少 Python 依赖 {missing}", end="  ")
        if not repo_ok:
            print("未找到本地 KEGNI 仓库(--kegni-repo 指向含 train.py 的目录)", end="")
        print()
        print("    [KEGNI] 上游为命令行研究代码,无 pip 包、无 requirements.txt。")
        print("        git clone https://github.com/Lipxiao/KEGNI")
        # 依赖清单 = 对上游全部 *.py 的 import 做穷举得到,不是凭印象写的;版本上游未给。
        print("        pip install torch dgl transformers pandas numpy scikit-learn matplotlib")
        print("        # 数据(表达/ground truth/细胞类型知识图)见 https://zenodo.org/records/14028211")
        print(f"        cd KEGNI && {report['command']}")
        print("        # ↑ 训练完 embedding 落在 outputs/<input去扩展名>_<时间戳>-scg_embedding.csv")
        print("        # 上游默认 **不导出边表**(trainer.py:247 的 to_csv 被注释掉了),")
        print("        # 用本模块 --kegni-embedding <该csv> 按上游 recon(tanh+norm=-1) 复原边表并入榜。")
        report["status"] = "skipped"
        return report

    if not run:
        print("    [KEGNI] 依赖与仓库齐备。加 --run-kegni 才会真正训练(耗时,建议 GPU)。")
        report["status"] = "available_not_run"
        return report

    report["eval_flag"] = use_eval
    print(f"    [KEGNI] 执行:{report['command']}  (cwd={repo})")
    proc = subprocess.run(cmd, cwd=repo, capture_output=True, text=True)
    report.update(status="ran", returncode=proc.returncode, ran=True,
                  stdout_tail=proc.stdout[-2000:], stderr_tail=proc.stderr[-2000:])
    with open(os.path.join(outdir, "kegni_stdout.log"), "w", encoding="utf-8") as fh:
        fh.write(proc.stdout + "\n----- stderr -----\n" + proc.stderr)
    return report


def score_from_kegni_embedding(emb_file: str, pairs):
    """从上游导出的基因嵌入复原边权,复刻 train/trainer.py 的 recon(norm=-1)。

    源码依据 train/trainer.py:336-360(逐行核对过,**含 tanh**)——
        def recon(self, z):
            norm = self.args.norm
            ...
            z = torch.tanh(z)                              # :343  ★先过 tanh
            if norm == -1:  pred = torch.mm(z, z.t())      # :346-347 纯点积,不归一化
            pred = pd.DataFrame(pred.numpy(), index=list(self.model.scg2id.keys()),
                                columns=list(self.model.scg2id.keys()))   # :353
            ... melt 成 Gene1/Gene2/EdgeWeight,并去掉 Gene1 == Gene2      # :355-358
    这里用 numpy 做同一件事(tanh 后 z @ z.T),因此**不需要 torch**。
    落盘的 scg_embedding.csv 存的是 tanh **之前**的原始 z(trainer.py:234,240),
    所以复原时必须自己补上 tanh —— 漏掉会得到另一套(单调性也不同的)边权。
    嵌入文件即上游 trainer.py:240 落盘的
    ``outputs/<log-basename>-scg_embedding.csv``,其中 log-basename =
    ``<input 文件名去扩展名>_<YYYY-MM-DD-HH-MM-SS>``(trainer.py:89 取自日志文件名,
    train.py:29 带了时间戳),故实际文件名含时间戳。第一列为基因名索引。
    norm != -1 时上游走的是 cosine_similarity(trainer.py:338-342),本函数不覆盖。
    """
    z = pd.read_csv(emb_file, index_col=0, header=0)
    z.index = z.index.astype(str).str.upper()
    Zt = np.tanh(z.values.astype(float))          # 对齐 trainer.py:343
    M = Zt @ Zt.T
    idx = {g: i for i, g in enumerate(z.index)}
    return np.array([M[idx[a], idx[b]] if (a in idx and b in idx) else 0.0
                     for a, b in pairs])


# --------------------------------------------------------------------------- #
# 五、出图(无条形图:dot plot / PR 曲线 / 散点 / slopegraph)
# --------------------------------------------------------------------------- #
def fig_metric_dotplot(res: pd.DataFrame, out):
    """方法 × 指标 dot plot:点大小=组内归一化,颜色=原值。取代分组条形图。"""
    metrics = ["EPR", "AUPRC", "AUROC"]
    methods = list(res.index)
    fig, axes = plt.subplots(1, 3, figsize=(NATURE_W2, 0.42 * len(methods) + 1.7),
                             sharey=True)
    for ax, m in zip(axes, metrics):
        v = res[m].values
        rel = (v - v.min()) / (v.max() - v.min() + 1e-12)
        ax.scatter(v, np.arange(len(methods)), s=60 + 340 * rel, c=v,
                   cmap=CMAP_CONT, edgecolor="black", linewidth=0.7, zorder=3)
        for i, val in enumerate(v):
            ax.hlines(i, v.min() - 0.02 * (v.max() - v.min() + 1e-9), val,
                      color="#BBBBBB", linewidth=0.8, zorder=1)
        ax.set_title(m)
        ax.set_xlabel(m)
        ax.margins(x=0.22, y=0.12)
    axes[0].set_yticks(np.arange(len(methods)))
    axes[0].set_yticklabels(methods)
    fig.suptitle("GRN inference benchmark (synthetic demo)", y=1.03,
                 fontweight="bold")
    fig.tight_layout()
    save_fig(fig, out)
    plt.close(fig)


def fig_pr_curves(y, scores: dict, res: pd.DataFrame, out):
    fig, ax = plt.subplots(figsize=(NATURE_W1 + 1.4, NATURE_W1 + 0.4))
    colors = pal(len(scores))
    for (name, s), c in zip(scores.items(), colors):
        pr, rc, _ = precision_recall_curve(y, s)
        ax.plot(rc, pr, color=c, linewidth=1.7,
                label=f"{name} (AUPRC {res.loc[name, 'AUPRC']:.3f})")
    ax.axhline(y.mean(), color="black", linestyle="--", linewidth=1.0,
               label=f"random ({y.mean():.3f})")
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_title("Precision-recall on held-out ground-truth edges")
    ax.legend(loc="upper right", fontsize=7)
    fig.tight_layout()
    save_fig(fig, out)
    plt.close(fig)


def fig_complementarity(y, s_expr, s_kg, out):
    """表达证据 vs 知识证据 散点:真实边为何需要两路互补,一眼看清。"""
    rx = rankdata(s_expr) / len(s_expr)
    ry = rankdata(s_kg) / len(s_kg)
    jit = np.random.default_rng(SEED).normal(0, 0.008, len(ry))
    fig, ax = plt.subplots(figsize=(NATURE_W1 + 1.2, NATURE_W1 + 0.6))
    ax.scatter(rx[y == 0], ry[y == 0] + jit[y == 0], s=7, c="#C9CDD4",
               alpha=0.55, linewidth=0, label="non-edge")
    ax.scatter(rx[y == 1], ry[y == 1] + jit[y == 1], s=34, c=pal(3)[0],
               edgecolor="black", linewidth=0.5, label="true edge", zorder=3)
    ax.set_xlabel("Expression evidence (percentile rank)")
    ax.set_ylabel("Knowledge-graph proximity (percentile rank)")
    ax.set_title("Two evidence axes are complementary")
    ax.legend(loc="lower left", fontsize=8)
    fig.tight_layout()
    save_fig(fig, out)
    plt.close(fig)


def fig_slopegraph(res: pd.DataFrame, left: str, right: str, out):
    """slopegraph:表达单路 → 知识融合,三个指标各自怎么走。

    三个指标量纲不同,各占一条独立泳道(避免归一化后标签重叠);
    斜率 = 相对变化幅度,方向即"知识有没有帮上忙",两端标原始数值。
    """
    metrics = ["EPR", "AUPRC", "AUROC"]
    fig, ax = plt.subplots(figsize=(NATURE_W1 + 1.6, NATURE_W1 + 0.5))
    colors = pal(len(metrics))
    for i, (m, c) in enumerate(zip(metrics, colors)):
        a, b = float(res.loc[left, m]), float(res.loc[right, m])
        lane = len(metrics) - 1 - i                       # 自上而下 EPR/AUPRC/AUROC
        rel = (b - a) / max(abs(a), abs(b), 1e-9)         # 相对变化,-1..1
        y0, y1 = lane, lane + 0.62 * float(np.clip(rel, -1, 1))
        ax.plot([0, 1], [y0, y1], "-o", color=c, linewidth=2.0, markersize=6,
                markeredgecolor="black", markeredgewidth=0.6, zorder=3)
        ax.text(-0.05, y0, f"{m}  {a:.3f}", ha="right", va="center",
                fontsize=8.5, color=c)
        ax.text(1.05, y1, f"{b:.3f}  ({rel * 100:+.0f}%)", ha="left", va="center",
                fontsize=8.5, color=c)
        ax.hlines(lane, 0, 1, color="#E4E4E4", linewidth=0.8, zorder=1)
    ax.set_xlim(-0.62, 1.72)
    ax.set_ylim(-0.9, len(metrics) - 0.2)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([left, right], fontsize=8.5)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.tick_params(axis="x", length=0)
    ax.set_title("Adding knowledge to expression evidence")
    fig.tight_layout()
    save_fig(fig, out)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# 六、主流程
# --------------------------------------------------------------------------- #
def main():
    p = argparse.ArgumentParser(description="583 KEGNI 知识图增强 GRN 推断(基线 + 守卫封装)")
    d = os.path.join(HERE, "example_data")
    p.add_argument("--expr", default=os.path.join(d, "expression.csv"),
                   help="表达矩阵 CSV,行=基因 列=细胞(KEGNI --input 同格式)")
    p.add_argument("--kg", default=os.path.join(d, "knowledge_graph.tsv"),
                   help="知识图三元组 TSV,无表头 head/relation/tail(KEGNI --data_path 同格式)")
    p.add_argument("--truth", default=os.path.join(d, "ground_truth_network.csv"),
                   help="ground truth 网络 CSV,表头 Gene1,Gene2")
    p.add_argument("--outdir", default=os.path.join(HERE, "results"))
    p.add_argument("--assets", default=os.path.join(HERE, "assets"))
    p.add_argument("--fusion-weight", type=float, default=0.5,
                   help="秩融合中知识先验的权重 w(0=纯表达,1=纯知识)")
    p.add_argument("--n-comp", type=int, default=10, help="PCA 基因嵌入维度")
    p.add_argument("--kegni-repo", default=None, help="本地 KEGNI 仓库路径(含 train.py)")
    p.add_argument("--run-kegni", action="store_true", help="依赖齐备时真正训练 KEGNI")
    p.add_argument("--kegni-pred", default=None,
                   help="已有的 KEGNI 预测边 CSV(Gene1,Gene2,EdgeWeight),并入同一张榜单")
    p.add_argument("--kegni-embedding", default=None,
                   help="上游 outputs/<basename>-scg_embedding.csv;按 recon(norm=-1) 点积复原边表并入榜")
    p.add_argument("--kegni-eval", action="store_true",
                   help="给上游命令行加 --eval(仅对 BEELINE 命名的数据有效,见 kegni_guard 文档)")
    p.add_argument("--genes", type=int, default=None, help="传给 KEGNI 的 --genes")
    p.add_argument("--max-steps", type=int, default=KEGNI_DEFAULTS["max_steps"])
    p.add_argument("--seed", type=int, default=SEED)
    args = p.parse_args()

    np.random.seed(args.seed)
    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(args.assets, exist_ok=True)
    set_pub_style(base_size=9)

    print("Step 1  读取表达矩阵 / 知识图 / ground truth")
    expr, kg, gt = load_inputs(args.expr, args.kg, args.truth)
    regs, pairs, y = candidate_pairs(expr, gt)
    print(f"        基因 {expr.shape[0]} × 细胞 {expr.shape[1]} · 知识图三元组 {len(kg)} "
          f"· 调控子 {len(regs)} · 候选边 {len(pairs)} · 真实边 {int(y.sum())}")

    print("Step 2  基线打分(相关性 / 秩相关 / PCA 嵌入点积 / 纯知识先验 / 知识-表达融合)")
    scores: dict[str, np.ndarray] = {}
    scores["Pearson |r|"] = score_correlation(expr, pairs, "pearson")
    scores["Spearman |rho|"] = score_correlation(expr, pairs, "spearman")
    scores["PCA embedding dot"] = score_pca_dot(expr, pairs, args.n_comp, args.seed)
    scores["KG prior only"] = score_kg_prior(kg, pairs)
    # 融合用的是表达一侧 **最强** 的基线(Pearson),不是最弱的那个 ——
    # 否则"知识带来增益"会是靠挑软柿子挑出来的。
    scores["KG + expression (rank fusion)"] = score_fusion(
        scores["Pearson |r|"], scores["KG prior only"], args.fusion_weight)

    if args.kegni_pred:
        print(f"        并入外部 KEGNI 预测:{args.kegni_pred}")
        pred = pd.read_csv(args.kegni_pred)
        pred.columns = ["Gene1", "Gene2", "EdgeWeight"] + list(pred.columns[3:])
        pred["Gene1"] = pred["Gene1"].str.upper()
        pred["Gene2"] = pred["Gene2"].str.upper()
        lut = dict(zip(zip(pred["Gene1"], pred["Gene2"]), pred["EdgeWeight"]))
        scores["KEGNI (upstream)"] = np.array([float(lut.get(pr, 0.0)) for pr in pairs])

    if args.kegni_embedding:
        print(f"        由 KEGNI 基因嵌入复原边表(recon norm=-1):{args.kegni_embedding}")
        scores["KEGNI (from embedding)"] = score_from_kegni_embedding(
            args.kegni_embedding, pairs)

    print("Step 3  评估(EPR / AUPRC / AUROC;BEELINE 语义,与上游 eval.py 的定义差异见源码注释)")
    res = pd.DataFrame({k: evaluate(y, v) for k, v in scores.items()}).T
    res.index.name = "method"
    print(res.round(4).to_string())
    res.to_csv(os.path.join(args.outdir, "benchmark_metrics.csv"))

    # 预测边表按上游格式落盘:Gene1,Gene2,EdgeWeight —— 可直接喂 KEGNI 的 eval.py
    best = res["AUPRC"].idxmax()
    edge_tab = pd.DataFrame({"Gene1": [a for a, _ in pairs],
                             "Gene2": [b for _, b in pairs],
                             "EdgeWeight": scores[best]})
    edge_tab.sort_values("EdgeWeight", ascending=False).to_csv(
        os.path.join(args.outdir, "predicted_edges_best.csv"), index=False)
    for name, s in scores.items():
        # 整名 slug 化(只取首词会让 "KG prior only" 与 "KG + expression" 撞同一个文件名)
        slug = "".join(ch if ch.isalnum() else "_" for ch in name.lower()).strip("_")
        while "__" in slug:
            slug = slug.replace("__", "_")
        pd.DataFrame({"Gene1": [a for a, _ in pairs], "Gene2": [b for _, b in pairs],
                      "EdgeWeight": s}).sort_values(
            "EdgeWeight", ascending=False).to_csv(
            os.path.join(args.outdir, f"edges_{slug}.csv"), index=False)

    print("Step 4  KEGNI 上游守卫检查")
    kg_report = kegni_guard(args.kegni_repo, args.expr, args.kg, args.truth,
                            args.outdir, args.run_kegni,
                            args.genes if args.genes else expr.shape[0],
                            args.max_steps, args.kegni_eval)

    print("Step 5  出图(dot plot / PR 曲线 / 互补性散点 / slopegraph)")
    figs = [("fig1_benchmark_dotplot", lambda o: fig_metric_dotplot(res, o)),
            ("fig2_pr_curves", lambda o: fig_pr_curves(y, scores, res, o)),
            ("fig3_evidence_complementarity",
             lambda o: fig_complementarity(y, scores["Pearson |r|"],
                                           scores["KG prior only"], o)),
            ("fig4_knowledge_gain_slopegraph",
             lambda o: fig_slopegraph(res, "Pearson |r|",
                                      "KG + expression (rank fusion)", o))]
    for stem, fn in figs:
        fn(os.path.join(args.outdir, stem))
        for ext in (".png", ".pdf"):
            src = os.path.join(args.outdir, stem + ext)
            if os.path.exists(src) and ext == ".png":
                shutil.copyfile(src, os.path.join(args.assets, stem + ext))
        print(f"        {stem}.png / .pdf")

    # 依赖版本快照(等价于 R 的 sessionInfo,落盘保证可复现)
    import platform
    import sklearn
    import scipy
    import matplotlib
    session = {"python": platform.python_version(), "numpy": np.__version__,
               "pandas": pd.__version__, "scipy": scipy.__version__,
               "scikit-learn": sklearn.__version__, "networkx": nx.__version__,
               "matplotlib": matplotlib.__version__}

    summary = {"session_info": session,
               "seed": args.seed, "n_genes": int(expr.shape[0]),
               "n_cells": int(expr.shape[1]), "n_kg_triples": int(len(kg)),
               "n_candidate_pairs": int(len(pairs)), "n_true_edges": int(y.sum()),
               "best_by_auprc": best, "metrics": res.round(6).to_dict(orient="index"),
               "kegni": kg_report,
               "citation": ("Li P et al. KEGNI: knowledge graph enhanced framework for "
                            "gene regulatory network inference. Genome Biol 2025;26:294. "
                            "PMID 40983951; doi:10.1186/s13059-025-03780-7")}
    with open(os.path.join(args.outdir, "583_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)
    print(f"完成 → {args.outdir}")


if __name__ == "__main__":
    main()
