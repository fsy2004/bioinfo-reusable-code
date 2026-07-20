# -*- coding: utf-8 -*-
# =============================================================================
# 586 · PSGRN — 单细胞扰动数据的基因调控网络推断(自训练 + 合成金标准)
# -----------------------------------------------------------------------------
# PSGRN 是 GSK CausalBench Challenge 的优胜方案(Song X, Deng K, Chen M, Guan Y.
# Sci Adv 2026;12(18):eaeb3376; doi:10.1126/sciadv.aeb3376; PMID 42054465)。
# 措辞注:上游 README 原文自述 "our winning solution on the CausalBench Challenge";
#         论文摘要原文措辞为 "a top-performing method in the CausalBench Challenge"。
# 上游仓库 https://github.com/GuanLab/PSGRN
#
# ★核心思想(读上游 src/main.py 原文得到,非臆造):
#   ① 用【相关性 > T】的基因对当作"合成金标准"正样本(noisy label);
#   ② 为每个有序基因对 (gene1 -> gene2) 造 4 个【扰动特征】:
#      非靶向对照下 gene1 均值、gene2 均值、gene1 敲除后 gene1 自身均值、
#      gene1 敲除后 gene2 均值;
#   ③ 用 GBDT 在这些噪声标签上自训练,再对全部基因对打分取 top-N 作为网络。
#   直觉:相关性只提供关联,分类器要用"干预后的响应"才能把它重排成因果方向。
#
# ★本机运行状态(诚实标注,见 README):
#   - 上游官方入口是 CausalBench 的 causalscbench/apps/main_app.py,本机【未装 causalscbench】
#     → run_psgrn_upstream() 是【守卫式封装】,缺包时优雅退出并打印真实安装/运行命令,
#       绝不伪造 causalscbench 的调用结果。
#   - 本模块另带一条【本地可跑的忠实复现】(psgrn_selftrain),逐函数对照上游 src/main.py
#     重写(get_topK_pairs / create_dataset / train_lgb 同名同结构),LightGBM 超参
#     LGB_PARAMS 逐字照抄上游;本机已装 lightgbm 4.6.0,故走原路而非替代实现。
#   - 以及两条【朴素基线】:共表达相关性、单变量扰动效应(库规矩:任何"更好"都要有对照)。
#   - ⚠ 复现的是【算法】,不是上游在 CausalBench 上的评测数字;别拿本模块结果替代原文报告值。
#
# Turnkey: python 586_psgrn_grn_inference.py
# 换数据 : python 586_psgrn_grn_inference.py --input data/你的.csv --outdir results/run1
# 复用 _framework/pubstyle.py;图中文字英文,注释中文。
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

try:                                    # Windows 控制台默认 GBK,避免中文进度条乱码
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

SEED = 0

# ---- 载入顶刊绘图框架(向上搜 _framework) ----------------------------------
HERE = Path(__file__).resolve().parent
for _up in [HERE, *HERE.parents]:
    if (_up / "_framework" / "pubstyle.py").exists():
        sys.path.insert(0, str(_up / "_framework"))
        break
try:
    from pubstyle import set_pub_style, save_fig, pal, CMAP_CONT, NATURE_W1, NATURE_W2
except Exception:                                          # 框架缺失时的最小降级
    CMAP_CONT, NATURE_W1, NATURE_W2 = "viridis", 3.5, 7.0

    def set_pub_style(*a, **k):
        pass

    def pal(n=None, name="npg"):
        base = ["#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F"]
        return base[:n] if n else base

    def save_fig(fig, f, dpi=300):
        Path(f).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(str(f) + ".pdf")
        fig.savefig(str(f) + ".png", dpi=dpi)

EXAMPLE = HERE / "example_data"
RESULTS = HERE / "results"
ASSETS = HERE / "assets"
OBS_LABEL = "non-targeting"     # 上游 src/main.py 用的对照标签,保持一致


# =============================================================================
# 0. 合成示例数据(synthetic, for demo only)
# =============================================================================
def make_example(n_genes: int = 40, n_perturbed: int = 24, n_edges: int = 60,
                 n_obs_cells: int = 300, n_cells_per_ko: int = 40, seed: int = SEED):
    """造一份小型 CRISPRi 风格扰动矩阵 + 已知真值网络。

    设计要点(为了让"基线 vs PSGRN"的比较有意义):
      - 存在共享潜因子 → 制造大量【相关但无调控关系】的基因对(相关性基线的陷阱);
      - 真实调控边只在【敲除 gene1 的细胞里】才让 gene2 发生偏移(干预信号);
      - 部分真实边的两端在对照细胞中几乎不相关 → 相关性基线必然漏掉。
    """
    rng = np.random.default_rng(seed)
    genes = [f"G{i:02d}" for i in range(n_genes)]
    perturbed = genes[:n_perturbed]

    # 真值网络:from 必须是被敲除过的基因(否则干预信号无从观测)
    edges = set()
    while len(edges) < n_edges:
        a = perturbed[rng.integers(n_perturbed)]
        b = genes[rng.integers(n_genes)]
        if a != b:
            edges.add((a, b))
    edges = sorted(edges)
    W = pd.DataFrame(0.0, index=genes, columns=genes)
    for a, b in edges:
        W.loc[a, b] = rng.choice([-1, 1]) * rng.uniform(0.8, 2.0)

    base = rng.uniform(2.0, 6.0, n_genes)                  # 基线表达
    n_factor = 4
    loading = rng.normal(0, 0.25, (n_genes, n_factor))     # 共享潜因子载荷 → 混杂相关

    def _sample(k, ko_gene=None):
        F = rng.normal(0, 1.0, (k, n_factor))
        X = base[None, :] + F @ loading.T + rng.normal(0, 1.5, (k, n_genes))
        if ko_gene is not None:
            j = genes.index(ko_gene)
            X[:, j] *= 0.15                                # 敲除:自身表达塌掉
            X += W.loc[ko_gene].values[None, :]            # 下游按边权偏移
        return np.clip(X, 0, None)

    blocks, labels = [_sample(n_obs_cells)], [OBS_LABEL] * n_obs_cells
    for g in perturbed:
        blocks.append(_sample(n_cells_per_ko, ko_gene=g))
        labels += [g] * n_cells_per_ko

    expr = pd.DataFrame(np.vstack(blocks), columns=genes)
    expr.insert(0, "intervention", labels)

    EXAMPLE.mkdir(parents=True, exist_ok=True)
    head = ("# synthetic, for demo only — 非真实数据,仅用于冒烟测试与展示图\n"
            "# 行=细胞;intervention=该细胞被敲除的基因,'non-targeting'=对照;其余列=基因表达\n")
    p1 = EXAMPLE / "perturb_expression.csv"
    with open(p1, "w", newline="", encoding="utf-8") as fh:
        fh.write(head)
        expr.to_csv(fh, index=False)
    p2 = EXAMPLE / "ground_truth_edges.csv"
    with open(p2, "w", newline="", encoding="utf-8") as fh:
        fh.write("# synthetic, for demo only — 生成本示例时使用的真值调控边\n")
        pd.DataFrame(edges, columns=["From", "To"]).to_csv(fh, index=False)
    return p1, p2


# =============================================================================
# 1. PSGRN 忠实复现(本地版)
#    结构逐步对照上游 https://github.com/GuanLab/PSGRN/blob/main/src/main.py
# =============================================================================
def get_topK_pairs(expression_matrix: pd.DataFrame, T: float = 0.1):
    """合成金标准:细胞级相关性 > T 的有序基因对当正样本。

    对照上游 `get_topK_pairs(expression_matrix, T=0.5)`(实际调用时 T=0.1):
    gene1 被敲除过时,把对照细胞下采样到与干预细胞等量后拼接再算相关,
    这样相关性里带上了干预造成的变异。行索引 = intervention 标签。
    """
    cols = list(expression_matrix.columns)
    obs = {g: expression_matrix.loc[OBS_LABEL, g] for g in cols}
    self_inv = {g: expression_matrix.loc[g, g] for g in cols if g in expression_matrix.index}

    rows = []
    for g1 in cols:
        for g2 in cols:
            if g1 == g2:
                continue
            if g1 in expression_matrix.index:
                e1i, e2i = self_inv[g1], expression_matrix.loc[g1, g2]
                e1 = pd.concat([obs[g1].sample(e1i.shape[0], random_state=0), e1i])
                e2 = pd.concat([obs[g2].sample(e2i.shape[0], random_state=0), e2i])
            else:
                e1, e2 = obs[g1], obs[g2]
            rows.append([g1, g2, np.abs(e1.corr(e2))])
    corrs = pd.DataFrame(rows, columns=["From", "To", "weights"]).sort_values(
        "weights", ascending=False)
    topK = [(i, j) for i, j in corrs.loc[corrs["weights"] > T, ["From", "To"]].values]
    return corrs, topK


def create_dataset(expression_matrix: pd.DataFrame, pairs) -> pd.DataFrame:
    """造 4 特征 + label 的成对表(与上游 `create_dataset` 同构)。

    特征(按 intervention 分组求均值后取标量):
      f0 对照下 gene1 均值 · f1 对照下 gene2 均值
      f2 gene1 敲除后 gene1 自身均值(缺失记 0) · f3 gene1 敲除后 gene2 均值(缺失记 NaN)
    """
    genes = list(expression_matrix.columns)
    summ = expression_matrix.groupby(expression_matrix.index).mean()
    obs = {g: summ.loc[OBS_LABEL, g] for g in genes}
    self_inv = {g: summ.loc[g, g] for g in genes if g in summ.index}
    pairs = set(pairs)

    rec, idx = [], []
    for g1 in genes:
        for g2 in genes:
            if g1 == g2:
                continue
            idx.append(f"{g1}_{g2}")
            rec.append([obs[g1], obs[g2], self_inv.get(g1, 0),
                        summ.loc[g1, g2] if g1 in summ.index else np.nan,
                        1 if (g1, g2) in pairs else 0])
    ds = pd.DataFrame(rec, index=idx, columns=["f0", "f1", "f2", "f3", "label"])
    return ds


# 上游 src/main.py 中 Custom.__call__ 里写死的 LightGBM 超参(逐字照抄)
LGB_PARAMS = {
    "boosting_type": "gbdt",
    "objective": "binary",
    "metric": "binary_logloss",
    "num_leaves": 5,
    "max_depth": 2,
    "min_data_in_leaf": 5,
    "learning_rate": 0.05,
    "min_gain_to_split": 0.01,
    "num_iterations": 1000,
    "num_threads": 8,
    "verbose": -1,          # 上游为 0;这里压掉 LightGBM 的刷屏日志,不影响模型
}


def _fit_gbdt(ds: pd.DataFrame, seed: int = SEED):
    """自训练用的 GBDT。

    本机已装 lightgbm 4.6.0 → 走上游 `train_lgb()` 的原路:
    `lgb.Dataset(X, y)` → `lgb.train(params, train_set, keep_training_booster=True)`,
    超参用 LGB_PARAMS(上游原值)。lightgbm 缺失时退回 sklearn
    HistGradientBoostingClassifier(参数对齐到同量级,属替代实现)。
    """
    X, y = ds.drop(columns="label"), ds["label"]
    try:
        import lightgbm as lgb
    except ImportError:
        from sklearn.ensemble import HistGradientBoostingClassifier
        clf = HistGradientBoostingClassifier(
            max_depth=2, max_leaf_nodes=5, learning_rate=0.05, max_iter=300,
            min_samples_leaf=5, random_state=seed)
        clf.fit(X, y)
        return ("sklearn-fallback", lambda M: clf.predict_proba(M)[:, 1])
    gbm = lgb.train(params=LGB_PARAMS, train_set=lgb.Dataset(X, y),
                    keep_training_booster=True)
    return ("lightgbm", gbm.predict)


def psgrn_selftrain(expr: pd.DataFrame, T: float = 0.1, seed: int = SEED):
    """PSGRN 主流程(本地复现):相关性金标准 → z-score → 4 特征 → GBDT 自训练 → 打分。"""
    corrs, topK = get_topK_pairs(expr, T=T)
    z = (expr - expr.mean(axis=0)) / expr.std(axis=0)       # 上游在造特征前做整体标准化
    ds = create_dataset(z, topK)
    backend, predict = _fit_gbdt(ds.sample(frac=1, random_state=seed), seed=seed)
    score = np.asarray(predict(ds.drop(columns="label")))
    out = pd.DataFrame({"pair": ds.index, "score": score})
    out[["From", "To"]] = out["pair"].str.split("_", n=1, expand=True)
    return out[["From", "To", "score"]], corrs, len(topK), backend


# =============================================================================
# 2. 朴素基线(库规矩:任何"更好"都必须有可跑对照)
# =============================================================================
def baseline_coexpression(expr: pd.DataFrame) -> pd.DataFrame:
    """基线 A:仅用对照细胞的 |Pearson| 共表达。经典 GRN 起点,无因果方向。"""
    obs = expr.loc[OBS_LABEL]
    C = obs.corr().abs()
    np.fill_diagonal(C.values, np.nan)
    s = C.stack().reset_index()
    s.columns = ["From", "To", "score"]
    return s


def baseline_perturbation_effect(expr: pd.DataFrame) -> pd.DataFrame:
    """基线 B:单变量干预效应 |mean(g2 | KO g1) - mean(g2 | ctrl)| / sd(ctrl)。

    这是最"显然"的因果基线;PSGRN 若不能超过它,就不该声称更好。
    """
    genes = list(expr.columns)
    ctrl = expr.loc[OBS_LABEL]
    mu, sd = ctrl.mean(), ctrl.std().replace(0, np.nan)
    ko = expr.groupby(expr.index).mean()
    rec = []
    for g1 in genes:
        if g1 not in ko.index or g1 == OBS_LABEL:
            eff = pd.Series(0.0, index=genes)
        else:
            eff = ((ko.loc[g1] - mu) / sd).abs()
        for g2 in genes:
            if g1 != g2:
                rec.append([g1, g2, float(eff[g2])])
    return pd.DataFrame(rec, columns=["From", "To", "score"])


# =============================================================================
# 3. 上游真包的守卫式封装(不装包、不伪造)
# =============================================================================
def run_psgrn_upstream(expression_matrix=None, interventions=None, gene_names=None):
    """调用上游 PSGRN 官方实现。缺依赖时优雅退出并打印真实命令。

    上游真实接口(读自 https://github.com/GuanLab/PSGRN/blob/main/src/main.py):
        class Custom(AbstractInferenceModel):
            __call__(self, expression_matrix: np.array, interventions: List[str],
                     gene_names: List[str], training_regime: TrainingRegime,
                     seed: int = 0) -> List[Tuple[str, str]]
    官方推荐入口是 CausalBench 的 app,而不是直接 import(README 原文命令):
        export PYTHONPATH="./"
        python causalscbench/apps/main_app.py --dataset_name "weissmann_rpe1" \\
            --model_name "custom" --inference_function_file_path "./src/main.py" ...
    """
    missing = []
    for mod in ("causalscbench", "lightgbm"):
        try:
            __import__(mod)
        except Exception:
            missing.append(mod)
    if missing:
        return {
            "status": "skipped",
            "missing": missing,
            "install": ('conda create -n causal python=3.10 && conda activate causal && '
                        'pip install causalbench==1.1.2 lightgbm && pip uninstall causalbench -y && '
                        'git clone https://github.com/GuanLab/PSGRN'),
            "run": ('export PYTHONPATH="./"; python causalscbench/apps/main_app.py '
                    '--dataset_name weissmann_rpe1 --model_name custom '
                    '--inference_function_file_path ./src/main.py --do_filter'),
            "note": ("上游 Custom.__call__ 签名已按 src/main.py 原文记录;此处不模拟其返回值。"
                     "N(输出边数)在 src/main.py 内改:PSGRN 1K = N=1000, 5K = N=5000。"),
        }
    return {"status": "deps_present",
            "next": ("依赖齐全,但官方入口是 causalscbench/apps/main_app.py,需要 CausalBench "
                     "数据目录;请按上游 README/run.sh 运行,不要绕过 benchmark 评测层。")}


# =============================================================================
# 4. 评估与出图
# =============================================================================
def evaluate(scores: pd.DataFrame, truth: set, ks=(10, 25, 50, 100)) -> dict:
    """全基因对排序评估:AUPRC + precision@K(有真值时才有意义)。"""
    from sklearn.metrics import average_precision_score
    s = scores.sort_values("score", ascending=False).reset_index(drop=True)
    y = np.array([(a, b) in truth for a, b in zip(s["From"], s["To"])], dtype=int)
    res = {"AUPRC": float(average_precision_score(y, s["score"].values)),
           "n_pairs": int(len(s)), "n_true": int(y.sum()),
           "random_AUPRC": float(y.mean())}
    for k in ks:
        res[f"P@{k}"] = float(y[:k].mean())
    res["_y"] = y
    res["_s"] = s
    return res


def _pr_curve(y):
    """按排序累积算 precision/recall(用于 PR 曲线)。"""
    tp = np.cumsum(y)
    prec = tp / np.arange(1, len(y) + 1)
    rec = tp / max(tp[-1], 1)
    return rec, prec


def figures(ev: dict, ks, outdir: Path):
    """三张图:PR 曲线 / precision@K dumbbell / PSGRN 打分热图。★无条形图。"""
    import matplotlib.pyplot as plt
    set_pub_style(base_size=10)
    cols = pal(3, "npg")
    names = list(ev.keys())
    outs = []

    # --- fig1: PR 曲线(排序质量的完整视图) ---
    fig, ax = plt.subplots(figsize=(NATURE_W1, NATURE_W1 * 0.85))
    for c, n in zip(cols, names):
        rec, prec = _pr_curve(ev[n]["_y"])
        ax.plot(rec, prec, lw=1.6, color=c, label=f"{n} (AUPRC={ev[n]['AUPRC']:.3f})")
    ax.axhline(ev[names[0]]["random_AUPRC"], ls=":", lw=1.0, color="grey")
    ax.text(0.02, ev[names[0]]["random_AUPRC"], " random", va="bottom", fontsize=7, color="grey")
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_title("Edge-ranking performance")
    ax.legend(loc="lower left", bbox_to_anchor=(0.28, 0.14), fontsize=6.5)
    save_fig(fig, outdir / "fig1_pr_curves")
    plt.close(fig)
    outs.append("fig1_pr_curves.png")

    # --- fig2: precision@K dumbbell(最好基线 → PSGRN 的位移) ---
    ref = "PSGRN (self-training)"
    others = [n for n in names if n != ref]
    fig, ax = plt.subplots(figsize=(NATURE_W1 * 1.25, NATURE_W1 * 0.9))
    ypos = np.arange(len(ks))
    best = [max(ev[o][f"P@{k}"] for o in others) for k in ks]
    mine = [ev[ref][f"P@{k}"] for k in ks]
    for i, (b, m) in enumerate(zip(best, mine)):
        ax.plot([b, m], [i, i], color="#BBBBBB", lw=2.0, zorder=1, solid_capstyle="round")
    ax.scatter(best, ypos, s=52, color=cols[1], zorder=3, label="best baseline", edgecolor="white")
    ax.scatter(mine, ypos, s=52, color=cols[0], zorder=3, label="PSGRN", edgecolor="white")
    for i, o in enumerate(others):
        ax.scatter([ev[o][f"P@{k}"] for k in ks], ypos, s=16, color="#888888",
                   zorder=2, alpha=0.7,
                   label="each baseline" if i == 0 else None)
    ax.set_yticks(ypos)
    ax.set_yticklabels([f"top {k}" for k in ks])
    ax.set_xlabel("Precision @ K")
    ax.set_title("PSGRN vs naive baselines")
    ax.set_xlim(-0.03, 1.03)
    ax.legend(loc="lower right", fontsize=7)
    save_fig(fig, outdir / "fig2_precision_at_k")
    plt.close(fig)
    outs.append("fig2_precision_at_k.png")

    # --- fig3: PSGRN 打分矩阵热图 + 真值边散点覆盖 ---
    s = ev[ref]["_s"]
    genes = sorted(set(s["From"]) | set(s["To"]))
    M = s.pivot_table(index="From", columns="To", values="score").reindex(
        index=genes, columns=genes)
    fig, ax = plt.subplots(figsize=(NATURE_W2 * 0.62, NATURE_W2 * 0.58))
    im = ax.imshow(M.values, cmap=CMAP_CONT, aspect="auto", vmin=0, vmax=1)
    truth_pts = [(genes.index(a), genes.index(b))
                 for (a, b), t in zip(zip(s["From"], s["To"]), ev[ref]["_y"]) if t]
    if truth_pts:
        ax.scatter([b for a, b in truth_pts], [a for a, b in truth_pts],
                   s=9, facecolors="none", edgecolors="#DC0000", linewidths=0.7,
                   label="true edge")
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), fontsize=7,
                  framealpha=0.0)
    ax.set_xlabel("Target gene")
    ax.set_ylabel("Regulator gene")
    ax.set_title("PSGRN edge score matrix")
    ax.set_xticks(range(0, len(genes), 5))
    ax.set_xticklabels(genes[::5], rotation=90, fontsize=6)
    ax.set_yticks(range(0, len(genes), 5))
    ax.set_yticklabels(genes[::5], fontsize=6)
    fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02, label="predicted edge score")
    save_fig(fig, outdir / "fig3_score_matrix")
    plt.close(fig)
    outs.append("fig3_score_matrix.png")
    return outs


# =============================================================================
# 5. 主流程
# =============================================================================
def _session_info() -> dict:
    """依赖版本快照,落进 summary json 便于复现。"""
    import platform
    info = {"python": platform.python_version(), "platform": platform.platform()}
    for m in ("numpy", "pandas", "sklearn", "lightgbm", "matplotlib"):
        try:
            info[m] = __import__(m).__version__
        except Exception:
            info[m] = "not installed"
    return info



def main():
    ap = argparse.ArgumentParser(description="586 PSGRN — perturbational GRN inference")
    ap.add_argument("--input", default=str(EXAMPLE / "perturb_expression.csv"),
                    help="细胞×基因表达 csv,首列 intervention('non-targeting'=对照)")
    ap.add_argument("--truth", default=str(EXAMPLE / "ground_truth_edges.csv"),
                    help="可选真值边 csv(From,To);无真值则只出打分,不做评估")
    ap.add_argument("--outdir", default=str(RESULTS))
    ap.add_argument("--T", type=float, default=0.1, help="合成金标准的相关性阈值(上游默认 0.1)")
    ap.add_argument("--topN", type=int, default=100, help="导出的 top-N 边(上游 PSGRN 1K=1000)")
    ap.add_argument("--run-upstream", action="store_true", help="尝试上游官方 PSGRN(需装包)")
    ap.add_argument("--seed", type=int, default=SEED)
    a = ap.parse_args()

    outdir = Path(a.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    ASSETS.mkdir(parents=True, exist_ok=True)

    inp = Path(a.input)
    if not inp.exists() and inp == EXAMPLE / "perturb_expression.csv":
        print("[586] example_data 缺失,现场生成合成扰动数据")
        make_example(seed=a.seed)

    print(f"[586] Step 1 读入 {inp}")
    expr = pd.read_csv(inp, comment="#", encoding="utf-8")
    expr = expr.set_index("intervention")
    print(f"       cells={expr.shape[0]} genes={expr.shape[1]} "
          f"perturbed={expr.index.nunique() - 1} ctrl={(expr.index == OBS_LABEL).sum()}")
    if OBS_LABEL not in expr.index:
        sys.exit(f"输入缺少对照细胞(intervention == '{OBS_LABEL}'),PSGRN 无法构造特征")

    print("[586] Step 2 朴素基线 A:对照细胞共表达 |Pearson|")
    sc_co = baseline_coexpression(expr)
    print("[586] Step 3 朴素基线 B:单变量干预效应")
    sc_pe = baseline_perturbation_effect(expr)
    print(f"[586] Step 4 PSGRN 自训练(T={a.T})")
    sc_ps, corrs, n_pos, backend = psgrn_selftrain(expr, T=a.T, seed=a.seed)
    print(f"       GBDT backend = {backend}")
    print(f"       合成金标准正样本 {n_pos} 对 / 共 {len(sc_ps)} 有序基因对")

    methods = {"Co-expression |r|": sc_co,
               "Perturbation effect": sc_pe,
               "PSGRN (self-training)": sc_ps}
    for n, df in methods.items():
        df.sort_values("score", ascending=False).head(a.topN).to_csv(
            outdir / f"edges_{n.split()[0].lower().strip('|')}.csv", index=False)

    summary = {
        "input": str(inp), "n_cells": int(expr.shape[0]), "n_genes": int(expr.shape[1]),
        "T": a.T, "seed": a.seed, "n_synthetic_gold_positives": int(n_pos),
        "gbdt_backend": backend, "lgb_params": LGB_PARAMS if backend == "lightgbm" else None,
        "session": _session_info(),          # 依赖版本快照(铁律6:可复现)
    }

    ks = (10, 25, 50, 100)
    truth_p = Path(a.truth)
    if truth_p.exists():
        print(f"[586] Step 5 用真值 {truth_p.name} 评估")
        tdf = pd.read_csv(truth_p, comment="#", encoding="utf-8")
        truth = set(map(tuple, tdf[["From", "To"]].values))
        ev = {n: evaluate(df, truth, ks) for n, df in methods.items()}
        for n, r in ev.items():
            print(f"       {n:24s} AUPRC={r['AUPRC']:.3f}  " +
                  "  ".join(f"P@{k}={r[f'P@{k}']:.2f}" for k in ks))
        summary["metrics"] = {n: {k: v for k, v in r.items() if not k.startswith("_")}
                              for n, r in ev.items()}
        print("[586] Step 6 出图")
        figs = figures(ev, ks, ASSETS)
        for f in figs:
            print(f"       assets/{f}")
        summary["figures"] = figs
    else:
        print("[586] 无真值文件 → 跳过评估与出图,仅导出边列表")

    print("[586] Step 7 上游官方 PSGRN")
    up = run_psgrn_upstream() if a.run_upstream else {"status": "not requested (--run-upstream)"}
    for k, v in up.items():
        print(f"       {k}: {v}")
    summary["upstream"] = up

    with open(outdir / "586_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=1, ensure_ascii=False, default=str)
    print(f"[586] wrote {outdir / '586_summary.json'}")


if __name__ == "__main__":
    main()
