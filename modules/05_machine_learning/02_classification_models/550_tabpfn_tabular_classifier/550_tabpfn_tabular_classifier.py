# -*- coding: utf-8 -*-
# =============================================================================
# 模块 550 · TabPFN 表格基础模型诊断分类 (TabPFN tabular foundation model)
# -----------------------------------------------------------------------------
# 分类目录 : 05_diagnostic_models
# 语言     : Python (turnkey, CPU 即可)
# 做什么   : 在「小样本 × 表达签名」二分类(case/control)诊断任务上,用表格基础
#            模型 TabPFN 做零训练 in-context 预测,并与两个经典基线【三方对照】,
#            用重复分层交叉验证报告各自 AUROC/AUPRC,客观回答
#            「TabPFN 在你的小数据上是否真的比线性/树模型更好(ΔAUROC 是否真存在)」。
#            出图: 多模型叠加 ROC+PR 曲线、校准(可靠性)曲线、混淆矩阵 heatmap、
#                  TabPFN 置换重要性 lollipop。
#
# ★诚实基线(本模块灵魂 — 对应记忆 reference_dl_ai_strategy「两把屠刀」):
#   2025 年的硬证据是「DL/基础模型常常打不过简单线性基线 / GBDT」。所以预测类
#   任务必须配诚实对照,绝不只报一个好看模型。本模块强制【同管线三方对照】:
#     · TabPFN            (Prior Labs 表格基础模型 v2, Nature 2025; in-context, 零训练)
#     · LASSO-logistic    (L1 稀疏逻辑回归 — 高维小样本经典强基线)
#     · Gradient Boosting (XGBoost 若已装, 否则 sklearn HistGradientBoosting — 树模型代表)
#   三者用同一重复分层 CV、同一【折内】预筛与标准化、同口径报 AUROC/AUPRC,
#   并落盘 ΔAUROC + 配对 t 检验,显式回答增量是否真存在、是否显著。
#
#   ★防数据泄漏(铁律第 7 类): 差异基因预筛 (SelectKBest ANOVA-F top-k) 与
#     标准化(StandardScaler)全部封进 sklearn Pipeline —— 只在每个 CV 训练折上
#     fit,再 transform 验证折,选择器/缩放器永不见验证标签。绝不在划分前全局选基因。
#
# ★TabPFN 权重(诚实说明,实测): tabpfn>=2.x 的【默认 v3 / v2.5 权重是 gated】,
#   首跑需浏览器接受许可(无头环境会卡死)。本模块改用【公开免授权的 v2 权重】:
#   设 TABPFN_MODEL_VERSION=v2 走 GCS 直链下载(_version_has_direct_download_option),
#   无需 token,CPU 实跑。v2 即 TabPFN 原版 Nature 2025 基础模型。首跑联网下载权重
#   (约几十 MB,缓存于 ~/.cache/tabpfn),之后离线可跑。★本模块不含任何"代理/stub"
#   顶替 —— TabPFN 行就是真实基础模型推理;取不到权重则该行如实报错而非伪造数字。
#
# 依赖     : tabpfn(>=2)  scikit-learn  numpy  pandas  matplotlib  (可选 xgboost)
# 运行     : python 550_tabpfn_tabular_classifier.py                 # 合成示例, 零改动跑
#            python 550_tabpfn_tabular_classifier.py --input my.csv  # 换数据
# 输入     : CSV — 首列样本名, 末列 label∈{0,1}, 其余列为基因表达 (样本×基因)
# =============================================================================
from __future__ import annotations

import argparse
import os
import sys
import warnings
from pathlib import Path

# ★必须在 import tabpfn 之前设好:用公开免授权的 v2 权重(GCS 直链, 无需 token)
os.environ.setdefault("TABPFN_MODEL_VERSION", "v2")
os.environ.setdefault("TABPFN_ALLOW_CPU_LARGE_DATASET", "1")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.calibration import calibration_curve
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.inspection import permutation_importance
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    confusion_matrix,
    precision_recall_curve,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import RepeatedStratifiedKFold, cross_val_predict, train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")

# --- 复用顶刊绘图框架 (向上定位 _framework/pubstyle.py) ----------------------
HERE = Path(__file__).resolve().parent
for up in [HERE, *HERE.parents]:
    if (up / "_framework" / "pubstyle.py").exists():
        sys.path.insert(0, str(up / "_framework")); break
try:
    from pubstyle import set_pub_style, save_fig, pal, panel_labels, CMAP_CONT, NATURE_W1, NATURE_W2
except Exception:  # 框架缺失时最小降级,不影响分析
    def set_pub_style(*a, **k): pass
    def save_fig(fig, f, dpi=300):
        from pathlib import Path as _P; _P(f).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(str(f) + ".pdf"); fig.savefig(str(f) + ".png", dpi=dpi)
    def pal(n=None, name="npg"):
        import matplotlib as _m; return list(_m.cm.tab10.colors)
    def panel_labels(*a, **k): pass
    CMAP_CONT = "viridis"; NATURE_W1 = 3.5; NATURE_W2 = 7.0

# --- 全局参数 ------------------------------------------------------------------
SEED = 42
TOPK = 12          # 折内差异基因预筛保留数 (p>n 时先在训练折降维, 防泄漏)
CV_SPLITS = 5      # 重复分层 CV 折数
CV_REPEATS = 4     # 重复次数 (-> 20 个折估计, 稳健 ΔAUROC + 配对检验)
TEST_FRAC = 0.30   # 单独留出测试集(仅用于画曲线/混淆矩阵, 不参与选模型)
DIR_DATA = HERE / "example_data"
DIR_RES = HERE / "results"
DIR_ASSETS = HERE / "assets"
for d in (DIR_DATA, DIR_RES, DIR_ASSETS):
    d.mkdir(parents=True, exist_ok=True)


# =============================================================================
# 1. 合成示例数据 (synthetic, for demo only)
# =============================================================================
def make_synthetic(n=180, n_genes=30, n_signal=6, seed=SEED) -> pd.DataFrame:
    """生成 case/control 表达矩阵 (样本×基因)。前 n_signal 个基因携带真实组间差异
    (强弱递减, 模拟 DEG), 其余为噪声。加一个共享潜变量制造温和共线性, 更像真表达谱。
    synthetic, for demo only —— 仅用于冒烟测试与展示图。"""
    rng = np.random.default_rng(seed)
    y = np.r_[np.zeros(n // 2, int), np.ones(n - n // 2, int)]
    latent = rng.normal(size=(n, 1))
    X = rng.normal(size=(n, n_genes)) + 0.4 * latent
    effects = np.linspace(1.3, 0.45, n_signal)         # 强→弱差异
    for j, eff in enumerate(effects):
        X[y == 1, j] += eff
    genes = [f"SIG{j+1}" for j in range(n_signal)] + \
            [f"NOISE{j+1}" for j in range(n_genes - n_signal)]
    df = pd.DataFrame(X, columns=genes)
    df.insert(0, "sample", [f"S{i+1:03d}" for i in range(n)])
    df["label"] = y
    return df


def load_data(input_path: str | None):
    if input_path:
        df = pd.read_csv(input_path)
        print(f"[load] 用户数据 {input_path}  shape={df.shape}")
    else:
        fp = DIR_DATA / "expr_demo.csv"
        if not fp.exists():
            make_synthetic().to_csv(fp, index=False)
            print(f"[data] 已生成合成示例 (synthetic, demo only) -> {fp}")
        df = pd.read_csv(fp)
        print(f"[load] 合成示例 shape={df.shape}")
    y = df["label"].astype(int).to_numpy()
    feat = df.drop(columns=[c for c in ("sample", "label") if c in df.columns])
    return feat, y


# =============================================================================
# 2. 三方模型工厂 (全部封进 Pipeline: 折内 scaler + 折内预筛 -> 防泄漏)
# =============================================================================
def make_tabpfn():
    """构建真实 TabPFNClassifier (公开 v2 权重, CPU)。返回 (estimator, note)。
    取不到权重时如实抛错 —— 本模块绝不用代理顶替伪造 TabPFN 数字。"""
    from tabpfn import TabPFNClassifier
    # n_estimators=4: in-context ensemble, CPU 上速度/精度折中
    est = TabPFNClassifier(device="cpu", n_estimators=4, random_state=SEED)
    return est, "real TabPFN v2 weights (公开免授权, GCS 直链)"


def make_gbdt():
    """GBDT 基线: XGBoost 若已装则优先(论文常用), 否则 sklearn HistGradientBoosting
    (现代直方图 GBDT, 与 XGB/LightGBM 同档)。返回 (estimator, name)。"""
    try:
        from xgboost import XGBClassifier
        est = XGBClassifier(
            n_estimators=300, max_depth=3, learning_rate=0.05,
            subsample=0.9, colsample_bytree=0.8, eval_metric="logloss",
            random_state=SEED, n_jobs=1, verbosity=0)
        return est, "XGBoost"
    except Exception:
        est = HistGradientBoostingClassifier(
            max_iter=300, max_depth=3, learning_rate=0.05,
            l2_regularization=1.0, random_state=SEED)
        return est, "HistGBDT"


def build_models():
    """返回 dict[name] = (pipeline, color)。
    ★每个模型的预筛(SelectKBest ANOVA-F)与标准化(StandardScaler)都封在 Pipeline
      内 —— CV 时只在训练折 fit, transform 验证折, 选择器/缩放器永不见验证标签。"""
    tab_est, tab_note = make_tabpfn()
    print(f"[tabpfn] {tab_note}")
    gbdt_est, gbdt_name = make_gbdt()
    print(f"[gbdt]   baseline = {gbdt_name}")

    def pipe(final):
        return Pipeline([
            ("scale", StandardScaler()),
            ("prefilter", SelectKBest(score_func=f_classif, k=TOPK)),
            ("clf", final),
        ])

    colors = pal(3, name="npg")
    models = {
        "TabPFN": (pipe(tab_est), colors[0]),
        "LASSO-logistic": (pipe(LogisticRegression(
            penalty="l1", solver="liblinear", C=0.5, max_iter=5000, random_state=SEED)), colors[1]),
        gbdt_name: (pipe(gbdt_est), colors[2]),
    }
    return models, "TabPFN", gbdt_name


# =============================================================================
# 3. 重复分层 CV: 同口径评估三方 + 折级 AUROC 供配对检验
# =============================================================================
def cv_evaluate(models, X, y):
    """对每个模型跑同一组重复分层 CV, 返回:
       fold_auc[name] -> 各折 AUROC 数组(配对, 同折同 split);
       oof[name]      -> out-of-fold 预测概率(供画曲线/校准, 无泄漏)。"""
    cv = RepeatedStratifiedKFold(n_splits=CV_SPLITS, n_repeats=CV_REPEATS, random_state=SEED)
    splits = list(cv.split(X, y))
    fold_auc, fold_ap, oof = {}, {}, {}
    for name, (est, _) in models.items():
        print(f"[cv] {name} ... ({len(splits)} folds)")
        aucs, aps = [], []
        oof_p = np.full(len(y), np.nan)
        for tr, va in splits:
            m = est
            m.fit(X[tr], y[tr])
            p = m.predict_proba(X[va])[:, 1]
            aucs.append(roc_auc_score(y[va], p))
            aps.append(average_precision_score(y[va], p))
            oof_p[va] = p     # 每个样本被多次预测(重复), 末次覆盖即可用于展示
        fold_auc[name] = np.array(aucs)
        fold_ap[name] = np.array(aps)
        oof[name] = oof_p
        print(f"     AUROC={np.mean(aucs):.3f}±{np.std(aucs):.3f}  "
              f"AUPRC={np.mean(aps):.3f}±{np.std(aps):.3f}")
    return fold_auc, fold_ap, oof, splits


def holdout_curves(models, X, y):
    """单独留出测试集, 训练折拟合(含折内预筛)后在测试集出概率 —— 专供画 ROC/PR/
    校准/混淆 的【单一干净评估】(避免重复 CV 的 oof 在曲线上自相关)。"""
    Xtr, Xte, ytr, yte = train_test_split(
        X, y, test_size=TEST_FRAC, random_state=SEED, stratify=y)
    out = {}
    for name, (est, color) in models.items():
        est.fit(Xtr, ytr)
        p = est.predict_proba(Xte)[:, 1]
        out[name] = dict(color=color, proba=p,
                         auroc=roc_auc_score(yte, p),
                         auprc=average_precision_score(yte, p),
                         est=est)
    return out, Xtr, Xte, ytr, yte


# =============================================================================
# 4. 出图 (全部非平凡条形图: 叠加曲线 / lollipop / heatmap)
# =============================================================================
def plot_roc_pr(curves, yte):
    fig, axes = plt.subplots(1, 2, figsize=(NATURE_W2, NATURE_W1 + 0.2))
    axA, axB = axes
    for name, r in curves.items():
        fpr, tpr, _ = roc_curve(yte, r["proba"])
        axA.plot(fpr, tpr, color=r["color"], lw=2, label=f"{name} (AUC={r['auroc']:.3f})")
    axA.plot([0, 1], [0, 1], ls="--", color="0.6", lw=1)
    axA.set_xlabel("False positive rate"); axA.set_ylabel("True positive rate")
    axA.set_title("ROC curves"); axA.set_xlim(-0.02, 1.02); axA.set_ylim(-0.02, 1.02)
    axA.legend(loc="lower right")
    base = float(np.mean(yte))
    for name, r in curves.items():
        prec, rec, _ = precision_recall_curve(yte, r["proba"])
        axB.plot(rec, prec, color=r["color"], lw=2, label=f"{name} (AP={r['auprc']:.3f})")
    axB.axhline(base, ls="--", color="0.6", lw=1, label=f"Prevalence={base:.2f}")
    axB.set_xlabel("Recall"); axB.set_ylabel("Precision")
    axB.set_title("Precision-Recall curves"); axB.set_ylim(-0.02, 1.02)
    axB.legend(loc="lower left")
    panel_labels(axes)
    save_fig(fig, str(DIR_ASSETS / "550_roc_pr")); plt.close(fig)
    print("[fig] 550_roc_pr.png")


def plot_calibration(curves, yte):
    fig, ax = plt.subplots(figsize=(NATURE_W1 + 0.6, NATURE_W1 + 0.4))
    ax.plot([0, 1], [0, 1], ls="--", color="0.5", lw=1, label="Perfectly calibrated")
    nb = min(8, max(3, int(len(yte) / 8)))
    for name, r in curves.items():
        frac_pos, mean_pred = calibration_curve(yte, r["proba"], n_bins=nb, strategy="quantile")
        ax.plot(mean_pred, frac_pos, "o-", color=r["color"], lw=1.8, ms=5, label=name)
    ax.set_xlabel("Mean predicted probability"); ax.set_ylabel("Observed fraction positive")
    ax.set_title("Calibration (reliability) curves")
    ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.02); ax.legend(loc="upper left")
    save_fig(fig, str(DIR_ASSETS / "550_calibration")); plt.close(fig)
    print("[fig] 550_calibration.png")


def plot_confusion(curves, yte):
    n = len(curves)
    fig, axes = plt.subplots(1, n, figsize=(2.5 * n, 2.7))
    for ax, (name, r) in zip(np.ravel(axes), curves.items()):
        pred = (r["proba"] >= 0.5).astype(int)
        cm = confusion_matrix(yte, pred)
        im = ax.imshow(cm, cmap=CMAP_CONT, vmin=0)
        for i in range(2):
            for j in range(2):
                ax.text(j, i, cm[i, j], ha="center", va="center",
                        color="white" if cm[i, j] < cm.max() * 0.6 else "black",
                        fontsize=13, fontweight="bold")
        ax.set_xticks([0, 1]); ax.set_yticks([0, 1])
        ax.set_xticklabels(["Ctrl", "Case"]); ax.set_yticklabels(["Ctrl", "Case"])
        ax.set_xlabel("Predicted"); ax.set_ylabel("True")
        ax.set_title(f"{name}\nAUROC={r['auroc']:.3f}", fontsize=9)
        ax.spines[:].set_visible(True)
    fig.tight_layout()
    save_fig(fig, str(DIR_ASSETS / "550_confusion")); plt.close(fig)
    print("[fig] 550_confusion.png")


def plot_permutation_importance(curves, Xte, yte, feat_names, focus="TabPFN"):
    """对 TabPFN 整 pipeline 做置换重要性(作用在原始基因空间), lollipop 展示。
    红=注入的真信号基因(SIG*), 蓝=噪声 —— 视觉验证模型确实抓到真信号。"""
    r = curves[focus]
    def _auroc_scorer(est, Xv, yv):
        return roc_auc_score(yv, est.predict_proba(Xv)[:, 1])
    pi = permutation_importance(r["est"], Xte, yte, n_repeats=20,
                                random_state=SEED, scoring=_auroc_scorer)
    imp = pd.Series(pi.importances_mean, index=feat_names).sort_values()
    top = imp.tail(15)
    err = pd.Series(pi.importances_std, index=feat_names).loc[top.index]

    fig, ax = plt.subplots(figsize=(NATURE_W1 + 1.4, 0.32 * len(top) + 0.8))
    ypos = np.arange(len(top))
    colors = ["#E64B35" if g.startswith("SIG") else "#3C5488" for g in top.index]
    ax.hlines(ypos, 0, top.values, color="0.7", lw=1.4)
    ax.errorbar(top.values, ypos, xerr=err.values, fmt="none",
                ecolor="0.6", elinewidth=1, capsize=2, zorder=2)
    ax.scatter(top.values, ypos, c=colors, s=55, zorder=3, edgecolor="white", lw=0.6)
    ax.set_yticks(ypos); ax.set_yticklabels(top.index)
    ax.axvline(0, color="0.5", lw=0.8, ls="--")
    ax.set_xlabel("Permutation importance (Δ AUROC)")
    ax.set_title(f"Feature importance · {focus}\n(red = injected signal gene)")
    save_fig(fig, str(DIR_ASSETS / "550_permutation_importance")); plt.close(fig)
    print("[fig] 550_permutation_importance.png")
    return imp


def plot_cv_strip(fold_auc):
    """各模型折级 AUROC 的 strip + 均值点(violin/dot 风, 非条形)——
    顶刊偏好展示分布与不确定性, 而非只画一根均值柱。"""
    names = list(fold_auc.keys())
    cols = pal(len(names), "npg")
    fig, ax = plt.subplots(figsize=(NATURE_W1 + 1.2, NATURE_W1))
    rng = np.random.default_rng(SEED)
    for i, nm in enumerate(names):
        v = fold_auc[nm]
        jitter = rng.uniform(-0.12, 0.12, size=len(v))
        ax.scatter(np.full(len(v), i) + jitter, v, s=26, color=cols[i],
                   alpha=0.55, edgecolor="white", lw=0.4, zorder=2)
        ax.scatter([i], [v.mean()], s=130, color=cols[i], edgecolor="black",
                   lw=1.1, zorder=3)
        ax.errorbar(i, v.mean(), yerr=v.std(), color="black", lw=1.2, capsize=4, zorder=3)
    ax.set_xticks(range(len(names))); ax.set_xticklabels(names, rotation=15, ha="right")
    ax.set_ylabel("Fold AUROC (repeated stratified CV)")
    ax.set_title(f"Cross-validated AUROC ({CV_SPLITS}×{CV_REPEATS} folds)")
    ax.axhline(0.5, ls="--", color="0.6", lw=1)
    save_fig(fig, str(DIR_ASSETS / "550_cv_auroc_strip")); plt.close(fig)
    print("[fig] 550_cv_auroc_strip.png")


# =============================================================================
# 5. 主流程
# =============================================================================
def main():
    ap = argparse.ArgumentParser(description="550 TabPFN vs LASSO/GBDT 诚实三方对照")
    ap.add_argument("--input", default=None,
                    help="CSV: 首列样本名, 末列 label∈{0,1}, 其余列基因表达; 缺省用合成示例")
    args = ap.parse_args()

    set_pub_style(base_size=11, palette="npg")
    print("[step] 加载数据")
    feat, y = load_data(args.input)
    feat_names = list(feat.columns)
    X = feat.to_numpy(float)
    print(f"[step] 样本={X.shape[0]}  候选基因={X.shape[1]}  case比例={y.mean():.2f}  "
          f"(p>n: {X.shape[1] > X.shape[0]})")

    print("[step] 构建三方模型 (折内预筛+标准化, 防泄漏)")
    models, tab_name, gbdt_name = build_models()

    print("[step] 重复分层 CV (同口径三方对照)")
    fold_auc, fold_ap, oof, splits = cv_evaluate(models, X, y)

    print("[step] 留出测试集 -> 曲线/混淆/重要性")
    curves, Xtr, Xte, ytr, yte = holdout_curves(models, X, y)

    print("[step] 绘图")
    plot_cv_strip(fold_auc)
    plot_roc_pr(curves, yte)
    plot_calibration(curves, yte)
    plot_confusion(curves, yte)
    plot_permutation_importance(curves, Xte, yte, feat_names, focus=tab_name)

    # --- CV 汇总表 ---
    print("[step] 写 results/")
    summ = pd.DataFrame({
        "model": list(fold_auc.keys()),
        "CV_AUROC_mean": [fold_auc[k].mean() for k in fold_auc],
        "CV_AUROC_std": [fold_auc[k].std() for k in fold_auc],
        "CV_AUPRC_mean": [fold_ap[k].mean() for k in fold_ap],
        "holdout_AUROC": [curves[k]["auroc"] for k in curves],
    }).sort_values("CV_AUROC_mean", ascending=False)
    summ.to_csv(DIR_RES / "model_comparison.csv", index=False)

    # --- 诚实裁决: TabPFN vs 最强基线, ΔAUROC + 配对 t 检验(同折配对) ---
    base_best = max([k for k in fold_auc if k != tab_name],
                    key=lambda k: fold_auc[k].mean())
    a_tab, a_base = fold_auc[tab_name], fold_auc[base_best]
    delta = a_tab.mean() - a_base.mean()
    tstat, pval = stats.ttest_rel(a_tab, a_base)   # 配对(同 split 同折)
    if delta > 0 and pval < 0.05:
        verdict = f"TabPFN 显著更优 (ΔAUROC={delta:+.3f}, paired t p={pval:.3g})"
    elif delta > 0:
        verdict = f"TabPFN 略优但不显著 (ΔAUROC={delta:+.3f}, p={pval:.3g}) —— 增量证据不足"
    else:
        verdict = f"TabPFN 未胜出 (ΔAUROC={delta:+.3f}, p={pval:.3g}) —— 简单基线已够强"

    lines = [
        "=== 550 TabPFN vs baselines · 诚实三方对照裁决 ===",
        f"TabPFN 路径: 真实 v2 公开权重 (无 proxy, 无 token)",
        f"GBDT 基线  : {gbdt_name}",
        f"数据       : n={X.shape[0]}, p={X.shape[1]}, case={y.mean():.2f}, "
        f"CV={CV_SPLITS}x{CV_REPEATS} 折内预筛 top-{TOPK}",
        "",
        summ.to_string(index=False),
        "",
        f"最强基线   : {base_best}  (CV_AUROC={a_base.mean():.3f})",
        f"★诚实结论  : {verdict}",
        "",
        "★防泄漏    : SelectKBest(ANOVA-F)+StandardScaler 封在 Pipeline 内, 仅训练折拟合。",
        "★配对检验  : 同 split 同折配对 t 检验 (scipy.stats.ttest_rel), 控制折间方差。",
    ]
    (DIR_RES / "verdict.txt").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))

    # --- 版本快照 (铁律6) ---
    import sklearn, matplotlib as mpl, scipy
    vlines = [f"python={sys.version.split()[0]}", f"numpy={np.__version__}",
              f"pandas={pd.__version__}", f"scipy={scipy.__version__}",
              f"scikit-learn={sklearn.__version__}", f"matplotlib={mpl.__version__}"]
    try:
        import tabpfn; vlines.append(f"tabpfn={tabpfn.__version__} (model_version=v2, real weights)")
    except Exception:
        vlines.append("tabpfn=not-installed")
    try:
        import xgboost; vlines.append(f"xgboost={xgboost.__version__}")
    except Exception:
        vlines.append("xgboost=not-installed (GBDT baseline = sklearn HistGradientBoosting)")
    (DIR_RES / "versions.txt").write_text("\n".join(vlines), encoding="utf-8")
    print("[done] 全部完成 -> assets/ + results/")


if __name__ == "__main__":
    main()
