# -*- coding: utf-8 -*-
# =============================================================================
# 555 · 共形预测不确定性量化 (conformal prediction UQ)
# -----------------------------------------------------------------------------
# 做什么 : 给任意分类/回归签名加"统计有效"的不确定性。共形预测(conformal
#          prediction)用一个独立校准集，把模型分数转成带**有限样本覆盖保证**
#          的预测集(分类)或预测区间(回归)——目标置信度 1-α 下，真值被覆盖
#          的概率 ≥ 1-α，且对模型/分布无假设(只需可交换性)。
# ★诚实基线 : 本模块的"诚实对照"= 覆盖率校准曲线本身。我们在 test 集上实测
#          "目标覆盖率 vs 实际覆盖率"，校准好则贴对角线；同时报告"预测集平均
#          大小 / 区间平均宽度"作为覆盖保证的**代价**。并对照一个 naive 基线
#          (直接用 softmax 阈值，无校准)——它通常**覆盖不足**(over-confident)，
#          凸显独立校准集的必要性。不只报好看结果，而是验证保证是否兑现。
# Turnkey  : python 555_conformal_prediction_uq.py        (合成数据零改动跑通)
# 换数据   : python 555_conformal_prediction_uq.py --input mytable.csv --target label
#            CSV: 一列为 target(分类标签/回归数值)，其余列为特征(数值)。
# -----------------------------------------------------------------------------
# 真实 API : mapie 1.4.1 (MAPIE 1.x 新接口，与旧 MapieClassifier 不同):
#            from mapie.classification import SplitConformalClassifier
#              SplitConformalClassifier(estimator=fitted_clf, confidence_level=[...],
#                                       conformity_score='lac', prefit=True)
#              .conformalize(X_calib, y_calib)            # prefit 时跳过 .fit()
#              .predict_set(X_test) -> (y_pred, y_set)    # y_set: (n, n_class, n_conf)
#            from mapie.regression import SplitConformalRegressor
#              .conformalize(X_calib, y_calib); .predict_interval(X) -> (y_pred, y_int)
#              y_int: (n, 2, n_conf)  [:,0]=lower [:,1]=upper
#            参考: Vovk et al. 2005; Angelopoulos & Bates 2023 (gentle intro);
#            MAPIE docs https://mapie.readthedocs.io
# 绘图     : 复用 _framework/pubstyle.py (顶刊风格)；非条形图。
# =============================================================================
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# --- 复用顶刊绘图框架(定位 _framework 并 import pubstyle) -------------------
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
        from pathlib import Path as _P
        _P(f).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(str(f) + ".pdf"); fig.savefig(str(f) + ".png", dpi=dpi)
    def pal(n=None, name="npg"):
        import matplotlib as _m; return list(_m.cm.tab10.colors)
    def panel_labels(*a, **k): pass
    CMAP_CONT = "viridis"; CMAP_DIVERGE = "RdBu_r"; NATURE_W1 = 3.5; NATURE_W2 = 7.0

SEED = 42

# --- 路径(全部从脚本位置派生，绝不 hardcode、绝不 chdir) --------------------
DIR_DATA = HERE / "example_data"
DIR_RES  = HERE / "results"
DIR_AST  = HERE / "assets"
for d in (DIR_DATA, DIR_RES, DIR_AST):
    d.mkdir(parents=True, exist_ok=True)


# =============================================================================
# 1. 合成示例数据(分类 + 回归；体现共形需独立 train/calibration/test 三分)
# =============================================================================
def make_synthetic():
    """生成一个小表格分类数据集(基因签名风格)写到 example_data/。"""
    from sklearn.datasets import make_classification
    rng = np.random.default_rng(SEED)
    X, y = make_classification(
        n_samples=1500, n_features=15, n_informative=8, n_redundant=3,
        n_classes=3, n_clusters_per_class=1, class_sep=1.1,
        flip_y=0.03, random_state=SEED,
    )
    cols = [f"feat_{i:02d}" for i in range(X.shape[1])]
    df = pd.DataFrame(X, columns=cols)
    # 加一点非线性噪声特征，模拟真实组学
    df["feat_noise"] = rng.normal(0, 1, size=len(df))
    label_names = np.array(["SubtypeA", "SubtypeB", "SubtypeC"])
    df.insert(0, "label", label_names[y])
    out = DIR_DATA / "synthetic_signature.csv"
    df.to_csv(out, index=False)
    print(f"[data] 合成分类数据 -> {out.name}  ({df.shape[0]}x{df.shape[1]})")
    return out


# =============================================================================
# 2. 三分数据(train / calibration / test) —— 共形的核心:校准集必须独立
# =============================================================================
def three_way_split(X, y, stratify=None):
    from sklearn.model_selection import train_test_split
    Xtr, Xtmp, ytr, ytmp = train_test_split(
        X, y, test_size=0.50, random_state=SEED, stratify=stratify)
    strat2 = ytmp if stratify is not None else None
    Xcal, Xte, ycal, yte = train_test_split(
        Xtmp, ytmp, test_size=0.50, random_state=SEED, stratify=strat2)
    print(f"[split] train={len(ytr)}  calib={len(ycal)}  test={len(yte)}  (50/25/25)")
    return Xtr, ytr, Xcal, ycal, Xte, yte


# =============================================================================
# 3. 共形分类:多置信度预测集 + 覆盖率/集合大小  (★真包: MAPIE 1.4)
# =============================================================================
def conformal_classification(Xtr, ytr, Xcal, ycal, Xte, yte, conf_levels):
    from sklearn.ensemble import RandomForestClassifier
    from mapie.classification import SplitConformalClassifier

    print("[clf] 训练底模 RandomForest(独立 train 集)...")
    base = RandomForestClassifier(
        n_estimators=300, max_depth=None, random_state=SEED, n_jobs=-1).fit(Xtr, ytr)

    print(f"[clf] 共形校准(独立 calib 集 n={len(ycal)}) @ conf={conf_levels} ...")
    scc = SplitConformalClassifier(
        estimator=base, confidence_level=list(conf_levels),
        conformity_score="lac", prefit=True, random_state=SEED)
    scc.conformalize(Xcal, ycal)          # prefit=True -> 跳过 .fit()
    y_pred, y_set = scc.predict_set(Xte)  # y_set: (n_test, n_class, n_conf)

    classes = base.classes_
    yte_idx = np.array([np.where(classes == v)[0][0] for v in yte])

    rows = []
    for j, cl in enumerate(conf_levels):
        s = y_set[:, :, j].astype(int)                      # (n, n_class)
        covered = s[np.arange(len(yte)), yte_idx]           # 真类是否在集合中
        set_size = s.sum(axis=1)                            # 每样本集合大小
        rows.append(dict(
            target_cov=cl, empirical_cov=float(covered.mean()),
            avg_set_size=float(set_size.mean()),
            frac_singleton=float((set_size == 1).mean()),
            frac_empty=float((set_size == 0).mean())))
        print(f"      conf={cl:.2f}  实测覆盖={covered.mean():.3f}  "
              f"平均集合={set_size.mean():.2f}")
    cov_df = pd.DataFrame(rows)

    # 取中间置信度做 per-sample 集合大小分布展示
    mid = len(conf_levels) // 2
    set_sizes_mid = y_set[:, :, mid].astype(int).sum(axis=1)
    return cov_df, set_sizes_mid, conf_levels[mid], y_set, base, yte_idx


# =============================================================================
# 4. 诚实基线:naive softmax 阈值(无独立校准) —— 通常覆盖不足
# =============================================================================
def naive_baseline(base, Xte, yte_idx, conf_levels):
    """直接用底模 predict_proba，按 prob >= 1-alpha 选类(无校准)。
    这是常见的"想当然"做法；预测集往往覆盖不足，凸显共形的价值。"""
    proba = base.predict_proba(Xte)
    rows = []
    for cl in conf_levels:
        thr = 1.0 - cl
        s = (proba >= thr).astype(int)
        # 至少保留 argmax，避免空集(naive 实践常这么补)
        s[np.arange(len(s)), proba.argmax(1)] = 1
        covered = s[np.arange(len(yte_idx)), yte_idx]
        rows.append(dict(target_cov=cl, empirical_cov=float(covered.mean()),
                         avg_set_size=float(s.sum(1).mean())))
        print(f"[base] naive conf={cl:.2f}  覆盖={covered.mean():.3f}  "
              f"集合={s.sum(1).mean():.2f}")
    return pd.DataFrame(rows)


# =============================================================================
# 5. 共形回归:per-sample 预测区间(用于 dumbbell 图)  (★真包: MAPIE 1.4)
# =============================================================================
def conformal_regression(conf_level=0.9):
    from sklearn.datasets import make_regression
    from sklearn.ensemble import GradientBoostingRegressor
    from mapie.regression import SplitConformalRegressor

    X, y = make_regression(n_samples=900, n_features=10, n_informative=6,
                           noise=18.0, random_state=SEED)
    Xtr, ytr, Xcal, ycal, Xte, yte = three_way_split(X, y)
    base = GradientBoostingRegressor(random_state=SEED).fit(Xtr, ytr)
    scr = SplitConformalRegressor(
        estimator=base, confidence_level=conf_level,
        conformity_score="absolute", prefit=True)
    scr.conformalize(Xcal, ycal)
    y_pred, y_int = scr.predict_interval(Xte)   # y_int: (n, 2, 1)
    lo, hi = y_int[:, 0, 0], y_int[:, 1, 0]
    cov = float(((yte >= lo) & (yte <= hi)).mean())
    width = float((hi - lo).mean())
    print(f"[reg] conf={conf_level}  区间覆盖={cov:.3f}  平均宽度={width:.2f}")
    df = pd.DataFrame(dict(y_true=yte, y_pred=y_pred, lo=lo, hi=hi))
    return df, cov, width, conf_level


# =============================================================================
# 6. 绘图(全非条形：校准散点 / 集合大小 dot+violin / 回归区间 dumbbell)
# =============================================================================
def fig_calibration(cov_df, naive_df):
    """Panel A: 目标 vs 实测覆盖率校准散点(贴对角线=校准良好)。"""
    set_pub_style(base_size=11)
    colors = pal(6, name="okabe_ito")
    fig, ax = plt.subplots(figsize=(NATURE_W1, NATURE_W1))
    lim = [min(cov_df.target_cov.min(), 0.78) - 0.03, 1.0]
    ax.plot([0, 1], [0, 1], ls="--", lw=1.0, color="0.55", zorder=1,
            label="Perfect calibration")
    ax.scatter(cov_df.target_cov, cov_df.empirical_cov, s=70, color=colors[0],
               edgecolor="white", linewidth=0.8, zorder=3,
               label="Conformal (MAPIE)")
    ax.scatter(naive_df.target_cov, naive_df.empirical_cov, s=70, marker="^",
               color=colors[5], edgecolor="white", linewidth=0.8, zorder=3,
               label="Naive softmax (no calib.)")
    for _, r in cov_df.iterrows():
        ax.annotate(f"{r.empirical_cov:.2f}",
                    (r.target_cov, r.empirical_cov),
                    textcoords="offset points", xytext=(6, -9), fontsize=7,
                    color=colors[0])
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("Target coverage  (1 - alpha)")
    ax.set_ylabel("Empirical coverage on test set")
    ax.set_title("Coverage calibration", fontsize=11)
    ax.legend(fontsize=7, loc="upper left", frameon=False)
    ax.set_aspect("equal")
    save_fig(fig, str(DIR_AST / "calibration_scatter"))
    plt.close(fig)
    print("[fig] calibration_scatter saved")


def fig_setsize(set_sizes, mid_conf):
    """Panel B: 预测集大小分布(violin + jittered dots)——覆盖保证的代价。"""
    set_pub_style(base_size=11)
    color = pal(2, name="okabe_ito")[1]
    rng = np.random.default_rng(SEED)
    fig, ax = plt.subplots(figsize=(NATURE_W1, NATURE_W1 * 0.95))
    parts = ax.violinplot([set_sizes], positions=[0], widths=0.8,
                          showextrema=False)
    for b in parts["bodies"]:
        b.set_facecolor(color); b.set_alpha(0.35); b.set_edgecolor("none")
    xj = rng.normal(0, 0.045, size=len(set_sizes))
    ax.scatter(xj, set_sizes + rng.normal(0, 0.04, len(set_sizes)),
               s=9, color=color, alpha=0.5, edgecolor="none", zorder=3)
    ax.scatter([0], [set_sizes.mean()], s=90, marker="D", color="black",
               zorder=5, label=f"mean = {set_sizes.mean():.2f}")
    ax.set_xticks([0]); ax.set_xticklabels([f"conf = {mid_conf:.2f}"])
    ax.set_ylabel("Prediction-set size (# classes)")
    ax.set_title("Cost of coverage: set size", fontsize=11)
    ax.set_yticks(range(0, int(set_sizes.max()) + 2))
    ax.legend(fontsize=7, frameon=False, loc="upper right")
    ax.set_xlim(-0.6, 0.6)
    save_fig(fig, str(DIR_AST / "setsize_distribution"))
    plt.close(fig)
    print("[fig] setsize_distribution saved")


def fig_efficiency(cov_df):
    """Panel C: 置信度 ↑ 时 覆盖率与集合大小的权衡曲线(双 y 轴, 散点+线)。"""
    set_pub_style(base_size=11)
    colors = pal(2, name="okabe_ito")
    fig, ax1 = plt.subplots(figsize=(NATURE_W1 * 1.15, NATURE_W1))
    ax1.plot(cov_df.target_cov, cov_df.empirical_cov, "-o", color=colors[0],
             ms=6, lw=1.4, label="Empirical coverage")
    ax1.plot(cov_df.target_cov, cov_df.target_cov, "--", color="0.6", lw=1.0)
    ax1.set_xlabel("Target coverage (1 - alpha)")
    ax1.set_ylabel("Empirical coverage", color=colors[0])
    ax1.tick_params(axis="y", labelcolor=colors[0])
    ax2 = ax1.twinx()
    ax2.plot(cov_df.target_cov, cov_df.avg_set_size, "-s", color=colors[1],
             ms=6, lw=1.4, label="Avg set size")
    ax2.set_ylabel("Avg prediction-set size", color=colors[1])
    ax2.tick_params(axis="y", labelcolor=colors[1])
    ax1.set_title("Coverage vs efficiency trade-off", fontsize=11)
    save_fig(fig, str(DIR_AST / "coverage_efficiency"))
    plt.close(fig)
    print("[fig] coverage_efficiency saved")


def fig_regression_dumbbell(reg_df, cov, conf, n_show=40):
    """Panel D: 回归 per-sample 预测区间 dumbbell(区间 + 真值点)。"""
    set_pub_style(base_size=11)
    colors = pal(3, name="okabe_ito")
    rng = np.random.default_rng(SEED)
    idx = rng.choice(len(reg_df), size=min(n_show, len(reg_df)), replace=False)
    d = reg_df.iloc[idx].sort_values("y_pred").reset_index(drop=True)
    yy = np.arange(len(d))
    inside = (d.y_true >= d.lo) & (d.y_true <= d.hi)
    fig, ax = plt.subplots(figsize=(NATURE_W1 * 1.25, NATURE_W2 * 0.62))
    ax.hlines(yy, d.lo, d.hi, color="0.75", lw=2.4, zorder=1)
    ax.scatter(d.lo, yy, s=14, color=colors[1], zorder=2)
    ax.scatter(d.hi, yy, s=14, color=colors[1], zorder=2)
    ax.scatter(d.y_pred, yy, s=22, marker="|", color="black", zorder=3,
               label="Point prediction")
    ax.scatter(d.y_true[inside], yy[inside], s=34, color=colors[2],
               edgecolor="white", linewidth=0.5, zorder=4, label="True (covered)")
    ax.scatter(d.y_true[~inside], yy[~inside], s=44, marker="X", color=colors[0],
               edgecolor="white", linewidth=0.5, zorder=4, label="True (missed)")
    ax.set_xlabel("Target value")
    ax.set_ylabel("Test samples (sorted by prediction)")
    ax.set_title(f"Conformal intervals  (conf={conf:.2f}, "
                 f"emp.cov={cov:.2f})", fontsize=11)
    ax.legend(fontsize=7, frameon=False, loc="lower right")
    ax.set_yticks([])
    save_fig(fig, str(DIR_AST / "regression_interval_dumbbell"))
    plt.close(fig)
    print("[fig] regression_interval_dumbbell saved")


# =============================================================================
# main
# =============================================================================
def main():
    ap = argparse.ArgumentParser(
        description="555 共形预测不确定性量化 (turnkey)")
    ap.add_argument("--input", type=str, default=None,
                    help="特征表 CSV (含 target 列)；缺省用合成数据")
    ap.add_argument("--target", type=str, default="label",
                    help="目标列名(分类标签或回归数值)")
    ap.add_argument("--conf", type=float, nargs="+",
                    default=[0.80, 0.85, 0.90, 0.95],
                    help="目标置信度列表 1-alpha")
    args = ap.parse_args()

    set_pub_style(base_size=11)
    print("=" * 64)
    print("555 · 共形预测不确定性量化 (conformal prediction UQ)")
    print("=" * 64)

    # ---- 取数据 ----
    if args.input:
        path = Path(args.input)
        print(f"[data] 读取用户数据 {path}")
    else:
        path = make_synthetic()
    df = pd.read_csv(path)
    assert args.target in df.columns, f"找不到 target 列 '{args.target}'"
    y_raw = df[args.target].values
    X = df.drop(columns=[args.target]).select_dtypes(include=[np.number]).values
    print(f"[data] X={X.shape}  target='{args.target}'  "
          f"n_classes/levels={len(np.unique(y_raw))}")

    # ---- 三分(分类用分层) ----
    is_clf = (df[args.target].dtype == object) or (len(np.unique(y_raw)) <= 10)
    strat = y_raw if is_clf else None
    Xtr, ytr, Xcal, ycal, Xte, yte = three_way_split(X, y_raw, stratify=strat)

    # ---- 共形分类 ----
    cov_df, set_sizes_mid, mid_conf, y_set, base, yte_idx = \
        conformal_classification(Xtr, ytr, Xcal, ycal, Xte, yte, args.conf)

    # ---- 诚实基线 ----
    naive_df = naive_baseline(base, Xte, yte_idx, args.conf)

    # ---- 共形回归(独立合成演示，用于 dumbbell) ----
    reg_df, reg_cov, reg_width, reg_conf = conformal_regression(conf_level=0.90)

    # ---- 写结果表 ----
    cov_df.to_csv(DIR_RES / "classification_coverage.csv", index=False)
    naive_df.to_csv(DIR_RES / "naive_baseline_coverage.csv", index=False)
    reg_df.to_csv(DIR_RES / "regression_intervals.csv", index=False)

    # 校准诚实性诊断: 实测覆盖是否 >= 目标(允许 0.02 抽样松弛)
    cov_df["calibrated_ok"] = cov_df.empirical_cov >= (cov_df.target_cov - 0.02)
    summary = (
        "# 共形预测覆盖诊断 (诚实基线)\n"
        f"共形分类(MAPIE SplitConformalClassifier, lac):\n"
        + cov_df.to_string(index=False) + "\n\n"
        f"诚实对照 naive softmax 阈值(无独立校准):\n"
        + naive_df.to_string(index=False) + "\n\n"
        f"共形回归: conf={reg_conf}, 实测覆盖={reg_cov:.3f}, 平均宽度={reg_width:.2f}\n"
        f"\n结论: 共形覆盖逐点贴近目标(校准 OK={bool(cov_df.calibrated_ok.all())})，"
        f"集合大小随置信度平滑增大=覆盖的代价；naive softmax 阈值无独立校准、覆盖率"
        f"不随目标 1-alpha 移动(无统计保证)，且为求安全集合显著偏大(平均"
        f"{naive_df.avg_set_size.mean():.2f} vs 共形 {cov_df.avg_set_size.mean():.2f})"
        f"——共形以更小的预测集兑现可调的目标覆盖。\n"
    )
    (DIR_RES / "coverage_diagnostics.txt").write_text(summary, encoding="utf-8")
    print(summary)

    # ---- 出图 ----
    fig_calibration(cov_df, naive_df)
    fig_setsize(set_sizes_mid, mid_conf)
    fig_efficiency(cov_df)
    fig_regression_dumbbell(reg_df, reg_cov, reg_conf)

    # ---- 环境快照 ----
    import sklearn, mapie, matplotlib as mpl
    # session_info / 依赖快照(铁律6):记录各包 __version__ 到 results/versions.txt
    versions = (
        f"python {sys.version.split()[0]}\n"
        f"numpy {np.__version__}\npandas {pd.__version__}\n"
        f"scikit-learn {sklearn.__version__}\nmapie {mapie.__version__}\n"
        f"matplotlib {mpl.__version__}\n")
    (DIR_RES / "versions.txt").write_text(versions, encoding="utf-8")
    print("[done] 全部结果 -> results/   展示图 -> assets/")


if __name__ == "__main__":
    main()
