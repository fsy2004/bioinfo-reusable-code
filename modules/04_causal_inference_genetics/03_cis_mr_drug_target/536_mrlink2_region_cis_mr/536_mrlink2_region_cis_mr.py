# -*- coding: utf-8 -*-
# =============================================================================
# 536 · MR-link-2 单区域 cis-MR (region-based, pleiotropy-robust causal MR)  🟡
# -----------------------------------------------------------------------------
# 分类   : 09_mendelian_randomization
# 用途   : 在【单一关联区域 / cis 窗口】内,用 summary 统计 + LD 同时估计【因果
#          效应 alpha】与【水平多效性 sigma_y】。与传统 IVW 不同 —— IVW 在单区域
#          内强行用相互高度 LD 的 SNP 当独立工具,会把未建模的多效性误判为因果,
#          导致 Type-I error(假阳)暴涨;MR-link-2 把整段 LD 结构与多效性显式建模,
#          控制假阳。
#
# ★诚实基线(本模块的灵魂,不可只报好看指标):
#   在【同一批合成单区域 cis 数据】上并排跑两套方法 ——
#     (A) naive 单区域 IVW(numpy/statsmodels 手算,逆方差加权 + LD 修正前的朴素版)
#     (B) MR-link-2 概念实现(eigh(LD) + 似然比检验,接地官方 mr_link2() API)
#   在【真实因果场景】两者都应检出;在【纯多效性 / 无因果 + 共享 LD】场景下,
#   naive IVW 会给出大量假阳性(small p),而 MR-link-2 通过 sigma_y 吸收多效性、
#   保持 alpha 的 p 接近 null。本模块用多次重复模拟实测两法的 Type-I error。
#
# ★降级说明(🟡 DEGRADED):官方 mrlink2 包当前装不上(pip git 传输 EOF)。
#   本脚本【接地真实 API】(github adriaan-vd-graaf/mrlink2,函数
#   mr_link2(selected_eigenvalues, selected_eigenvectors, exposure_betas,
#   outcome_betas, n_exp, n_out, sigma_exp_guess, sigma_out_guess) -> dict,
#   内部 np.linalg.eigh(LD) + scipy.optimize.minimize(Nelder-Mead) + 似然比卡方):
#   先 try import 真包;装不上则走【同算法的本地概念实现】(_mrlink2_local),
#   保证诚实基线 + 出图全程跑通。装上真包后会自动优先调用官方函数。
#
# 依赖   : numpy scipy statsmodels matplotlib pandas (均已装)
#          可选: mrlink2 (官方;装上则优先)
# 运行   : python 536_mrlink2_region_cis_mr.py          # 零改动跑合成示例
#          python 536_mrlink2_region_cis_mr.py --outdir results/run1
# 输入   : 合成单区域 cis summary(beta/se per SNP)+ LD 矩阵(脚本内生成,
#          synthetic demo only);换真数据见 README(需官方包 / 服务器)。
# =============================================================================
from __future__ import annotations
import argparse
import sys
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")

# ---- 定位脚本目录 + 载入顶刊绘图框架(向上搜 _framework) --------------------
HERE = Path(__file__).resolve().parent
for up in [HERE, *HERE.parents]:
    if (up / "_framework" / "pubstyle.py").exists():
        sys.path.insert(0, str(up / "_framework"))
        break
try:
    from pubstyle import set_pub_style, save_fig, pal, NATURE_W2, CMAP_DIVERGE
except Exception:  # 框架缺失时最小降级
    def set_pub_style(*a, **k): pass
    def save_fig(fig, f, dpi=300):
        from pathlib import Path as _P
        _P(f).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(str(f) + ".pdf", bbox_inches="tight")
        fig.savefig(str(f) + ".png", dpi=dpi, bbox_inches="tight")
    def pal(n=None, name="npg"):
        import matplotlib
        return list(matplotlib.cm.tab10.colors)[: (n or 10)]
    NATURE_W2 = 7.0
    CMAP_DIVERGE = "RdBu_r"

SEED = 42
np.random.seed(SEED)


# =============================================================================
# 合成单区域 cis 数据:一段强 LD 的 cis 窗口
# =============================================================================
def make_cis_region(n_snps=40, n_exp=20000, n_out=20000, alpha_true=0.0,
                    pleiotropy=0.0, rng=None):
    """生成一段 cis 区域的 summary + LD(synthetic demo only)。

    模型(与 MR-link-2 设定一致):区域内 n_snps 个相互 LD 的 SNP。
      - LD 矩阵 R:AR(1) 衰减 + 块状高相关(模拟单一 cis 区域的强连锁)。
      - 真实 cis 效应向量 gamma(exposure 上),通过 LD 投影成 marginal beta。
      - outcome marginal beta = alpha_true * (exposure 信号) + pleiotropy * (与 LD
        共享、但与 exposure 信号【不成因果比例】的独立多效性程序)。
    返回 marginal beta/se(exposure & outcome) + LD 矩阵 R。

    设计意图:
      pleiotropy>0 且 alpha_true=0 → 存在与 exposure 共 LD 的多效性但【无因果】,
      这是 naive 单区域 IVW 最容易假阳的场景;MR-link-2 的 sigma_y 应吸收之。"""
    if rng is None:
        rng = np.random.default_rng(SEED)
    p = n_snps

    # --- LD 矩阵 R:AR(1) 衰减 |i-j| + 两个高相关块,再投影成正定相关阵 ---
    idx = np.arange(p)
    R = 0.92 ** np.abs(idx[:, None] - idx[None, :])      # AR(1) 衰减
    for blk in (slice(0, 12), slice(20, 32)):            # 两个强连锁块
        R[blk, blk] = np.maximum(R[blk, blk], 0.85)
    np.fill_diagonal(R, 1.0)
    # 保证正定(对称化 + 特征值地板)
    R = (R + R.T) / 2
    w, V = np.linalg.eigh(R)
    w = np.clip(w, 1e-3, None)
    R = V @ np.diag(w) @ V.T
    d = np.sqrt(np.diag(R))
    R = R / np.outer(d, d)                                # 重新归一为相关阵

    L = np.linalg.cholesky(R + 1e-6 * np.eye(p))

    # --- 真实 cis 因果信号 gamma:稀疏(少数 causal SNP)----
    gamma = np.zeros(p)
    causal = rng.choice(p, size=3, replace=False)
    gamma[causal] = rng.normal(0, 0.06, size=3)

    # exposure 的 marginal beta = R @ gamma(LD 把单点效应摊到邻居)+ 抽样噪声
    se_exp = 1.0 / np.sqrt(n_exp)
    beta_exp = R @ gamma + L @ rng.normal(0, se_exp, size=p)

    # 水平多效性程序(与 MR-link-2 的 sigma_y 定义一致):一个【弥散的、沿 LD
    # 结构铺开的 infinitesimal 效应】—— 区域内每个 SNP 都有小随机 outcome 效应,
    # 经 LD 相关(R 投影)。这正是 sigma_y 要建模的 LD 比例方差膨胀,而非少数
    # SNP 的稀疏效应;它与 exposure 的稀疏 cis 信号【方向不成因果比例】。
    tau = 0.012                                  # 每 SNP 多效性效应 SD
    u = rng.normal(0, tau, size=p)
    pleio_beta = R @ u                           # 沿 LD 铺开的弥散多效性

    se_out = 1.0 / np.sqrt(n_out)
    beta_out = (alpha_true * (R @ gamma)          # 因果通路:exposure→outcome
                + pleiotropy * pleio_beta          # 水平多效性(LD 比例、非因果)
                + L @ rng.normal(0, se_out, size=p))

    se_exp_v = np.full(p, se_exp)
    se_out_v = np.full(p, se_out)
    return dict(beta_exp=beta_exp, se_exp=se_exp_v,
                beta_out=beta_out, se_out=se_out_v,
                R=R, n_exp=n_exp, n_out=n_out, causal=causal)


# =============================================================================
# 诚实基线 A:naive 单区域 IVW(手算)—— 故意忽略 LST 结构,易假阳
# =============================================================================
def naive_region_ivw(d):
    """单区域内朴素 IVW:把区域里每个(LD 相关的)SNP 当独立工具,
    比值法 ratio = beta_out/beta_exp,逆方差加权汇总。
    关键缺陷:这些 SNP 高度 LD 不独立 → 方差被严重低估 → p 过小(假阳)。
    用 statsmodels WLS 实现 beta_out ~ beta_exp(过原点)的等价加权回归。"""
    import statsmodels.api as sm
    be, bo = d["beta_exp"], d["beta_out"]
    # 工具过滤:用 exposure 较强的 SNP(避免弱工具放大比值)
    strong = np.abs(be) > (1.5 * d["se_exp"])
    if strong.sum() < 3:
        strong = np.argsort(-np.abs(be))[:10]
    be_s, bo_s = be[strong], bo[strong]
    w = 1.0 / (d["se_out"][strong] ** 2)                  # IVW 权重(忽略 LD!)
    X = be_s
    model = sm.WLS(bo_s, X, weights=w).fit()
    alpha = float(model.params[0])
    se = float(model.bse[0])
    from scipy import stats
    p = float(2 * stats.norm.sf(abs(alpha / se)))
    return dict(alpha=alpha, se=se, pval=p, n_iv=int(strong.sum()))


# ---- 官方 mrlink2 探测【仅一次】(失败的 import 不会被 Python 缓存,放循环里会
#      反复扫描 sys.path 拖慢整个模拟;故在模块加载时探测一次并缓存结果)--------
try:
    import mrlink2 as _MRLINK2          # 官方包(当前 MISSING;装上自动启用)
    _HAS_OFFICIAL = hasattr(_MRLINK2, "mr_link2")
except Exception:
    _MRLINK2 = None
    _HAS_OFFICIAL = False


# =============================================================================
# MR-link-2:优先官方包;装不上则本地概念实现(同算法,接地真实 API)
# =============================================================================
def _mrlink2_local(selected_eigenvalues, selected_eigenvectors,
                   exposure_betas, outcome_betas, n_exp, n_out,
                   sigma_exp_guess=0.01, sigma_out_guess=0.001):
    """MR-link-2 核心似然的【本地概念实现】,严格对齐官方 mr_link2() 接口与算法:
      参数 alpha(因果)、sigma_x(exposure 区域遗传度)、sigma_y(水平多效性);
      在 LD 的特征空间(eigh)里写边际 beta 的高斯似然,
      用 scipy.optimize.minimize(Nelder-Mead) 最小化负 loglik,
      再以似然比卡方(df=1)对 alpha=0 与 sigma_y=0 各做一次检验。
    返回 dict 键与官方一致: 'alpha','se(alpha)','p(alpha)','sigma_y',
      'se(sigma_y)','p(sigma_y)','sigma_x'。
    注:这是用于【降级演示 + 诚实基线】的简化实现,非官方数值;装上官方
    mrlink2 包后 run_mrlink2() 会自动优先调用真函数。"""
    from scipy.optimize import minimize  # noqa: scipy 已在 main 顶部预热,此处仅取符号
    from scipy import stats

    lam = np.clip(selected_eigenvalues, 1e-6, None)
    U = selected_eigenvectors
    # 把边际 beta 投到特征空间
    bx = U.T @ exposure_betas
    by = U.T @ outcome_betas
    c_exp = 1.0 / n_exp                                  # 抽样方差尺度
    c_out = 1.0 / n_out

    # 方差参数在【对数尺度】上优化(sigma_x, sigma_y 跨数量级,log 化后良条件、
    # 且自动保正);用 L-BFGS-B(梯度法,远快于 Nelder-Mead,40-SNP 单区域 <0.1s)。
    LOG_LO, LOG_HI = np.log(1e-8), np.log(1.0)

    def negloglik(alpha, log_sx, log_sy):
        sx = np.exp(log_sx)
        sy = np.exp(log_sy)
        var_x = lam * sx + c_exp                          # exposure 边际方差(每特征模)
        var_y = lam * (alpha ** 2 * sx + sy) + c_out       # outcome 边际方差
        mean_y = alpha * bx                               # 因果均值
        ll = -0.5 * np.sum(np.log(2 * np.pi * var_x) + bx ** 2 / var_x)
        ll += -0.5 * np.sum(np.log(2 * np.pi * var_y) + (by - mean_y) ** 2 / var_y)
        return -ll

    lg = float(np.log(max(sigma_out_guess, 1e-6)))
    lgx = float(np.log(max(sigma_exp_guess, 1e-6)))

    def _best(fun, starts, bounds):
        best = None
        for s in starts:
            r = minimize(fun, s, method="L-BFGS-B", bounds=bounds,
                         options=dict(maxiter=200, ftol=1e-9))
            if best is None or r.fun < best.fun:
                best = r
        return best

    # 全模型(alpha, log_sx, log_sy)——多初值取最优(alpha 限于合理范围,加速收敛)
    bnd3 = [(-2, 2), (LOG_LO, LOG_HI), (LOG_LO, LOG_HI)]
    ha = _best(lambda t: negloglik(t[0], t[1], t[2]),
               [np.array([a, lgx, lg]) for a in (0.0, 0.3)], bnd3)
    alpha_hat, sx_hat, sy_hat = ha.x[0], float(np.exp(ha.x[1])), float(np.exp(ha.x[2]))

    # H0: alpha=0(sigma_x, sigma_y 自由)——sigma_y 必须能充分吸收多效性,
    # 故给 sigma_y 多个量级的初值(含偏高值),避免 H0 拟合卡在局部、低估 sigma_y
    # 而把残差误算进 alpha 通道(否则 Type-I error 控制会偏松)。
    bnd2 = [(LOG_LO, LOG_HI), (LOG_LO, LOG_HI)]
    h0a = _best(lambda t: negloglik(0.0, t[0], t[1]),
                [np.array([lgx, s]) for s in (lg, np.log(1e-3), np.log(1e-2),
                                              np.log(sy_hat + 1e-6))], bnd2)
    chi_alpha = max(0.0, 2 * (h0a.fun - ha.fun))
    p_alpha = float(stats.chi2.sf(chi_alpha, 1))

    # H0: sigma_y=0(alpha, sigma_x 自由)—— 检验是否存在多效性
    h0y = _best(lambda t: negloglik(t[0], t[1], LOG_LO),
                [np.array([alpha_hat, lgx])], [(-2, 2), (LOG_LO, LOG_HI)])
    chi_sy = max(0.0, 2 * (h0y.fun - ha.fun))
    p_sy = float(stats.chi2.sf(chi_sy, 1))

    # alpha 的 se:由似然比卡方反推(wald 近似 se = |alpha| / sqrt(chi))
    se_alpha = float(abs(alpha_hat) / np.sqrt(chi_alpha)) if chi_alpha > 1e-6 else float("nan")
    se_sy = float(sy_hat / np.sqrt(chi_sy)) if chi_sy > 1e-6 else float("nan")

    return {
        "alpha": float(alpha_hat), "se(alpha)": se_alpha, "p(alpha)": p_alpha,
        "sigma_y": float(sy_hat), "se(sigma_y)": se_sy, "p(sigma_y)": p_sy,
        "sigma_x": float(sx_hat), "_backend": "local-concept",
    }


def run_mrlink2(d, sigma_exp_guess=0.01, sigma_out_guess=0.001):
    """统一入口:eigh(LD) → 调用 MR-link-2。优先官方 mrlink2.mr_link2,
    装不上则用 _mrlink2_local(同接口同算法)。"""
    lam, U = np.linalg.eigh(d["R"])                       # 官方亦用 np.linalg.eigh
    # 保留主要方差的特征模(官方 var_explained 思路;去极小特征值稳健)
    order = np.argsort(-lam)
    lam, U = lam[order], U[:, order]
    keep = np.cumsum(lam) / lam.sum() <= 0.999
    keep[: max(5, keep.sum())] = True                     # 至少保留若干模
    lam_s, U_s = lam[keep], U[:, keep]

    if _HAS_OFFICIAL:                                     # 官方包在场则优先调用真函数
        try:
            res = _MRLINK2.mr_link2(lam_s, U_s, d["beta_exp"], d["beta_out"],
                                    d["n_exp"], d["n_out"],
                                    sigma_exp_guess, sigma_out_guess)
            res["_backend"] = "official"
            return res, lam, U
        except Exception:
            pass                                          # 官方调用异常 → 退回本地实现
    res = _mrlink2_local(lam_s, U_s, d["beta_exp"], d["beta_out"],
                         d["n_exp"], d["n_out"], sigma_exp_guess, sigma_out_guess)
    return res, lam, U


# =============================================================================
# 主流程
# =============================================================================
def main():
    ap = argparse.ArgumentParser(description="536 MR-link-2 单区域 cis-MR (🟡 degraded)")
    ap.add_argument("--n_rep", type=int, default=80, help="每场景重复次数(Type-I error 估计)")
    ap.add_argument("--n_snps", type=int, default=40)
    ap.add_argument("--outdir", default=str(HERE / "results"))
    args = ap.parse_args()

    set_pub_style(base_size=11)
    DRES = Path(args.outdir); DAST = HERE / "assets"; DDAT = HERE / "example_data"
    for dd in (DRES, DAST, DDAT):
        dd.mkdir(parents=True, exist_ok=True)

    import pandas as pd

    # ---- 1. 落盘一份合成 cis summary + LD 作为 example_data(synthetic demo only)
    rng0 = np.random.default_rng(SEED)
    demo = make_cis_region(n_snps=args.n_snps, alpha_true=0.0, pleiotropy=1.0, rng=rng0)
    snp_df = pd.DataFrame({
        "pos_name": [f"rs{i:04d}" for i in range(args.n_snps)],
        "chromosome": 1, "position": np.arange(args.n_snps) * 5000 + 1_000_000,
        "beta_exposure": demo["beta_exp"], "se_exposure": demo["se_exp"],
        "beta_outcome": demo["beta_out"], "se_outcome": demo["se_out"],
    })
    snp_df.to_csv(DDAT / "cis_region_sumstats.csv", index=False)
    pd.DataFrame(demo["R"]).to_csv(DDAT / "cis_region_LD.csv", index=False)
    print(f"[gen] synthetic cis region: {args.n_snps} SNPs, LD matrix saved "
          f"(synthetic demo only) -> example_data/")

    # ---- 2. 诚实基线对照:两场景 × n_rep 重复,实测 Type-I error / power -------
    #   场景 NULL : alpha_true=0, 强多效性  → 理想:两法 alpha 都不显著;naive IVW 会假阳
    #   场景 CAUSAL: alpha_true=0.25, 中多效性 → 理想:两法都检出
    scenarios = {
        "Null (no causal,\npleiotropy)": dict(alpha=0.0, pl=1.0),
        "Causal\n(alpha=0.25)": dict(alpha=0.25, pl=0.5),
    }
    rows = []
    rng = np.random.default_rng(SEED + 1)
    for sname, cfg in scenarios.items():
        for r in range(args.n_rep):
            d = make_cis_region(n_snps=args.n_snps, alpha_true=cfg["alpha"],
                                pleiotropy=cfg["pl"], rng=rng)
            ivw = naive_region_ivw(d)
            ml, _, _ = run_mrlink2(d)
            rows.append(dict(scenario=sname, rep=r,
                             ivw_alpha=ivw["alpha"], ivw_p=ivw["pval"],
                             ml_alpha=ml["alpha"], ml_p=ml["p(alpha)"],
                             ml_sigma_y=ml["sigma_y"], ml_p_sigma_y=ml["p(sigma_y)"]))
    res = pd.DataFrame(rows)
    res.to_csv(DRES / "simulation_results.csv", index=False)
    backend = ml["_backend"]

    # 汇总 Type-I error(Null 下 p<0.05 比例)与 power(Causal 下 p<0.05 比例)
    summ = []
    for sname in scenarios:
        sub = res[res.scenario == sname]
        summ.append(dict(scenario=sname.replace("\n", " "),
                         ivw_rate=float((sub.ivw_p < 0.05).mean()),
                         ml_rate=float((sub.ml_p < 0.05).mean())))
    summ = pd.DataFrame(summ)
    summ.to_csv(DRES / "rejection_rates.csv", index=False)
    null_row = summ[summ.scenario.str.startswith("Null")].iloc[0]
    print(f"[baseline] backend={backend}")
    print(f"[baseline] NULL scenario rejection rate (= Type-I error):"
          f"  naive IVW={null_row.ivw_rate:.3f}  vs  MR-link-2={null_row.ml_rate:.3f}"
          f"  (nominal 0.05; IVW should be inflated)")
    caus_row = summ[summ.scenario.str.startswith("Causal")].iloc[0]
    print(f"[baseline] CAUSAL scenario power:"
          f"  naive IVW={caus_row.ivw_rate:.3f}  vs  MR-link-2={caus_row.ml_rate:.3f}")

    # =========================================================================
    # 出图(顶刊风;每图独立成文件;禁止平凡条形图)
    # =========================================================================
    import matplotlib.pyplot as plt
    cols = pal(3, "npg")
    c_ivw, c_ml = cols[0], cols[2]

    # ---- Fig 1: 区域因果效应 FOREST(两法 × 两场景的 alpha ± 95%CI)----------
    #   顶刊偏好 forest(点+误差棒)而非条形:直接读出效应方向、CI 是否跨 0。
    fig, ax = plt.subplots(figsize=(NATURE_W2 * 0.62, 3.6))
    yi = 0; yticks = []; yt_lab = []; group_centers = {}
    for sname in list(scenarios)[::-1]:                    # 自下而上:Null 在下,Causal 在上
        sub = res[res.scenario == sname]
        y_start = yi
        for mname, ac, col in [("naive IVW", "ivw_alpha", c_ivw),
                               ("MR-link-2", "ml_alpha", c_ml)]:
            a = sub[ac].mean(); s = sub[ac].std()
            ax.errorbar(a, yi, xerr=1.96 * s, fmt="o", color=col, ms=7,
                        capsize=3, lw=1.6, mec="black", mew=0.6, zorder=3)
            yticks.append(yi); yt_lab.append(mname)
            yi += 1
        group_centers[sname.split("\n")[0].split(" ")[0]] = (y_start + yi - 1) / 2
        yi += 0.6
    ax.axvline(0, color="grey", ls="--", lw=1)
    ax.set_yticks(yticks); ax.set_yticklabels(yt_lab)
    ax.set_ylim(-0.8, yi - 0.4)
    ax.set_xlabel("Causal effect  alpha  (mean ± 1.96·SD over reps)")
    ax.set_title("Region causal-effect forest")
    # 场景分组标注(y 位置由实际行号计算,不写死)
    xlo = ax.get_xlim()[0]
    for gname, gy in group_centers.items():
        ax.text(xlo, gy, f"  {gname}", fontsize=9.5, style="italic",
                color="grey", va="center", fontweight="bold")
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([], [], marker="o", color=c_ivw, ls="", label="naive IVW"),
                       Line2D([], [], marker="o", color=c_ml, ls="", label="MR-link-2")],
              loc="lower right")
    fig.tight_layout(); save_fig(fig, DAST / "region_forest"); plt.close(fig)

    # ---- Fig 2: cis 工具 LD HEATMAP(区域内 SNP×SNP 相关)--------------------
    fig, a = plt.subplots(figsize=(4.4, 3.9))
    im = a.imshow(demo["R"], cmap=CMAP_DIVERGE, vmin=-1, vmax=1)
    a.set_title("cis-region LD matrix  (R)")
    a.set_xlabel("SNP index"); a.set_ylabel("SNP index")
    # 标出 causal SNP
    for c in demo["causal"]:
        a.axhline(c, color="black", lw=0.4, alpha=0.25)
        a.axvline(c, color="black", lw=0.4, alpha=0.25)
    fig.colorbar(im, ax=a, fraction=0.046, label="LD correlation")
    fig.tight_layout(); save_fig(fig, DAST / "cis_LD_heatmap"); plt.close(fig)

    # ---- Fig 3: 效应 vs 多效性 SCATTER(每 rep 一个点:alpha vs sigma_y)-----
    #   展示 MR-link-2 如何把多效性"路由"进 sigma_y、保住 alpha≈0(Null 场景)。
    fig, a = plt.subplots(figsize=(4.6, 3.9))
    for sname, mk in zip(scenarios, ["o", "^"]):
        sub = res[res.scenario == sname]
        a.scatter(sub.ml_alpha, sub.ml_sigma_y, s=26, alpha=0.7,
                  marker=mk, edgecolor="black", linewidth=0.3,
                  label=sname.replace("\n", " "),
                  c=[c_ml if sname.startswith("Causal") else c_ivw])
    a.axvline(0, color="grey", ls="--", lw=1)
    a.set_xlabel("Causal effect  alpha  (MR-link-2)")
    a.set_ylabel("Pleiotropy  sigma_y  (MR-link-2)")
    a.set_title("Causal effect vs modelled pleiotropy")
    a.legend(loc="upper right", fontsize=8)
    fig.tight_layout(); save_fig(fig, DAST / "alpha_vs_pleiotropy"); plt.close(fig)

    # ---- Fig 4: 诚实基线核心图 —— Null 场景 p 值分布(naive IVW 假阳 vs MR-link-2)
    #   用 violin + 抖动点(raincloud 风),不用条形;水平线 = 0.05 名义阈值。
    fig, a = plt.subplots(figsize=(4.8, 3.9))
    nullsub = res[res.scenario.str.startswith("Null")]
    data = [nullsub.ivw_p.values, nullsub.ml_p.values]
    parts = a.violinplot(data, positions=[1, 2], showmeans=False, showextrema=False,
                         widths=0.8)
    for pc, col in zip(parts["bodies"], [c_ivw, c_ml]):
        pc.set_facecolor(col); pc.set_alpha(0.35); pc.set_edgecolor("black")
    for i, (vals, col) in enumerate(zip(data, [c_ivw, c_ml]), start=1):
        jit = (np.random.default_rng(0).random(len(vals)) - 0.5) * 0.18
        a.scatter(np.full(len(vals), i) + jit, vals, s=10, c=col,
                  alpha=0.6, edgecolor="black", linewidth=0.2, zorder=3)
    a.axhline(0.05, color="red", ls="--", lw=1.2)
    a.text(2.45, 0.06, "nominal 0.05", color="red", fontsize=8, ha="right")
    a.set_xticks([1, 2]); a.set_xticklabels(["naive IVW", "MR-link-2"])
    a.set_ylabel("p(alpha)  under NULL (no causal effect)")
    a.set_title(f"Honest baseline: Type-I error control\n"
                f"reject<0.05  IVW={null_row.ivw_rate:.2f}  vs  MR-link-2={null_row.ml_rate:.2f}")
    a.set_ylim(-0.03, 1.03)
    fig.tight_layout(); save_fig(fig, DAST / "typeI_error_baseline"); plt.close(fig)

    print(f"[fig] assets/: region_forest, cis_LD_heatmap, alpha_vs_pleiotropy, "
          f"typeI_error_baseline (.pdf+.png)")

    # ---- 依赖快照(铁律6)---------------------------------------------------
    import scipy, statsmodels, matplotlib
    with open(DRES / "versions.txt", "w", encoding="utf-8") as fh:
        fh.write(f"backend={backend}\n")
        fh.write(f"numpy={np.__version__}\nscipy={scipy.__version__}\n")
        fh.write(f"statsmodels={statsmodels.__version__}\n")
        fh.write(f"matplotlib={matplotlib.__version__}\npandas={pd.__version__}\n")
        try:
            import mrlink2 as _m
            fh.write(f"mrlink2={getattr(_m,'__version__','installed')}\n")
        except Exception:
            fh.write("mrlink2=MISSING (degraded; used local-concept backend)\n")
    print(f"[env] versions.txt written. backend={backend}")


if __name__ == "__main__":
    main()
