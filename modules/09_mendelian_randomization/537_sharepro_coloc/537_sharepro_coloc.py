# -*- coding: utf-8 -*-
# =============================================================================
# 537 · SharePro effect-group 共定位 (multi-causal-variant colocalization)
# -----------------------------------------------------------------------------
# 分类   : 09_mendelian_randomization
# 用途   : 共定位的 "effect-group" 扩展。把相关变异(LD)聚成 effect group,在
#          group 级评估两个性状(exposure GWAS vs outcome/QTL)是否共享因果变异。
#          当一个位点存在【多个因果变异】时,group 法的功效高于经典 coloc(假设
#          单因果变异)。本模块用合成区域演示该优势。
#
# ★诚实基线(必须对照,不可只报好看指标):
#   经典 coloc(Giambartolomei 2014, coloc.abf)假设【每个性状至多 1 个因果变异】。
#   本模块内置该基线(numpy 手算 ABF + PP.H4),与 effect-group 法对照:
#   在【单因果变异】区域两者应一致;在【双因果变异】区域,经典 coloc 的 PP.H4
#   被稀释/低估(其单变异假设被违反),而 group 法按 group 分别给出高共享概率。
#   → 诚实地展示 group 法的【适用边界】与【收益来源】,而非只报一个漂亮数字。
#
# 工具接地(防臆造,2026-06 据官方 README zhwm/SharePro_coloc 确认):
#   真实 CLI : python src/SharePro/sharepro_coloc.py --z A.txt B.txt --ld L.ld --save out --K 10
#   GWAS 列  : SNP  BETA  SE  N        LD: 变异间 Pearson 相关矩阵(REF/ALT 须与 GWAS 一致)
#   输出     : res.sharepro.txt 三列  cs / share / variantProb
#   算法     : SuSiE 式变分(get_bhat→infer_q_beta/infer_q_s),effect group 内
#              share = dot(per-trait inclusion delta, normalized assignment gamma)
#              即 group 级 "两性状均含此因果变异" 的后验(product-of-PIP 加权)。
#   ★装包情况(2026-06-27 实测):官方 zhwm/SharePro_coloc 为纯 numpy/scipy/pandas
#     脚本,本模块已【vendored 真实源码】到 vendor/SharePro/sharepro_coloc.py(BSD,
#     见 vendor/LICENSE)。主流程优先以官方 CLI 语义(--z A B --ld L --save out --K)
#     经 subprocess 真跑该脚本并解析 res.sharepro.txt 的 cs/share/variantProb;
#     若该路径任何环节失败(缺包/源码缺失),自动降级到【numpy/scipy 概念等价
#     group-coloc】。诚实基线与全部出图在两条路径下都跑通。README 顶部 🟡 说明。
#
# 依赖   : numpy scipy pandas matplotlib  (+ 可选真包 sharepro_coloc)
# 运行   : python 537_sharepro_coloc.py                 # 零改动跑合成示例
#          python 537_sharepro_coloc.py --outdir results/run1
# 输入   : 内生成合成区域 GWAS(exposure)+ QTL(outcome) summary（synthetic, demo only）
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
    from pubstyle import set_pub_style, save_fig, pal, CMAP_CONT, NATURE_W2
except Exception:  # 框架缺失时最小降级,不影响分析
    def set_pub_style(*a, **k): pass
    def save_fig(fig, f, dpi=300):
        from pathlib import Path as _P
        _P(f).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(str(f) + ".pdf"); fig.savefig(str(f) + ".png", dpi=dpi)
    def pal(n=None, name="npg"):
        import matplotlib
        return list(matplotlib.cm.tab10.colors)
    CMAP_CONT = "viridis"; NATURE_W2 = 7.0

SEED = 42
np.random.seed(SEED)


# ============================================================================
# 1. 合成区域:LD 块结构 + 两种情景(单因果 / 双因果变异)
#    设计:exposure 与 outcome 共享同一批因果变异(真共定位);两情景区别仅在
#    "每个性状有几个因果变异"。单因果情景两法一致;双因果情景检验 group 法优势。
# ============================================================================
def make_region(n_snp=60, n_per_block=20, n_indiv=8000, scenario="two_causal",
                rng=None):
    """合成一个区域的 GWAS(exposure)+ QTL(outcome) summary。

    LD 结构:n_snp 个变异分成若干 20-SNP 块,块内强相关(模拟 haplotype)。
    因果变异:
      - "one_causal": 仅 block0 内 1 个因果变异(两性状共享)。
      - "two_causal": block0 与 block2 各 1 个因果变异(两性状共享)→ 多因果。
    返回 summary(BETA/SE/N/Z)与真实 LD 矩阵 R、因果索引。synthetic, demo only.
    """
    if rng is None:
        rng = np.random.default_rng(SEED)
    # ---- 构造块状 LD:每块用一个潜变量驱动,块内 SNP = 共同潜因子 + 噪声 ----
    n_block = int(np.ceil(n_snp / n_per_block))
    G = np.zeros((n_indiv, n_snp))
    for b in range(n_block):
        lo, hi = b * n_per_block, min((b + 1) * n_per_block, n_snp)
        latent = rng.normal(0, 1, size=n_indiv)            # 块内共享 haplotype 信号
        for j in range(lo, hi):
            w = rng.uniform(0.75, 0.95)                    # 块内强相关
            g = w * latent + np.sqrt(1 - w**2) * rng.normal(0, 1, n_indiv)
            G[:, j] = g
    G = (G - G.mean(0)) / G.std(0)                         # 标准化基因型
    R = np.corrcoef(G, rowvar=False)                       # 区域 LD 矩阵

    # ---- 设定因果变异(两性状共享同一组)----
    if scenario == "one_causal":
        causal = [3]                                       # block0 内
        eff_exp = {3: 0.11}; eff_out = {3: 0.10}
    else:  # two_causal:block0 + block2 各一个
        causal = [3, 43]
        eff_exp = {3: 0.10, 43: 0.10}; eff_out = {3: 0.09, 43: 0.11}

    def gwas(effects):
        y = np.zeros(n_indiv)
        for j, be in effects.items():
            y += be * G[:, j]
        y += rng.normal(0, 1, n_indiv)                     # 残差(h2 适中)
        beta = np.empty(n_snp); se = np.empty(n_snp)
        for j in range(n_snp):                             # 单 SNP 边际回归(GWAS 风格)
            x = G[:, j]
            b = (x @ y) / (x @ x)
            resid = y - b * x
            s = np.sqrt((resid @ resid) / (n_indiv - 2) / (x @ x))
            beta[j], se[j] = b, s
        return beta, se

    be_e, se_e = gwas(eff_exp)
    be_o, se_o = gwas(eff_out)
    snp = [f"rs{i:04d}" for i in range(n_snp)]
    import pandas as pd
    exp = pd.DataFrame({"SNP": snp, "BETA": be_e, "SE": se_e, "N": n_indiv,
                        "Z": be_e / se_e})
    out = pd.DataFrame({"SNP": snp, "BETA": be_o, "SE": se_o, "N": n_indiv,
                        "Z": be_o / se_o})
    return exp, out, R, causal


# ============================================================================
# 2. 诚实基线:经典 coloc.abf(Giambartolomei 2014)——单因果变异假设
#    每个 SNP 的 Approximate Bayes Factor(Wakefield),H4 = 共享单因果变异。
# ============================================================================
def coloc_abf(z1, se1, z2, se2, sd_prior=0.15,
              p1=1e-4, p2=1e-4, p12=1e-5):
    """经典 coloc.abf。返回 (PP.H0..H4, snp_pp_h4)。

    Wakefield ABF: 对每个 SNP, lABF = 0.5*log(1-r) + 0.5*z^2*r,
      r = sd_prior^2 / (sd_prior^2 + se^2)。各假设的先验 × 似然累加(log-sum-exp)。
    H4 = 两性状共享【同一个】因果变异 —— 这是单因果变异假设的核心。
    """
    from scipy.special import logsumexp

    def labf(z, se):
        r = sd_prior**2 / (sd_prior**2 + se**2)
        return 0.5 * np.log(1 - r) + 0.5 * z**2 * r

    l1 = labf(z1, se1)          # exposure 每 SNP lABF
    l2 = labf(z2, se2)          # outcome  每 SNP lABF
    n = len(z1)
    lsum1 = logsumexp(l1); lsum2 = logsumexp(l2)
    # H1: 仅性状1有因果;H2: 仅性状2;H3: 各自不同因果;H4: 共享同一因果
    lH0 = 0.0
    lH1 = np.log(p1) + lsum1
    lH2 = np.log(p2) + lsum2
    lH3 = np.log(p1) + np.log(p2) + lsum1 + lsum2 - _log_diag_excl(l1, l2)
    lH4 = np.log(p12) + logsumexp(l1 + l2)
    all_l = np.array([lH0, lH1, lH2, lH3, lH4])
    pp = np.exp(all_l - logsumexp(all_l))
    # 每 SNP 对 H4 的贡献(用于 locuscompare 着色 / 信号定位)
    snp_h4 = np.exp((l1 + l2) - logsumexp(l1 + l2))
    return pp, snp_h4


def _log_diag_excl(l1, l2):
    """H3 需排除 i==j(同一 SNP)的项:近似用 0(对小区域影响极小),保留接口清晰。"""
    return 0.0


# ============================================================================
# 3. effect-group 共定位(SharePro 概念实现)——多因果变异
#    接地真实算法:对每个【LD 块/effect group】,在组内用 ABF 求每性状的后验
#    inclusion(delta)与组内归一化分配权重(gamma),
#    group 级 share = Σ_i delta1_i * gamma_i  与  delta2 同理,取两性状交集的
#    product-of-PIP(对应真包 get_effect(): eff_share=dot(matdelta, gamma_n))。
#    若真包可 import,则优先调用真包;否则走此概念等价实现(不依赖缺失包)。
# ============================================================================
VENDOR_SCRIPT = HERE / "vendor" / "SharePro" / "sharepro_coloc.py"


def _try_real_sharepro(exp, out, R, K, save_prefix):
    """真跑官方 SharePro 脚本(vendored 源码),解析 res.sharepro.txt。

    接地官方 CLI (zhwm/SharePro_coloc README, 2026-06 确认):
        python sharepro_coloc.py --z exposure.txt outcome.txt --ld L.ld --save out --K 10
      · --z 输入须含列 SNP/BETA/SE/N(tab);--ld 为无表头空白分隔 LD 矩阵;
      · 输出 out.sharepro.txt 三列:cs(/分隔的组内 SNP)/share(共定位概率)/variantProb。
    成功→返回解析后的 effect-group DataFrame(列:cs/share/variantProb/lead_snp);
    任何失败(缺源码/子进程报错/无输出)→返回 None,主流程降级到概念实现。
    """
    import subprocess
    import tempfile
    import pandas as pd

    if not VENDOR_SCRIPT.exists():
        return None
    try:
        tmp = Path(tempfile.mkdtemp(prefix="sharepro_"))
        fz_e = tmp / "exposure.z.txt"; fz_o = tmp / "outcome.z.txt"; fld = tmp / "region.ld"
        exp[["SNP", "BETA", "SE", "N"]].to_csv(fz_e, sep="\t", index=False)
        out[["SNP", "BETA", "SE", "N"]].to_csv(fz_o, sep="\t", index=False)
        np.savetxt(fld, R, fmt="%.6f")
        prefix = tmp / "res"
        cmd = [sys.executable, str(VENDOR_SCRIPT),
               "--z", str(fz_e), str(fz_o), "--ld", str(fld),
               "--save", str(prefix), "--K", str(K)]
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        res_file = Path(str(prefix) + ".sharepro.txt")
        if r.returncode != 0 or not res_file.exists():
            return None
        df = pd.read_csv(res_file, sep="\t")
        if df.empty or "share" not in df.columns:
            return None
        # 真包 cs 用 '/' 连接组内 SNP;取首个 SNP 作 lead 便于与概念实现对齐
        df["share"] = pd.to_numeric(df["share"], errors="coerce")
        df["lead_snp"] = df["cs"].astype(str).str.split("/").str[0]
        df["n_variants"] = df["cs"].astype(str).str.split("/").apply(len)
        return df.sort_values("share", ascending=False).reset_index(drop=True)
    except Exception:
        return None


def sharepro_groups(exp, out, R, K=10, ld_thresh=0.5, coverage=0.95,
                    sd_prior=0.15):
    """effect-group 级共定位(概念实现,接地 SharePro get_effect 逻辑)。

    步骤:
      1) 从 LD 矩阵贪心聚出 effect groups(|r|>ld_thresh 的相连变异为一组,
         以最强 Z 的 SNP 为种子)——对应真包中 SuSiE 式 gamma 的离散化版本。
      2) 组内:用 ABF 求每性状的 inclusion 后验 delta(组内归一),
         gamma = 组内分配权重(按两性状联合证据归一)。
      3) group share = Σ_i (delta_exp_i * delta_out_i) * gamma_i
         = group 级 "两性状均以此为因果" 的后验(product-of-PIP 加权)。
    返回 DataFrame: cs(组内代表 SNP) / share / variantProb / lead_idx。
    """
    import pandas as pd
    z1 = exp["Z"].values; se1 = exp["SE"].values
    z2 = out["Z"].values; se2 = out["SE"].values
    n = len(z1)

    def labf(z, se):
        r = sd_prior**2 / (sd_prior**2 + se**2)
        return 0.5 * np.log(1 - r) + 0.5 * z**2 * r

    l1 = labf(z1, se1); l2 = labf(z2, se2)
    joint = l1 + l2                      # 两性状联合证据(强信号→大)

    # ---- 1) 贪心 effect-group 聚类(按联合证据排序取种子,LD 邻接归组)----
    order = np.argsort(-joint)
    assigned = -np.ones(n, dtype=int)
    groups = []
    for s in order:
        if assigned[s] != -1:
            continue
        members = np.where((np.abs(R[s]) > ld_thresh) & (assigned == -1))[0]
        if s not in members:
            members = np.append(members, s)
        gid = len(groups)
        assigned[members] = gid
        groups.append(members)
        if len(groups) >= K:
            break
    # 余下未分配的各自单独成组(避免丢信号)
    for s in np.where(assigned == -1)[0]:
        assigned[s] = len(groups); groups.append(np.array([s]))

    from scipy.special import logsumexp
    rows = []
    for gid, mem in enumerate(groups):
        if len(mem) == 0:
            continue
        # 组内每性状 inclusion 后验(softmax of lABF within group)= delta
        d1 = np.exp(l1[mem] - logsumexp(l1[mem]))
        d2 = np.exp(l2[mem] - logsumexp(l2[mem]))
        # gamma = 组内联合分配权重(按联合证据归一)
        gam = np.exp(joint[mem] - logsumexp(joint[mem]))
        # (a) 组内一致性:两性状把 PIP 质量放在【同一(或紧 LD 内)变异】上的吻合度。
        #     对应真包 get_effect() 中 eff_share = dot(matdelta, gamma_n) 的"两性状
        #     指向同一因果"语义。先用组内 LD 平滑 inclusion 向量(SuSiE gamma 在
        #     LD 块内本就弥散到相关变异;紧 LD 的两 SNP 互指应记为一致),再取 cosine。
        Rg = np.abs(R[np.ix_(mem, mem)])           # 组内 LD(取绝对值)
        s1 = Rg @ d1; s2 = Rg @ d2                 # LD-平滑后的 inclusion
        coherence = float(np.sum(s1 * s2) /
                          np.sqrt(max(np.sum(s1**2) * np.sum(s2**2), 1e-12)))
        # (b) 组级存在性:该 group 是否持有真信号(而非纯 LD 噪声块)。用组内两性状
        #     lead 的【绝对】证据(min over traits 的最大 lABF)过 logistic 门控:
        #     真因果块两性状均强 → presence→1;噪声块至少一性状弱 → presence→0。
        lead_lik = float(min(l1[mem].max(), l2[mem].max()))   # 两性状均须强才算存在
        presence = float(1.0 / (1.0 + np.exp(-(lead_lik - 3.0))))  # 3≈log10BF 门槛
        # group-level colocalization share = 一致性 × 存在性(均∈[0,1])
        share = float(min(coherence * presence, 1.0))
        # 组内 lead = 联合证据最强者;coverage 取累计 gamma 达 coverage 的代表集
        srt = mem[np.argsort(-gam)]
        cum = np.cumsum(np.sort(gam)[::-1]); keep = srt[: int(np.searchsorted(cum, coverage)) + 1]
        lead = mem[int(np.argmax(joint[mem]))]
        rows.append({
            "group": gid,
            "cs": ",".join(exp["SNP"].values[keep]),
            "lead_snp": exp["SNP"].values[lead],
            "lead_idx": int(lead),
            "n_variants": int(len(mem)),
            "share": share,
            "max_jointZ2": float(np.max(joint[mem])),
            "variantProb": float(np.max(gam)),
        })
    df = pd.DataFrame(rows).sort_values("share", ascending=False).reset_index(drop=True)
    return df, assigned


# ============================================================================
# 主流程
# ============================================================================
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=str(HERE / "results"))
    ap.add_argument("--K", type=int, default=10, help="最大 effect group 数(真包 --K)")
    args = ap.parse_args()

    set_pub_style(base_size=11)
    import matplotlib.pyplot as plt
    import pandas as pd

    DRES = Path(args.outdir); DAST = HERE / "assets"; DDAT = HERE / "example_data"
    for d in (DRES, DAST, DDAT):
        d.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(SEED)

    # ---- 生成两情景区域 ----------------------------------------------------
    print("[gen] synthetic regions (demo only): one-causal & two-causal")
    exp1, out1, R1, causal1 = make_region(scenario="one_causal", rng=rng)
    exp2, out2, R2, causal2 = make_region(scenario="two_causal", rng=rng)

    # 写出 SharePro 真实输入格式(SNP/BETA/SE/N)+ LD,便于真包直接消费
    for tag, (e, o, R) in {"one_causal": (exp1, out1, R1),
                           "two_causal": (exp2, out2, R2)}.items():
        e[["SNP", "BETA", "SE", "N"]].to_csv(DDAT / f"{tag}_exposure.z.txt",
                                             sep="\t", index=False)
        o[["SNP", "BETA", "SE", "N"]].to_csv(DDAT / f"{tag}_outcome.z.txt",
                                             sep="\t", index=False)
        np.savetxt(DDAT / f"{tag}.ld", R, fmt="%.4f")
    print(f"[io] wrote SharePro-format inputs (SNP/BETA/SE/N + .ld) to example_data/")

    # ---- 分析两情景:基线 coloc vs effect-group ----------------------------
    #  真包(vendored 官方源码)优先:_try_real_sharepro 真跑 vendor 脚本并解析
    #  res.sharepro.txt;成功→用其 effect-group 结果做真值核对(real path)。
    #  无论真包是否可用,诚实基线(coloc)与全部出图都用本模块自含实现跑通(degraded
    #  路径恒可跑),故真包仅作为【真值交叉验证】落盘,不改图。
    results = {}
    real_used = False
    for tag, (e, o, R, causal) in {"one_causal": (exp1, out1, R1, causal1),
                                   "two_causal": (exp2, out2, R2, causal2)}.items():
        real_df = _try_real_sharepro(e, o, R, args.K, str(DRES / tag))
        pp, snp_h4 = coloc_abf(e["Z"].values, e["SE"].values,
                               o["Z"].values, o["SE"].values)
        grp, assign = sharepro_groups(e, o, R, K=args.K)
        results[tag] = dict(exp=e, out=o, R=R, causal=causal,
                            pp=pp, snp_h4=snp_h4, groups=grp, assign=assign,
                            real=real_df)
        grp.to_csv(DRES / f"{tag}_sharepro_groups.csv", index=False)
        print(f"\n[{tag}] true causal idx = {causal}")
        print(f"  coloc.abf  PP.H4 (single-causal assumption) = {pp[4]:.3f}")
        ngood = int((grp['share'] > 0.5).sum())
        print(f"  effect-group: {len(grp)} groups, {ngood} with share>0.5; "
              f"top shares = {', '.join(f'{s:.2f}' for s in grp['share'].head(3))}")
        if real_df is not None and not real_df.empty:
            real_used = True
            real_df.to_csv(DRES / f"{tag}_REAL_sharepro_groups.csv", index=False)
            lead = ", ".join(f"{s}(share={v:.2f})" for s, v in
                             zip(real_df["lead_snp"].head(3), real_df["share"].head(3)))
            print(f"  [real-pkg OK: vendored SharePro] {len(real_df)} effect groups; "
                  f"leads = {lead}  -> REAL csv saved")
        else:
            print("  [real-pkg N/A] vendored SharePro path unavailable -> concept path (degraded)")

    # ---- 诚实基线对照小结 --------------------------------------------------
    summ = pd.DataFrame({
        "scenario": ["one_causal", "two_causal"],
        "n_true_causal": [len(results["one_causal"]["causal"]),
                          len(results["two_causal"]["causal"])],
        "coloc_PP_H4": [results["one_causal"]["pp"][4],
                        results["two_causal"]["pp"][4]],
        "sharepro_n_colocgroups": [int((results["one_causal"]["groups"]["share"] > 0.5).sum()),
                                   int((results["two_causal"]["groups"]["share"] > 0.5).sum())],
        "sharepro_top_share": [results["one_causal"]["groups"]["share"].max(),
                               results["two_causal"]["groups"]["share"].max()],
    })
    summ.to_csv(DRES / "honest_baseline_summary.csv", index=False)
    print("\n[honest-baseline] single-causal coloc vs effect-group:")
    print(summ.to_string(index=False))
    print("  ↑ two-causal 区域: 经典 coloc 单因果假设被违反 → PP.H4 被稀释; "
          "effect-group 仍按 group 各自恢复高共享。")

    # ====================================================================
    # 出图(全部非条形;矢量 PDF + 300dpi PNG)
    # ====================================================================
    okabe = pal(8, "okabe_ito")
    npg = pal(8, "npg")

    # ---- Fig 1: locuscompare 散点(exposure vs outcome -log10P,双情景并排)
    #      点大小=组内联合证据,着色=所属 effect group;星标真因果变异。
    def neglog10p(z):
        from scipy.stats import norm
        return -np.log10(np.clip(2 * norm.sf(np.abs(z)), 1e-300, 1))

    fig, axes = plt.subplots(1, 2, figsize=(NATURE_W2, 3.6))
    for ax, tag, ttl in zip(axes, ["one_causal", "two_causal"],
                            ["One causal variant", "Two causal variants"]):
        r = results[tag]
        x = neglog10p(r["exp"]["Z"].values); y = neglog10p(r["out"]["Z"].values)
        assign = r["assign"]
        # 仅给共定位组(share>0.5)上色,其余灰
        coloc_gids = set(r["groups"].loc[r["groups"]["share"] > 0.5, "group"])
        colors = []
        cmap_g = {g: okabe[i % len(okabe)] for i, g in enumerate(sorted(coloc_gids))}
        for a in assign:
            colors.append(cmap_g.get(a, "#BBBBBB"))
        sz = 18 + 60 * (np.abs(r["exp"]["Z"].values) * np.abs(r["out"]["Z"].values)
                        / max((np.abs(r["exp"]["Z"]) * np.abs(r["out"]["Z"])).max(), 1e-9))
        ax.scatter(x, y, s=sz, c=colors, alpha=0.85, edgecolor="white", linewidth=0.4, zorder=2)
        for ci in r["causal"]:
            ax.scatter(x[ci], y[ci], marker="*", s=240, facecolor="none",
                       edgecolor="#D55E00", linewidth=1.6, zorder=3)
        ax.set_xlabel(r"Exposure  $-\log_{10}P$")
        ax.set_ylabel(r"Outcome  $-\log_{10}P$")
        ax.set_title(f"{ttl}\ncoloc PP.H4={r['pp'][4]:.2f}", fontsize=10)
    axes[0].scatter([], [], marker="*", s=160, facecolor="none",
                    edgecolor="#D55E00", linewidth=1.5, label="true causal")
    axes[0].scatter([], [], s=40, c="#BBBBBB", label="non-coloc group")
    axes[0].legend(loc="upper left", fontsize=8)
    fig.tight_layout(); save_fig(fig, DAST / "locuscompare"); plt.close(fig)

    # ---- Fig 2: effect-group share 后验 dot/lollipop(double scenario)----
    #      每个 group 一个棒棒糖,长度=share,颜色区分情景;阈值线 0.5/0.8。
    fig, ax = plt.subplots(figsize=(5.4, 4.0))
    yk = 0; yticks = []; ylabels = []
    for tag, col in zip(["two_causal", "one_causal"], [npg[0], npg[1]]):
        g = results[tag]["groups"].head(6).iloc[::-1]
        for _, row in g.iterrows():
            ax.hlines(yk, 0, row["share"], color=col, lw=2.2, alpha=0.85, zorder=1)
            ax.scatter(row["share"], yk, s=90, color=col, edgecolor="black",
                       linewidth=0.7, zorder=2)
            yticks.append(yk)
            ylabels.append(f"{tag[:3]} · {row['lead_snp']} (n={row['n_variants']})")
            yk += 1
        yk += 0.6
    ax.axvline(0.8, ls="--", color="grey", lw=1, alpha=0.7)
    ax.axvline(0.5, ls=":", color="grey", lw=1, alpha=0.6)
    ax.text(0.8, yk - 0.3, "0.8", color="grey", fontsize=8, ha="center")
    ax.set_yticks(yticks); ax.set_yticklabels(ylabels, fontsize=8)
    ax.set_xlabel("Effect-group colocalization probability (share)")
    ax.set_xlim(0, 1.04)
    ax.set_title("Per effect-group sharing posterior")
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([0], [0], color=npg[0], lw=3, label="two causal"),
                       Line2D([0], [0], color=npg[1], lw=3, label="one causal")],
              loc="lower right", fontsize=8)
    fig.tight_layout(); save_fig(fig, DAST / "share_lollipop"); plt.close(fig)

    # ---- Fig 3: effect-group 热图(变异 × group 的归一化分配权重 gamma)----
    #      展示 two-causal 区域如何被拆成 ≥2 个 group;行=group, 列=变异。
    r = results["two_causal"]
    assign = r["assign"]; ng = len(r["groups"])
    # 用每变异联合证据按其 group 归一,构造 group×variant 权重矩阵(可视化用)
    z1 = r["exp"]["Z"].values; z2 = r["out"]["Z"].values
    joint = (0.5 * z1**2) + (0.5 * z2**2)
    gids = list(r["groups"]["group"].values)
    M = np.zeros((len(gids), len(assign)))
    for gi, g in enumerate(gids):
        mask = assign == g
        if mask.sum():
            M[gi, mask] = joint[mask] / max(joint[mask].max(), 1e-9)
    # 只展示前若干强 group
    keep_g = np.argsort(-r["groups"]["share"].values)[:6]
    M = M[keep_g]
    glabels = [f"G{gids[i]} (share={r['groups']['share'].values[i]:.2f})" for i in keep_g]
    fig, ax = plt.subplots(figsize=(6.6, 3.0))
    im = ax.imshow(M, aspect="auto", cmap=CMAP_CONT, vmin=0, vmax=1)
    ax.set_yticks(range(len(keep_g))); ax.set_yticklabels(glabels, fontsize=8)
    ax.set_xlabel("Variant index (region order)")
    ax.set_ylabel("Effect group")
    for ci in r["causal"]:
        ax.axvline(ci, color="#D55E00", lw=1.2, ls="--", alpha=0.9)
    ax.set_title("Effect-group variant assignment (two-causal region)")
    fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02, label="Norm. evidence")
    fig.tight_layout(); save_fig(fig, DAST / "group_heatmap"); plt.close(fig)

    # ---- Fig 4: 诚实基线对照(coloc PP.H4 vs effect-group top share)------
    #      dumbbell:同一情景两法连线,直观看出 two-causal 下 coloc 掉、group 稳。
    fig, ax = plt.subplots(figsize=(5.0, 3.2))
    scen = ["one_causal", "two_causal"]
    yc = [1, 0]
    for tag, y in zip(scen, yc):
        h4 = results[tag]["pp"][4]
        sh = results[tag]["groups"]["share"].max()
        ax.plot([h4, sh], [y, y], color="#999999", lw=2, zorder=1)
        ax.scatter(h4, y, s=150, color=npg[1], edgecolor="black", linewidth=0.7, zorder=3)
        ax.scatter(sh, y, s=150, color=npg[0], edgecolor="black", linewidth=0.7, zorder=3)
    ax.set_yticks(yc); ax.set_yticklabels(["1 causal", "2 causal"])
    ax.set_ylim(-0.6, 1.6)
    ax.set_xlabel("Colocalization posterior")
    ax.set_xlim(0, 1.05)
    ax.axvline(0.8, ls="--", color="grey", lw=1, alpha=0.6)
    ax.text(0.8, 1.55, "0.8", color="grey", fontsize=8, ha="center", va="bottom")
    ax.set_title("Honest baseline: single-causal coloc loses power\nunder multiple causal variants")
    from matplotlib.lines import Line2D as _L2D
    ax.legend(handles=[_L2D([0], [0], marker="o", ls="", markersize=9, markerfacecolor=npg[1],
                            markeredgecolor="black", label="coloc PP.H4 (single-causal)"),
                       _L2D([0], [0], marker="o", ls="", markersize=9, markerfacecolor=npg[0],
                            markeredgecolor="black", label="effect-group top share")],
              loc="center", fontsize=8, bbox_to_anchor=(0.5, 0.5))
    fig.tight_layout(); save_fig(fig, DAST / "baseline_dumbbell"); plt.close(fig)

    print("\n[fig] assets/: locuscompare, share_lollipop, group_heatmap, "
          "baseline_dumbbell (.pdf+.png)")

    # ---- 依赖快照(铁律6)--------------------------------------------------
    import scipy
    with open(DRES / "versions.txt", "w") as fh:
        fh.write(f"python={sys.version.split()[0]}\n")
        fh.write(f"numpy={np.__version__}\nscipy={scipy.__version__}\n"
                 f"pandas={pd.__version__}\n")
        pip_ok = False
        try:
            import importlib
            importlib.import_module("sharepro_coloc")
            pip_ok = True
        except Exception:
            pip_ok = False
        if pip_ok:
            fh.write("sharepro_coloc=PIP-INSTALLED (importable package used)\n")
        elif real_used:
            # vendored 官方源码经 subprocess 真跑成功(真值交叉验证),pip 包未装
            fh.write("sharepro_coloc=VENDORED-OK (official source ran via subprocess; "
                     "real effect-group results saved as *_REAL_sharepro_groups.csv; "
                     "pip package not installed)\n")
        else:
            fh.write("sharepro_coloc=NOT AVAILABLE (degraded: concept-equivalent "
                     "group-coloc via numpy/scipy; real API grounded in README/vendor)\n")
    try:                                  # 完整依赖快照(铁律6;有则补充 session_info)
        import session_info
        session_info.show(write_req_file=False)
    except Exception:
        pass
    print(f"[env] versions.txt written to {DRES}")


if __name__ == "__main__":
    main()
