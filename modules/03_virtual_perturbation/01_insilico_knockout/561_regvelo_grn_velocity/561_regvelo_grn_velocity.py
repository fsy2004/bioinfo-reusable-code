"""561 · RegVelo — 基因调控网络驱动的 RNA velocity 与 regulon 虚拟敲除.

上游: RegVelo (Wang W, Hu Z, Weiler P, Mayes S, Lange M, Fountain DM, Haug JO, Wang J,
Xue Z, Sauka-Spengler T, Theis FJ. Cell 2026;189:3773-3800.e44.
doi:10.1016/j.cell.2026.04.022 · PMID 42119563) · https://github.com/theislab/regvelo

本模块两条路径:
  ① RegVelo 官方路径 (--run-regvelo): 需 `pip install regvelo` + GPU。守卫式封装,
     函数名/参数全部核对自本地克隆源码 (见 README「源码接地」表)。
  ② GRN-linear 代理基线 (默认, CPU 可跑): 复刻 RegVelo 的**函数形式与下游打分逻辑**
     —— alpha = softplus(W @ s + b); du/dt = alpha - beta*u; ds/dt = beta*u - gamma*s
     (接地 src/regvelo/_module.py velocity_encoder.forward 第 293-313 行),
     W 用先验骨架约束的岭回归拟合而非 VAE 训练; TF 敲除 = 把 W 中该 TF 的列清零
     (接地 tools/_in_silico_block_simulation.py), 再走 CellRank 重算 fate probability,
     用上游同款 abundance_test (metrics/_abundance_test.py: ranksums + ROC-AUC + BH)。
     另出**无 GRN 的纯 scVelo 对照**,GRN 路径必须与之并排报告,不单独下结论。

零改动即跑: python 561_regvelo_grn_velocity.py
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import warnings

warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
EXAMPLE = os.path.join(HERE, "example_data")
RESULTS = os.path.join(HERE, "results")
ASSETS = os.path.join(HERE, "assets")
FRAMEWORK = os.path.abspath(os.path.join(HERE, "..", "..", "..", "_framework"))
if FRAMEWORK not in sys.path:
    sys.path.insert(0, FRAMEWORK)

SEED = 0

# ------------------------------------------------------------------ 参数区
P_N_CELLS = 500          # 合成细胞数
P_N_TF = 6               # 合成 TF 数(= GRN 的调控子)
P_N_TARGET = 54          # 合成靶基因数
P_RIDGE = 1.0            # GRN 岭回归正则强度
P_N_STATES = 3           # GPCCA macrostate 数(1 祖细胞 + 2 终端)


# ============================================================ ① 合成示例数据
def make_example_data(outdir: str, n_cells: int = P_N_CELLS, n_tf: int = P_N_TF,
                      n_target: int = P_N_TARGET, seed: int = SEED):
    """生成分叉分化的 spliced/unspliced 合成数据 + 先验 GRN 骨架,落盘到 example_data/。

    synthetic, for demo only. 设计:
      - 潜时 t∈[0,1];t>0.4 后分叉为 Branch_A / Branch_B 两个终端命运;
      - 6 个 TF:TF0 早期、TF1 中期、TF2/TF3 = Branch_A 驱动、TF4/TF5 = Branch_B 驱动;
      - 靶基因程序 = softplus(Σ w·TF程序 + b),即真实存在一条 TF→target 的因果链;
      - unspliced 超前 spliced(u = p + tau·dp/dt),这正是 velocity 读取的信号。
    这样 TF2/TF3 敲除应当降低 Branch_A 的 fate probability —— 内置阳性对照。
    """
    import numpy as np
    import pandas as pd
    import anndata as ad

    rng = np.random.default_rng(seed)
    t = np.sort(rng.uniform(0.0, 1.0, n_cells))
    branch = rng.integers(0, 2, n_cells)               # 0=A, 1=B
    tf_names = [f"TF{i}" for i in range(n_tf)]
    tg_names = [f"TG{i:02d}" for i in range(n_target)]

    def tf_program(tt, bb):
        """TF 活性程序:shape (len(tt), n_tf)。分叉前两支相同。"""
        tt = np.atleast_1d(tt).astype(float)
        bb = np.atleast_1d(bb)
        split = np.clip((tt - 0.4) / 0.6, 0.0, 1.0)     # 分叉权重
        out = np.empty((tt.size, n_tf))
        out[:, 0] = 2.5 * np.exp(-((tt - 0.05) ** 2) / 0.02)          # 早期
        out[:, 1] = 2.0 * np.exp(-((tt - 0.50) ** 2) / 0.03)          # 中期
        a_on, b_on = (bb == 0).astype(float), (bb == 1).astype(float)
        out[:, 2] = 3.0 * split * a_on
        out[:, 3] = 2.2 * split * a_on + 0.3 * tt
        out[:, 4] = 3.0 * split * b_on
        out[:, 5] = 2.2 * split * b_on + 0.3 * tt
        return out

    # 先验 GRN:每个靶基因随机挂 1-2 个 TF(有符号权重)
    W_true = np.zeros((n_target, n_tf))
    for g in range(n_target):
        regs = rng.choice(n_tf, size=int(rng.integers(1, 3)), replace=False)
        W_true[g, regs] = rng.uniform(0.6, 1.6, regs.size) * rng.choice([1.0, -1.0], regs.size)
    b_true = rng.uniform(-0.5, 0.5, n_target)

    def target_program(tt, bb):
        return np.log1p(np.exp(np.clip(tf_program(tt, bb) @ W_true.T + b_true, -30, 30)))

    def full_program(tt, bb):
        return np.hstack([tf_program(tt, bb), target_program(tt, bb)])

    prog = full_program(t, branch)
    h = 1e-3
    dprog = (full_program(np.clip(t + h, 0, 1), branch)
             - full_program(np.clip(t - h, 0, 1), branch)) / (2 * h)

    tau = 0.06                                          # unspliced 超前量
    scale = 12.0
    s_rate = np.clip(prog, 0.05, None) * scale
    u_rate = np.clip(prog + tau * dprog, 0.05, None) * scale * 0.5
    spliced = rng.poisson(s_rate).astype("float32")
    unspliced = rng.poisson(u_rate).astype("float32")

    A = ad.AnnData(spliced.copy())
    A.layers["spliced"] = spliced
    A.layers["unspliced"] = unspliced
    A.var_names = tf_names + tg_names
    A.obs_names = [f"C{i:04d}" for i in range(n_cells)]
    A.obs["true_time"] = t
    A.obs["cell_state"] = pd.Categorical(
        np.where(t < 0.4, "Progenitor", np.where(branch == 0, "Branch_A", "Branch_B")),
        categories=["Progenitor", "Branch_A", "Branch_B"])
    A.uns["synthetic"] = "synthetic, for demo only"

    os.makedirs(outdir, exist_ok=True)
    h5 = os.path.join(outdir, "synthetic_velocity.h5ad")
    A.write_h5ad(h5)

    # 先验 GRN 骨架:行=target, 列=regulator —— 与上游 pp.set_prior_grn(gt_net=) 约定一致
    grn = pd.DataFrame((W_true != 0).astype(int), index=tg_names, columns=tf_names)
    grn.index.name = "target"
    csv = os.path.join(outdir, "synthetic_prior_grn.csv")
    grn.to_csv(csv)

    with open(os.path.join(outdir, "README.txt"), "w", encoding="utf-8") as fh:
        fh.write("synthetic, for demo only — generated by 561_regvelo_grn_velocity.py "
                 f"(seed={seed}).\n"
                 "synthetic_velocity.h5ad : AnnData, layers['spliced'|'unspliced'], "
                 "obs['true_time','cell_state'].\n"
                 "synthetic_prior_grn.csv : prior GRN skeleton, rows=targets, cols=regulators "
                 "(0/1), matching regvelo.pp.set_prior_grn(gt_net=...).\n")
    print(f"[561] wrote example data -> {h5}")
    print(f"[561] wrote example data -> {csv}")
    return h5, csv


# ============================================================ 评估辅助
def _session_info():
    """记录依赖版本快照,写进 summary,便于复现(铁律6)。"""
    import platform
    import importlib
    out = {"python": platform.python_version(), "platform": platform.platform()}
    for m in ("numpy", "pandas", "scipy", "sklearn", "matplotlib",
              "anndata", "scanpy", "scvelo", "cellrank", "statsmodels"):
        try:
            out[m] = getattr(importlib.import_module(m), "__version__", "?")
        except Exception:
            out[m] = "not installed"
    return out


def lineage_rho(pt, true_t, states):
    """按**谱系**算 pseudotime 与真实潜时的 Spearman |rho|。

    ★ 分叉数据上把两支混在一起算一个 pooled rho 是错的:velocity_pseudotime 以单一
      root 为起点,两支各自单调,但合并后互相抵消 → 数值被人为压低。正确做法是
      每条谱系(祖细胞 + 该终端分支)各算一次再取均值。
    """
    import numpy as np
    from scipy.stats import spearmanr

    pt, true_t, states = np.asarray(pt, float), np.asarray(true_t, float), np.asarray(states)
    branches = [s for s in np.unique(states) if s != "Progenitor"]
    per = {}
    for b in branches:
        m = (states == b) | (states == "Progenitor")
        if m.sum() > 10:
            per[b] = round(float(abs(spearmanr(pt[m], true_t[m]).correlation)), 3)
    per["mean"] = round(float(np.mean(list(per.values()))), 3) if per else float("nan")
    return per


# ============================================================ ② 无 GRN 对照
def scvelo_baseline(adata):
    """诚实对照:标准 scVelo(不含任何 GRN 信息)。GRN 路径必须与它并排报告。"""
    import scvelo as scv
    import scanpy as sc
    import numpy as np
    from scipy.stats import spearmanr

    scv.pp.filter_and_normalize(adata, min_shared_counts=0)
    scv.pp.moments(adata, n_pcs=20, n_neighbors=30)   # 30 而非默认 15:合成数据基因少,
    #   15 邻居时晚期 Branch_A 会在 KNN 图上断成孤立分量,UMAP 画出假的"离岛"
    mode_used = None
    for mode in ("stochastic", "deterministic"):
        try:
            scv.tl.velocity(adata, mode=mode)
            mode_used = mode
            break
        except Exception as e:
            print(f"       scVelo mode='{mode}' failed ({type(e).__name__}); next")
    if mode_used is None:
        raise RuntimeError("scVelo could not fit a velocity field on this input")
    scv.tl.velocity_graph(adata, n_jobs=1)
    adata.layers["velocity_scvelo"] = np.asarray(adata.layers["velocity"]).copy()

    # 展示用二维基:**扩散图**而非 UMAP。分叉轨迹上 UMAP 常把晚期分支切成孤立"离岛"
    # (KNN 图断成多个连通分量),流线图会因此失真;diffusion map 保连续性,是速度场
    # 叠加的更合适基底(scvelo/CellRank 教程亦常用)。
    sc.tl.diffmap(adata, n_comps=6)
    adata.obsm["X_traj"] = np.asarray(adata.obsm["X_diffmap"])[:, 1:3]
    # ★ 用完即清:scvelo 的 velocity_pseudotime 继承 scanpy DPT,若发现已存在的
    #   `X_diffmap` / `diffmap_evals` 会拿它当自己的基并因维数不符抛 ValueError。
    adata.obsm.pop("X_diffmap", None)
    adata.uns.pop("diffmap_evals", None)

    out = {"mode": mode_used}
    try:
        scv.tl.velocity_pseudotime(adata)
        adata.obs["pt_scvelo"] = adata.obs["velocity_pseudotime"].values
        if "true_time" in adata.obs and "cell_state" in adata.obs:
            out["scvelo_lineage_rho"] = lineage_rho(
                adata.obs["pt_scvelo"], adata.obs["true_time"],
                adata.obs["cell_state"].astype(str))
    except Exception as e:
        out["pseudotime"] = f"skipped: {type(e).__name__}"
    return out


# ============================================================ ③ GRN-linear 代理
def fit_grn_velocity(adata, skeleton, ridge: float = P_RIDGE):
    """在先验骨架约束下拟合 alpha = softplus(W·s + b),复刻 RegVelo 的转录率函数形式。

    接地 src/regvelo/_module.py::velocity_encoder.forward
        alpha = clamp(softplus(fc1(s)), 0, 50); du = alpha - beta*u; ds = beta*u - gamma*s
    这里 beta 取 scVelo 惯例 1,gamma 取 scVelo 稳态回归得到的 var['velocity_gamma'];
    W 用**岭回归**拟合(不是 VAE 变分训练)—— 这是代理,不是 RegVelo 本身。
    """
    import numpy as np

    Ms = np.asarray(adata.layers["Ms"])
    Mu = np.asarray(adata.layers["Mu"])
    genes = list(adata.var_names)
    tfs = [g for g in skeleton.columns if g in genes]
    if not tfs:
        raise ValueError("prior GRN columns (regulators) do not intersect adata.var_names")
    tf_idx = [genes.index(g) for g in tfs]
    S_tf = Ms[:, tf_idx]

    gamma = np.asarray(adata.var["velocity_gamma"]).astype(float)
    gamma = np.nan_to_num(gamma, nan=1.0, posinf=1.0, neginf=1.0)

    # alpha 的观测代理:准稳态下 alpha ≈ beta*u(beta=1)
    alpha_obs = np.clip(Mu, 1e-3, 50.0)
    y = np.log(np.expm1(alpha_obs))                      # softplus 反函数

    W = np.zeros((adata.n_vars, len(tfs)))
    b = np.zeros(adata.n_vars)
    Xc = S_tf - S_tf.mean(0, keepdims=True)
    for gi, gene in enumerate(genes):
        if gene in skeleton.index:
            regs = [k for k, tf in enumerate(tfs) if skeleton.loc[gene, tf] != 0]
        else:                                   # TF 自身:允许被其他 TF 调控
            regs = [k for k, tf in enumerate(tfs) if tf != gene]
        if not regs:
            b[gi] = float(y[:, gi].mean())
            continue
        Xg = Xc[:, regs]
        yg = y[:, gi] - y[:, gi].mean()
        A = Xg.T @ Xg + ridge * np.eye(len(regs))
        w = np.linalg.solve(A, Xg.T @ yg)
        W[gi, regs] = w
        b[gi] = float(y[:, gi].mean() - S_tf[:, regs].mean(0) @ w)

    def velocity(W_use):
        a = np.clip(np.log1p(np.exp(np.clip(S_tf @ W_use.T + b, -30, 30))), 0, 50)
        return (a - gamma * Ms).astype("float32")           # beta = 1

    return {"W": W, "b": b, "gamma": gamma, "tfs": tfs, "velocity": velocity}


def _fate_probs(adata, vkey: str, terminal, n_states=None, cluster_key=None):
    """CellRank: VelocityKernel -> GPCCA -> terminal states -> fate probabilities.

    ★ API 核对(cellrank 2.3.2):`predict_terminal_states` 是正确方法名;
      `set_terminal_states_from_macrostates` **在 2.x 中不存在**(旧版本模块的报错来源)。
      调用序列接地上游 tools/_TFScanning_func.py 第 61-86 / 111-121 行。
    """
    import cellrank as cr
    import pandas as pd
    import numpy as np

    vk = cr.kernels.VelocityKernel(adata, vkey=vkey)
    vk.compute_transition_matrix(show_progress_bar=False)
    g = cr.estimators.GPCCA(vk)
    if terminal is None:
        g.compute_macrostates(n_states=n_states, n_cells=30, cluster_key=cluster_key)
        g.predict_terminal_states()                       # 上游 _TFScanning_func.py:79
        names = g.terminal_states.cat.categories.tolist()
        g.set_terminal_states(names)
        ts = g.terminal_states
        terminal = {ct: ts[ts == ct].index.tolist() for ct in names}
    else:
        g.set_terminal_states(terminal)
    g.compute_fate_probabilities(solver="direct", use_petsc=False, show_progress_bar=False)
    fp = pd.DataFrame(np.asarray(g.fate_probabilities),
                      index=adata.obs_names.tolist(),
                      columns=[str(c) for c in g.fate_probabilities.names])
    return fp, terminal


def abundance_test(prob_raw, prob_pert, method: str = "likelihood"):
    """上游 metrics/_abundance_test.py 的忠实复刻(ranksums + ROC-AUC + BH-FDR)。

    coefficient = ROC-AUC(原始 vs 扰动的 fate probability);0.5 = 无效应,
    >0.5 = 该终端命运在扰动后被削弱(与上游同向定义)。
    """
    import numpy as np
    import pandas as pd
    from scipy.stats import ranksums, ttest_ind
    from sklearn.metrics import roc_auc_score
    from statsmodels.stats.multitest import multipletests

    y = np.array([1] * prob_raw.shape[0] + [0] * prob_pert.shape[0])
    X = pd.concat([prob_raw, prob_pert], axis=0)
    rows = []
    for i in range(prob_raw.shape[1]):
        pred = np.asarray(X.iloc[:, i], dtype=float)
        if np.sum(pred) == 0:
            rows.append([np.nan, np.nan])
            continue
        pval = ranksums(pred[y == 0], pred[y == 1], alternative="less")[1]
        score = (roc_auc_score(y, pred) if method == "likelihood"
                 else ttest_ind(pred[y == 0], pred[y == 1])[0])
        rows.append([score, pval])
    tab = pd.DataFrame(rows, index=prob_raw.columns, columns=["coefficient", "p-value"])
    ok = tab["p-value"].notna()
    fdr = np.full(len(tab), np.nan)
    if ok.any():
        fdr[ok.values] = multipletests(tab.loc[ok, "p-value"], method="fdr_bh")[1]
    tab["FDR adjusted p-value"] = fdr
    return tab


def tf_ko_screen(adata, fit, outdir: str, n_states: int = P_N_STATES):
    """逐 TF 把 W 中该列清零 → 重算 velocity → CellRank fate → abundance_test。

    清零逻辑接地 tools/_in_silico_block_simulation.py(`perturb_GRN[row_mask, tf_mask] = effects`)
    与 tools/_TFScanning_func.py(逐 TF 循环 + 固定终端细胞集 + abundance_test)。
    """
    import numpy as np
    import pandas as pd
    import scvelo as scv

    W, tfs = fit["W"], fit["tfs"]
    adata.layers["velocity_grn"] = fit["velocity"](W)
    scv.tl.velocity_graph(adata, vkey="velocity_grn", n_jobs=1)
    fp0, terminal = _fate_probs(adata, "velocity_grn", None,
                                n_states=n_states, cluster_key="cell_state")
    print(f"       terminal states: {list(terminal)} "
          f"(n cells: {[len(v) for v in terminal.values()]})")

    recs = []
    for tf in tfs:
        Wko = W.copy()
        Wko[:, tfs.index(tf)] = 0.0                    # 静默该 TF 的整个 regulon
        adata.layers["velocity_ko"] = fit["velocity"](Wko)
        scv.tl.velocity_graph(adata, vkey="velocity_ko", n_jobs=1)
        fp1, _ = _fate_probs(adata, "velocity_ko", terminal)
        fp1 = fp1.reindex(columns=fp0.columns).fillna(0.0)
        fp1.index = [i + "_perturb" for i in fp1.index]
        res = abundance_test(fp0, fp1)
        for state in fp0.columns:
            recs.append({"TF": tf, "terminal_state": state,
                         "coefficient": float(res.loc[state, "coefficient"]),
                         "pvalue": float(res.loc[state, "p-value"]),
                         "FDR": float(res.loc[state, "FDR adjusted p-value"]),
                         "fate_prob_wt": float(fp0[state].mean()),
                         "fate_prob_ko": float(fp1[state].mean())})
        print(f"       KO {tf}: " + "  ".join(
            f"{s} {fp0[s].mean():.3f}->{fp1[s].mean():.3f}" for s in fp0.columns))

    df = pd.DataFrame(recs)
    df.to_csv(os.path.join(outdir, "tf_ko_fate_shift.csv"), index=False)
    fp0.to_csv(os.path.join(outdir, "fate_probabilities_wt.csv"))
    return df, fp0, terminal


# ============================================================ ④ RegVelo 官方路径
def run_regvelo(adata, skeleton=None):
    """RegVelo 官方路径。守卫式:未装包 / 无 GPU 时不跑,只报明确原因,绝不静默降级。

    以下 API 全部核对自本地克隆源码(Desktop/upstream-sources/561_regvelo):
      pp.set_prior_grn(adata, gt_net, keep_dim=False, cor_filter=True)   preprocessing/_set_prior_grn.py:7
      pp.sanity_check(adata)                                            preprocessing/_sanity_check.py:7
      pp.preprocess_data(adata, spliced_layer="Ms", unspliced_layer="Mu",
                         min_max_scale=True, filter_on_r2=True)          preprocessing/_preprocess_data.py:8
      REGVELOVI.setup_anndata(adata, spliced_layer=, unspliced_layer=)   _model.py:778
      REGVELOVI(adata, W=, regulators=, soft_constraint=True, lam=1,
                lam2=0, **model_kwargs)                                  _model.py:100
      .train(max_epochs=1500, lr=1e-2, ...)                              _model.py:208
      tl.set_output(adata, vae, n_samples=30, batch_size=None)           tools/_set_output.py:11
      tl.in_silico_block_simulation(model, adata, TF, effects=0,
                                    cutoff=1e-3, ...)                    tools/_in_silico_block_simulation.py:6
      tl.TFScanning_func(model, adata, cluster_label=, terminal_states=,
                         KO_list=, n_states=, method="likelihood")       tools/_TFScanning_func.py:14
      mt.abundance_test(prob_raw, prob_pert, method)                     metrics/_abundance_test.py:10
    """
    try:
        import regvelo as rgv
        import torch
    except ImportError as e:
        return {"status": "skipped",
                "reason": f"regvelo not importable ({getattr(e, 'name', e)}); pip install regvelo"}
    if not torch.cuda.is_available():
        return {"status": "skipped",
                "reason": "no CUDA GPU — REGVELOVI.train(max_epochs=1500) is impractical on CPU"}
    exported = [n for n in ("REGVELOVI", "VELOVAE", "ModelComparison") if hasattr(rgv, n)]
    tools = [n for n in ("in_silico_block_simulation", "TFScanning_func", "set_output",
                         "inferred_grn", "perturbation_effect") if hasattr(rgv.tl, n)]
    return {"status": "ready", "regvelo_version": getattr(rgv, "__version__", "?"),
            "exported": exported, "tools": tools,
            "next": "pp.set_prior_grn -> pp.sanity_check -> pp.preprocess_data -> "
                    "REGVELOVI.setup_anndata -> REGVELOVI(adata, W=, regulators=).train() -> "
                    "tl.set_output -> tl.TFScanning_func (regulon KO screen via CellRank)",
            "skeleton_supplied": skeleton is not None}


# ============================================================ ⑤ 出图
def make_figures(adata, ko_df):
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from pubstyle import set_pub_style, save_fig, pal

    set_pub_style(base_size=11)
    os.makedirs(ASSETS, exist_ok=True)
    cols = pal(3, "npg")
    emb = np.asarray(adata.obsm["X_traj"])

    # --- fig1: GRN-informed velocity stream on UMAP -------------------------
    import scvelo as scv
    scv.tl.velocity_embedding(adata, basis="traj", vkey="velocity_grn")
    V = np.asarray(adata.obsm["velocity_grn_traj"])
    fig, ax = plt.subplots(figsize=(5.4, 4.6))
    states = list(adata.obs["cell_state"].cat.categories)
    for c, st in zip(cols, states):
        m = (adata.obs["cell_state"] == st).values
        ax.scatter(emb[m, 0], emb[m, 1], s=16, c=c, lw=0, alpha=.75, label=st)
    gx = np.linspace(emb[:, 0].min(), emb[:, 0].max(), 34)
    gy = np.linspace(emb[:, 1].min(), emb[:, 1].max(), 34)
    GX, GY = np.meshgrid(gx, gy)
    span = max(np.ptp(emb[:, 0]), np.ptp(emb[:, 1]))
    sig = 0.05 * span
    d2 = ((GX.ravel()[:, None] - emb[:, 0][None]) ** 2 +
          (GY.ravel()[:, None] - emb[:, 1][None]) ** 2)
    Wg = np.exp(-d2 / (2 * sig ** 2))
    dens = Wg.sum(1)
    UU = (Wg @ V[:, 0]) / np.maximum(dens, 1e-9)
    VV = (Wg @ V[:, 1]) / np.maximum(dens, 1e-9)
    # ★ 只在**真有细胞**的地方画流线:网格点到最近细胞的距离超过阈值就置 NaN,
    #   否则 streamplot 会在 UMAP 的空白区外推出并不存在的箭头(误导)。
    nn_dist = np.sqrt(d2.min(axis=1))
    off_manifold = nn_dist > 0.045 * span
    UU[off_manifold] = np.nan
    VV[off_manifold] = np.nan
    ax.streamplot(GX, GY, UU.reshape(GX.shape), VV.reshape(GX.shape),
                  color="0.15", linewidth=.7, density=1.3, arrowsize=.9)
    ax.set_xlabel("Diffusion component 1"); ax.set_ylabel("Diffusion component 2")
    ax.set_title("GRN-informed velocity field")
    ax.legend(frameon=False, fontsize=8, loc="best")
    save_fig(fig, os.path.join(ASSETS, "fig1_grn_velocity_stream")); plt.close(fig)

    # --- fig2: TF KO -> fate probability dumbbell ---------------------------
    tstates = sorted(ko_df["terminal_state"].unique())
    tfs = sorted(ko_df["TF"].unique())
    fig, axes = plt.subplots(1, len(tstates), figsize=(3.2 * len(tstates), 3.8), sharey=True)
    axes = np.atleast_1d(axes)
    for ax, st in zip(axes, tstates):
        sub = ko_df[ko_df["terminal_state"] == st].set_index("TF").loc[tfs]
        yy = np.arange(len(tfs))
        for i, tf in enumerate(tfs):
            ax.plot([sub.loc[tf, "fate_prob_wt"], sub.loc[tf, "fate_prob_ko"]],
                    [i, i], color="0.75", lw=1.8, zorder=1)
        ax.scatter(sub["fate_prob_wt"], yy, s=46, color="0.35", zorder=3)
        sig_m = (sub["FDR"] < 0.05).values
        ax.scatter(sub["fate_prob_ko"], yy, s=54, zorder=3,
                   color=[cols[0] if s else cols[2] for s in sig_m],
                   edgecolor="k", linewidth=[0.9 if s else 0.0 for s in sig_m])
        ax.set_yticks(yy); ax.set_yticklabels(tfs)
        ax.set_ylim(-0.8, len(tfs) - 0.2)
        ax.set_xlim(-0.05, 1.05)
        ax.set_xlabel("mean fate probability")
        ax.set_title(st)
    handles = [Line2D([], [], ls="", marker="o", ms=7, color="0.35", label="wild type"),
               Line2D([], [], ls="", marker="o", ms=7, color=cols[0], mec="k", mew=.9,
                      label="regulon KO, FDR < 0.05"),
               Line2D([], [], ls="", marker="o", ms=7, color=cols[2], label="regulon KO, n.s.")]
    # 图例放**图外底部**,避免压住数据点
    fig.legend(handles=handles, frameon=False, fontsize=9, ncol=3,
               loc="lower center", bbox_to_anchor=(0.5, -0.035))
    fig.suptitle("In-silico regulon knockout shifts cell-fate probability", y=1.0, fontsize=11)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    save_fig(fig, os.path.join(ASSETS, "fig2_tf_ko_fate_dumbbell")); plt.close(fig)

    # --- fig3: TF x terminal-state effect heatmap ---------------------------
    M = ko_df.pivot(index="TF", columns="terminal_state", values="coefficient").loc[tfs, tstates]
    F = ko_df.pivot(index="TF", columns="terminal_state", values="FDR").loc[tfs, tstates]
    fig, ax = plt.subplots(figsize=(1.6 * len(tstates) + 2.4, 0.46 * len(tfs) + 2.0))
    v = max(float(np.nanmax(np.abs(M.values - 0.5))), 1e-3)
    im = ax.imshow(M.values, cmap="RdBu_r", vmin=0.5 - v, vmax=0.5 + v, aspect="auto")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if F.values[i, j] < 0.05:
                ax.scatter(j, i, marker="*", s=80, color="k", zorder=3)
    ax.set_xticks(range(len(tstates))); ax.set_xticklabels(tstates, rotation=20, ha="right")
    ax.set_yticks(range(len(tfs))); ax.set_yticklabels(tfs)
    ax.set_title("Regulon-KO effect (AUROC)\n$\\star$ FDR < 0.05")
    cb = fig.colorbar(im, ax=ax, shrink=.85)
    cb.set_label("AUROC  (0.5 = no effect)")
    fig.tight_layout()
    save_fig(fig, os.path.join(ASSETS, "fig3_ko_effect_heatmap")); plt.close(fig)

    # --- fig4: honest baseline — GRN-informed vs plain scVelo ---------------
    tt = np.asarray(adata.obs["true_time"])
    st = np.asarray(adata.obs["cell_state"].astype(str))
    panels = [("pt_scvelo", "scVelo (no GRN)"), ("pt_grn", "GRN-informed (RegVelo form)")]
    fig, axes = plt.subplots(1, 2, figsize=(7.8, 3.9), sharey=False)
    for ax, (key, lab) in zip(axes, panels):
        if key not in adata.obs:
            ax.set_axis_off(); continue
        pt = np.asarray(adata.obs[key], dtype=float)
        for c, s in zip(cols, list(adata.obs["cell_state"].cat.categories)):
            m = st == s
            ax.scatter(tt[m], pt[m], s=13, c=c, lw=0, alpha=.65, label=s)
        r = lineage_rho(pt, tt, st)
        ax.set_xlabel("true latent time (synthetic ground truth)")
        ax.set_title(f"{lab}\nlineage-wise Spearman |rho| = {r['mean']:.3f}", fontsize=10)
    axes[0].set_ylabel("inferred velocity pseudotime")
    axes[0].legend(frameon=False, fontsize=8, loc="upper left")
    fig.suptitle("Velocity pseudotime vs planted ground truth "
                 "(Spearman rho computed within each lineage)",
                 y=1.04, fontsize=9)
    fig.tight_layout()
    save_fig(fig, os.path.join(ASSETS, "fig4_pseudotime_vs_truth")); plt.close(fig)

    return ["fig1_grn_velocity_stream.png", "fig2_tf_ko_fate_dumbbell.png",
            "fig3_ko_effect_heatmap.png", "fig4_pseudotime_vs_truth.png"]


# ============================================================ main
def main():
    ap = argparse.ArgumentParser(description="561 · RegVelo GRN-informed velocity & regulon KO")
    ap.add_argument("--h5ad", default=os.path.join(EXAMPLE, "synthetic_velocity.h5ad"),
                    help="AnnData with layers['spliced'|'unspliced']")
    ap.add_argument("--grn", default=os.path.join(EXAMPLE, "synthetic_prior_grn.csv"),
                    help="prior GRN skeleton csv, rows=targets, cols=regulators")
    ap.add_argument("--outdir", default=RESULTS)
    ap.add_argument("--ridge", type=float, default=P_RIDGE)
    ap.add_argument("--n-states", type=int, default=P_N_STATES)
    ap.add_argument("--run-regvelo", action="store_true",
                    help="attempt the official RegVelo path (needs pip install regvelo + GPU)")
    ap.add_argument("--regen-example", action="store_true", help="force-regenerate example_data/")
    a = ap.parse_args()

    import numpy as np
    np.random.seed(SEED)
    import anndata as ad
    import pandas as pd
    import scvelo as scv
    from scipy.stats import spearmanr

    os.makedirs(a.outdir, exist_ok=True)
    os.makedirs(ASSETS, exist_ok=True)

    if a.regen_example or not (os.path.exists(a.h5ad) and os.path.exists(a.grn)):
        print("[561] Step 0  generating synthetic example data")
        a.h5ad, a.grn = make_example_data(EXAMPLE)

    print("[561] Step 1  loading input")
    adata = ad.read_h5ad(a.h5ad)
    for L in ("spliced", "unspliced"):
        if L not in adata.layers:
            sys.exit(f"input lacks layers['{L}'] — RegVelo and scVelo both need spliced/unspliced")
    skeleton = pd.read_csv(a.grn, index_col=0)
    print(f"       {adata.n_obs} cells x {adata.n_vars} genes; "
          f"prior GRN {skeleton.shape[0]} targets x {skeleton.shape[1]} regulators")

    print("[561] Step 2  honest comparator: plain scVelo (no GRN)")
    base = scvelo_baseline(adata)
    for k, v in base.items():
        print(f"       {k}: {v}")

    print("[561] Step 3  GRN-linear velocity (RegVelo functional form, ridge-fitted W)")
    fit = fit_grn_velocity(adata, skeleton, ridge=a.ridge)
    print(f"       regulators used: {fit['tfs']}")

    print("[561] Step 4  in-silico regulon knockout screen via CellRank")
    ko_df, fp0, terminal = tf_ko_screen(adata, fit, a.outdir, n_states=a.n_states)

    try:
        scv.tl.velocity_pseudotime(adata, vkey="velocity_grn")
        adata.obs["pt_grn"] = adata.obs["velocity_grn_pseudotime"].values
        base["grn_lineage_rho"] = lineage_rho(
            adata.obs["pt_grn"], adata.obs["true_time"], adata.obs["cell_state"].astype(str))
        print(f"       scvelo_lineage_rho: {base.get('scvelo_lineage_rho')}")
        print(f"       grn_lineage_rho:    {base['grn_lineage_rho']}")
    except Exception as e:
        print(f"       GRN pseudotime skipped: {type(e).__name__}")

    print("[561] Step 5  sanity-check (positive control: planted Branch_A drivers TF2/TF3)")
    sanity = {}
    if "Branch_A" in set(ko_df["terminal_state"]):
        sub = ko_df[ko_df["terminal_state"] == "Branch_A"].set_index("TF")
        drivers = [t for t in ("TF2", "TF3") if t in sub.index]
        others = [t for t in sub.index if t not in drivers]
        if drivers and others:
            da = float(np.mean([sub.loc[t, "fate_prob_wt"] - sub.loc[t, "fate_prob_ko"] for t in drivers]))
            do = float(np.mean([sub.loc[t, "fate_prob_wt"] - sub.loc[t, "fate_prob_ko"] for t in others]))
            sanity = {"Branch_A_loss_by_planted_drivers": round(da, 4),
                      "Branch_A_loss_by_other_TFs": round(do, 4),
                      "pipeline_recovers_planted_drivers": bool(da > do)}
            for k, v in sanity.items():
                print(f"       {k}: {v}")

    reg = None
    if a.run_regvelo:
        print("[561] Step 6  official RegVelo path")
        reg = run_regvelo(adata, skeleton)
        for k, v in reg.items():
            print(f"       {k}: {v}")
    else:
        print("[561] Step 6  official RegVelo path not requested (--run-regvelo); proxy only")

    print("[561] Step 7  figures")
    figs = make_figures(adata, ko_df)
    for f in figs:
        print(f"       assets/{f}")

    top = ko_df.assign(effect=(ko_df["coefficient"] - 0.5).abs()) \
               .sort_values("effect", ascending=False).head(5)
    summary = {
        "module": "561_regvelo_grn_velocity",
        "n_cells": int(adata.n_obs), "n_genes": int(adata.n_vars),
        "regulators": fit["tfs"],
        "terminal_states": {k: len(v) for k, v in terminal.items()},
        "comparator": base,
        "sanity_check": sanity,
        "top_ko_effects": top[["TF", "terminal_state", "coefficient", "FDR"]].to_dict("records"),
        "regvelo_official": reg,
        "figures": figs,
        "seed": SEED,
        "session": _session_info(),
    }
    with open(os.path.join(a.outdir, "561_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=1, ensure_ascii=False, default=str)
    print(f"[561] done -> {os.path.join(a.outdir, '561_summary.json')}")


if __name__ == "__main__":
    main()
