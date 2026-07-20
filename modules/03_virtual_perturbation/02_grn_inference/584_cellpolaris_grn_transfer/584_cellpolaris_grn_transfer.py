# -*- coding: utf-8 -*-
# =============================================================================
# 584 · CellPolaris —— GRN 迁移 + 概率图模型 TF 敲除模拟 + 主控 TF 打分
# -----------------------------------------------------------------------------
# 上游   : Feng G, et al. CellPolaris: Transfer Learning for Gene Regulatory
#          Network Construction to Guide Cell State Transitions.
#          Adv Sci (Weinh). 2026 Feb;13(12):e08697.
#          PMID 41498638 · DOI 10.1002/advs.202508697 · PMC12948241
#          仓库 https://github.com/xCompass-AI/CellPolaris
# -----------------------------------------------------------------------------
# 本模块做什么
#   CellPolaris 上游有三段:
#     (1) transfer_learning/  用 NCF link-prediction + 域迁移损失,从 PECA2 GRN 库
#         迁移出新语境的 GRN。需要 sci-db 下载的 PECA2 数据集 + 训练好的 .pt 权重
#         + torch_geometric + GPU —— 三者仓库里都没有 → 本模块做「守卫式封装」,
#         只检查环境、打印真实命令,绝不假装能跑。
#     (2) model_PGM/          在 GRN 上建高斯概率图模型,拟合边强度 k,再把 TF 置零
#         推出邻居基因表达变化 deltaX —— 本模块**忠实向量化重写**,本机 CPU 可跑。
#     (3) plot/CellCruise.py  用 deltaX 与真实相邻细胞状态表达差的余弦相似度,给 TF
#         在每段轨迹上打分,排出主控 TF —— 本模块忠实重写(不含 velocyto 箭头图部分)。
#
# 为什么是「重写」而不是「调用上游」(诚实说明,已逐行核对源码):
#   · model_PGM/model_PGM.py:280+ 的 __main__ 里硬编码了作者机器上的绝对路径
#     '/home/share/jingzi_sample/RS1o2_sample.csv',换机即崩;
#   · model_PGM/model_PGM.py:217-221 `Model.forward` 末尾对空列表 alpha_list 迭代后
#     `return -(p + q + k + g)`,变量 g 从未定义 → 原样运行会 NameError;
#     同文件 222 行被注释掉的 `return -(p + q + k)` 才是可运行形式,本模块采用它;
#   · model_PGM/model_PGM.py:167-168 `Gauss_condition.forward` 的 TF_high 分支引用了
#     该分支未定义的 gTF_mu / gTF_sigma(只定义了 gTF_high_mu / gTF_high_sigma)
#     → 同样 NameError。本模块按显然的意图用 gTF_high_mu / gTF_high_sigma。
#   这些都是上游源码的实际状态,不是猜测;数学形式(边界条件、relu、eps=0.01、
#   father_num-1 权重、-5*k*alpha*mu/sigma^2 的 deltaX)一律照抄。
#
# 依赖   : numpy pandas scipy matplotlib torch(本机均已具备,CPU 即可)
# 状态   : 🟡 PGM + 主控 TF 打分本机零改动跑通出图;GRN 迁移段需上游数据集+权重+GPU
# =============================================================================
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch

ROOT = Path(__file__).resolve().parent
FRAMEWORK = ROOT.parents[2] / "_framework"
sys.path.insert(0, str(FRAMEWORK))
from pubstyle import (  # noqa: E402
    CMAP_DIVERGE,
    NATURE_W1,
    NATURE_W2,
    pal,
    save_fig,
    set_pub_style,
)

import matplotlib.pyplot as plt  # noqa: E402

EPS = 0.01  # 上游 model_PGM.py:313 eps = 0.01


# =============================================================================
# 0. 合成示例数据(synthetic, for demo only)
# =============================================================================
def simulate_example_data(outdir: Path, seed: int = 0, n_metacell: int = 50) -> None:
    """造一套三状态分化轨迹的小型合成数据,文件格式严格对齐上游 model_PGM/data/。

    设计:两层 GRN。root TF(只当 TF 不当 TG)→ mid TF → target gene。
    真值:DriverA 驱动 S1→S2 的上调模块,DriverB 驱动 S2→S3 的上调模块,
    Decoy1/Decoy2 各带一个在轨迹上不变的模块。主控 TF 打分应把 Driver 排出来。
    """
    rng = np.random.default_rng(seed)
    outdir.mkdir(parents=True, exist_ok=True)

    root_tfs = ["DriverA", "DriverB", "Decoy1", "Decoy2"]
    mid_tfs = [f"MidTF{i}" for i in range(1, 5)]
    targets = [f"Gene{i:02d}" for i in range(1, 25)]

    # root → mid:每个 root 带一个 mid TF
    edges = [(r, m) for r, m in zip(root_tfs, mid_tfs)]
    # root → target 与 mid → target:每个模块 6 个下游基因
    module = {r: targets[i * 6:(i + 1) * 6] for i, r in enumerate(root_tfs)}
    for r, m in zip(root_tfs, mid_tfs):
        for g in module[r]:
            edges.append((r, g))
            edges.append((m, g))

    edge_w = {e: rng.uniform(0.35, 0.9) for e in edges}          # 真值耦合强度
    edge_score = {e: rng.uniform(0.55, 0.95) for e in edges}     # 模拟迁移模型给的置信度

    # 轨迹上的 root TF 活性:DriverA 在 S1→S2 升,DriverB 在 S2→S3 升,Decoy 不变
    activity = {
        "DriverA": {"S1": 18.0, "S2": 34.0, "S3": 34.0},
        "DriverB": {"S1": 16.0, "S2": 16.0, "S3": 32.0},
        "Decoy1":  {"S1": 25.0, "S2": 25.0, "S3": 25.0},
        "Decoy2":  {"S1": 22.0, "S2": 22.0, "S3": 22.0},
    }
    base = {g: rng.uniform(20.0, 45.0) for g in mid_tfs + targets}

    states = ["S1", "S2", "S3"]
    for st in states:
        # root TF 的 metacell 层面活性(带跨 metacell 波动 → 供 PGM 估协方差)
        a = {r: rng.normal(activity[r][st], 3.0, n_metacell) for r in root_tfs}
        expr = {r: np.clip(a[r], 0.1, None) for r in root_tfs}
        for r, m in zip(root_tfs, mid_tfs):
            expr[m] = np.clip(base[m] + edge_w[(r, m)] * (a[r] - activity[r][st])
                              + rng.normal(0, 2.0, n_metacell) + 0.35 * a[r], 0.1, None)
        for r, m in zip(root_tfs, mid_tfs):
            for g in module[r]:
                expr[g] = np.clip(
                    base[g]
                    + edge_w[(r, g)] * a[r]
                    + edge_w[(m, g)] * (expr[m] - expr[m].mean())
                    + rng.normal(0, 2.0, n_metacell), 0.1, None)

        genes = root_tfs + mid_tfs + targets
        mat = pd.DataFrame({g: expr[g] for g in genes}).T
        mat.columns = [f"MC{j:02d}" for j in range(n_metacell)]
        mat.index.name = "Row.names"
        mc_path = outdir / f"expr_metacell_{st}.csv"
        with open(mc_path, "w", encoding="utf-8", newline="") as fh:
            fh.write("# synthetic, for demo only -- NOT real scRNA-seq metacells\n")
            fh.write("# rows = genes, columns = metacells (pseudobulked cell groups)\n")
            mat.to_csv(fh)

        # 伪 bulk:上游 data/RS*.txt 格式 = 两列无表头 gene<TAB>value
        pb = mat.mean(axis=1)
        with open(outdir / f"expr_pseudobulk_{st}.txt", "w", encoding="utf-8", newline="") as fh:
            for g, v in pb.items():
                fh.write(f"{g}\t{v:.6f}\n")

    # GRN:上游 generated_grn 格式 = TF/TG/Score 制表符,按 Score 降序
    grn = pd.DataFrame(
        [{"TF": t, "TG": g, "Score": edge_score[(t, g)]} for t, g in edges]
    ).sort_values("Score", ascending=False)
    grn.to_csv(outdir / "grn.txt", sep="\t", index=False)

    (outdir / "trajectory.txt").write_text(
        "# synthetic, for demo only -- one trajectory per line, states tab-separated\n"
        "S1\tS2\tS3\n", encoding="utf-8")


def ensure_example_data(outdir: Path, seed: int = 0) -> None:
    need = ["grn.txt", "trajectory.txt"] + [
        f"expr_{k}_{s}.{ext}" for s in ("S1", "S2", "S3")
        for k, ext in (("metacell", "csv"), ("pseudobulk", "txt"))]
    if not all((outdir / f).exists() for f in need):
        print("[example_data] 缺文件,重新生成合成示例数据 ...")
        simulate_example_data(outdir, seed=seed)


# =============================================================================
# 1. 数据准备 —— 忠实对应上游 model_PGM/model_PGM.py:23 data_pro1()
# =============================================================================
def prepare(grn_file: Path, metacell_csv: Path, edge_num: int,
            n_columns: int = 50, min_std: float = 1.0):
    """返回 (root_tfs, tgs, father_num, edge_table, X_root, X_tg)。

    与上游 data_pro1 的对应关系(逐行核对 model_PGM.py:23-97):
      · 只保留 TF 与 TG 都出现在表达矩阵里的边;
      · 取表达矩阵前 n_columns 列(上游写死 columns = 50);
      · 删掉跨列标准差 < min_std 的基因(上游 gene_std < 1);
      · 去自环;取前 edge_num 条边;
      · Score 线性缩放到 [0.1, 0.7](上游 (x-min)/(max-min)*0.6 + 0.1);
      · root TF(上游叫 top_TF) = set(TF) - set(TG),即入度为 0 的调控者;
      · father_num = 每个 TG 的入度;
      · 每条边算 TF 与 TG 在这 n_columns 列上的协方差 Cov。
    唯一偏离:上游 model_PGM.py:30/44 算了 network_sorted 却切了未排序的 network
    (靠输入文件本身已按 Score 降序才没出错),本模块显式按 Score 降序后再切。
    """
    net = pd.read_table(grn_file)[["TF", "TG", "Score"]]
    expr = pd.read_csv(metacell_csv, index_col=0, comment="#")
    expr = expr.iloc[:, :n_columns]
    expr = expr.groupby(expr.index).mean()

    net = net[net["TF"].isin(expr.index) & net["TG"].isin(expr.index)]
    gene_std = expr.std(axis=1, ddof=0)
    drop = set(gene_std[gene_std < min_std].index)
    net = net[~net["TF"].isin(drop) & ~net["TG"].isin(drop)]
    net = net[net["TF"] != net["TG"]]
    net = net.sort_values("Score", ascending=False).iloc[:edge_num].reset_index(drop=True)
    if net.empty:
        raise ValueError(f"{metacell_csv.name}: 过滤后没有剩下任何边,请放宽 --min-std / --edge-num")

    lo, hi = net["Score"].min(), net["Score"].max()
    net["Alpha"] = ((net["Score"] - lo) / (hi - lo) * 0.6 + 0.1) if hi > lo else 0.4

    root_tfs = sorted(set(net["TF"]) - set(net["TG"]))
    tgs = sorted(set(net["TG"]))
    father_num = net["TG"].value_counts().reindex(tgs).to_numpy(dtype=float)

    X_root = expr.loc[root_tfs].to_numpy(dtype=np.float64).T   # (n_metacell, n_root)
    X_tg = expr.loc[tgs].to_numpy(dtype=np.float64).T          # (n_metacell, n_tg)

    tf_vals = expr.loc[net["TF"]].to_numpy(dtype=np.float64)
    tg_vals = expr.loc[net["TG"]].to_numpy(dtype=np.float64)
    tf_c = tf_vals - tf_vals.mean(axis=1, keepdims=True)
    tg_c = tg_vals - tg_vals.mean(axis=1, keepdims=True)
    net["Cov"] = (tf_c * tg_c).sum(axis=1) / (tf_vals.shape[1] - 1)

    return root_tfs, tgs, father_num, net, X_root, X_tg


# =============================================================================
# 2. 高斯概率图模型 —— 忠实对应 model_PGM.py:113-222 的 Gauss_* / Model
# =============================================================================
def fit_pgm(root_tfs, tgs, father_num, net, X_root, X_tg, *,
            n_iter=49, lr_mu=0.01, lr_sigma=0.01, lr_k=0.1,
            mu_clamp=0.9, sigma_lo=0.8, sigma_hi=1.2, clip_grad=1.0, seed=0):
    """在 GRN 上拟合高斯 PGM,返回 (mu, sigma, k, loss 曲线)。

    目标函数(model_PGM.py:204-222,采用其可运行分支 `-(p + q + k)`):
        p = Σ_rootTF   logN(x_TF | mu_TF, sigma_TF)
        q = -Σ_TG (indeg-1) · logN(x_TG | mu_TG, sigma_TG)     # 扣掉重复计入的边际
        k = Σ_edge     logN(x_TG | loc, scale)
        loc   = relu(mu_TG + k_e · Cov_e · (x_TF - mu_TF) / sigma_TF²) + 0.01
        scale = sqrt(relu(sigma_TG² - Alpha_e² / sigma_TF²) + 0.01)
    优化器与学习率(model_PGM.py:351-371):SGD,mu/sigma 0.01,k 0.1,k 初值 0。
    每步后按 model_PGM.py:396-403 夹逼:mu → 0.9·初值(上下界相同 = 冻结),
    sigma → [0.8, 1.2]·初值。这两条是上游行为,保留但暴露成参数。
    梯度裁剪对应 model_PGM.py:393-394(上游写了但注释掉);本模块默认开启,
    因为原始表达尺度下 k 的梯度很容易发散,关掉请传 --clip-grad 0。
    """
    torch.manual_seed(seed)
    n_root, n_tg = len(root_tfs), len(tgs)
    ridx = {g: i for i, g in enumerate(root_tfs)}
    tidx = {g: i + n_root for i, g in enumerate(tgs)}

    init_mu = np.concatenate([X_root.mean(0), X_tg.mean(0)])
    init_sd = np.concatenate([X_root.std(0), X_tg.std(0)])
    init_sd = np.clip(init_sd, 1e-3, None)

    MU = torch.nn.Parameter(torch.tensor(init_mu, dtype=torch.float32))
    SIGMA = torch.nn.Parameter(torch.tensor(np.abs(init_sd), dtype=torch.float32))
    K = torch.nn.Parameter(torch.zeros(len(net), dtype=torch.float32))
    mu0, sd0 = MU.detach().clone(), SIGMA.detach().clone()

    e_tf = torch.tensor([ridx.get(t, tidx.get(t)) for t in net["TF"]], dtype=torch.long)
    e_tg = torch.tensor([tidx[t] for t in net["TG"]], dtype=torch.long)
    alpha = torch.tensor(net["Alpha"].to_numpy(), dtype=torch.float32)
    cov = torch.tensor(net["Cov"].to_numpy(), dtype=torch.float32)

    Xall = torch.tensor(np.concatenate([X_root, X_tg], axis=1), dtype=torch.float32)
    x_root_t = torch.tensor(X_root, dtype=torch.float32)
    x_tg_t = torch.tensor(X_tg, dtype=torch.float32)
    fnum = torch.tensor(father_num - 1.0, dtype=torch.float32)

    opt = torch.optim.SGD([{"params": [MU], "lr": lr_mu},
                           {"params": [SIGMA], "lr": lr_sigma},
                           {"params": [K], "lr": lr_k}])

    losses = []
    for it in range(n_iter):
        opt.zero_grad()
        p = torch.distributions.Normal(MU[:n_root], SIGMA[:n_root]).log_prob(x_root_t).sum()
        q = -(fnum * torch.distributions.Normal(MU[n_root:], SIGMA[n_root:])
              .log_prob(x_tg_t)).sum()
        x_tf_e, x_tg_e = Xall[:, e_tf], Xall[:, e_tg]
        loc = torch.relu(MU[e_tg] + K * cov * (x_tf_e - MU[e_tf])
                         / torch.square(SIGMA[e_tf])) + EPS
        scale = torch.sqrt(torch.relu(torch.square(SIGMA[e_tg])
                                      - torch.square(alpha) / torch.square(SIGMA[e_tf])) + EPS)
        kterm = torch.distributions.Normal(loc, scale).log_prob(x_tg_e).sum()

        loss = -(p + q + kterm)
        loss.backward()
        if clip_grad and clip_grad > 0:
            torch.nn.utils.clip_grad_norm_([MU, SIGMA, K], clip_grad)
        opt.step()
        with torch.no_grad():
            MU.copy_(torch.clamp(MU, mu_clamp * mu0, mu_clamp * mu0))
            SIGMA.copy_(torch.clamp(SIGMA, sigma_lo * sd0, sigma_hi * sd0))
        losses.append(float(loss.detach()))
        if it % 10 == 0 or it == n_iter - 1:
            print(f"    iter {it:>3d}  NLL = {losses[-1]:.2f}")

    return (MU.detach().numpy(), SIGMA.detach().numpy(), K.detach().numpy(),
            np.array(losses), {**ridx, **tidx})


# =============================================================================
# 3. deltaX —— 忠实对应 model_PGM.py:247-275 get_param_numpy()
# =============================================================================
def compute_deltax(net, mu, sigma, k, gidx, ko_factor=5.0) -> pd.DataFrame:
    """TF 置零后邻居 TG 的表达变化量。

    上游公式(model_PGM.py:272):
        delta_x = -5 * k * weight_value * tf_u / (tf_sigma ** 2)
    其中 weight_value 是该边缩放后的 Score(本模块列名 Alpha),
    tf_u / tf_sigma 是 TF 的高斯边际参数。系数 5 上游写死,这里暴露为 --ko-factor。
    """
    tf_i = np.array([gidx[t] for t in net["TF"]])
    dx = -ko_factor * k * net["Alpha"].to_numpy() * mu[tf_i] / np.square(sigma[tf_i])
    out = pd.DataFrame({"TF": net["TF"].to_numpy(), "TG": net["TG"].to_numpy(), "deltaX": dx})
    wide = out.pivot_table(index="TG", columns="TF", values="deltaX", aggfunc="sum")
    return wide.fillna(0.0)


# =============================================================================
# 4. 主控 TF 打分 —— 忠实对应 plot/CellCruise.py:900 identifying_important_TFs()
# =============================================================================
def master_tf_scores(deltax_by_state, pseudobulk, trajectories) -> pd.DataFrame:
    """每段状态转换上,算「观测表达差」与「模拟敲除 deltaX」的余弦相似度。

    上游 CellCruise.py:954-969:
        data1 = fpkm[next] - fpkm[start]      观测到的相邻细胞状态表达差
        data2 = deltaX[start]                 在 start 状态敲除该 TF 的模拟位移
        只保留 data2 != 0 的基因(即该 TF 真正连到的 TG),再算余弦相似度。
    读法:得分为负 = 敲除把细胞往转换的反方向推,该 TF 是这段转换的推动者;
    得分为正 = 敲除反而顺着转换推。本模块按得分排序并保留符号,不做二元判定。

    额外报一列 obs_shift_norm = ||data1||(该 TF 靶基因上真实发生的表达变化幅度)。
    余弦相似度是尺度无关的,靶基因在轨迹上根本没动的 TF 也能靠噪声拿到 |cos|≈1
    (本模块 demo 里 Decoy1 就是这样)。这一列不改变上游打分,只是让"方向对但幅度
    可忽略"的假阳性能被看出来;Fig3 用点的大小编码它。
    """
    rows = []
    for traj in trajectories:
        for start, nxt in zip(traj[:-1], traj[1:]):
            if start not in deltax_by_state:
                continue
            dxw = deltax_by_state[start]
            d1_full = pseudobulk[nxt] - pseudobulk[start]
            for tf in dxw.columns:
                d2 = dxw[tf]
                d2 = d2[d2 != 0]
                genes = [g for g in d2.index if g in d1_full.index]
                if len(genes) < 2:
                    continue
                a = d1_full.loc[genes].to_numpy(dtype=float)
                b = d2.loc[genes].to_numpy(dtype=float)
                na, nb = np.linalg.norm(a), np.linalg.norm(b)
                if na == 0 or nb == 0:
                    continue
                rows.append({"transition": f"{start}->{nxt}", "TF": tf,
                             "n_targets": len(genes),
                             "cosine": float(a @ b / (na * nb)),
                             "obs_shift_norm": float(na)})
    return pd.DataFrame(rows)


# =============================================================================
# 5. 守卫式封装:上游 transfer_learning 段(GRN 迁移)
# =============================================================================
def check_transfer_learning(repo: Path | None) -> dict:
    """检查能否真的跑上游 GRN 迁移;跑不了就说清缺什么,不静默降级、不假装。

    API 逐条读自本地克隆的上游源码(2026-07-21):
      transfer_learning/generate_external_grn.py     命令行入口与全部参数
      transfer_learning/models/base_model/ncf.py     NCF(num_nodes, hidden_dim=768)
      transfer_learning/models/transfer_model.py     TRModel(link_prediction_model,
                                                       loss_type, num_sourcedomains,
                                                       mixup_alpha, top_ratio, ...)
      transfer_learning/load_dataset/dataset.py:12   GRNPredictionDataset(root="dataset")
      transfer_learning/generate_grn.py              generate_grn / dump_generated_grn
    """
    st = {"available": False, "repo": str(repo) if repo else None, "missing": []}
    if repo is None:
        st["missing"].append("未提供 --cellpolaris-repo(上游仓库路径)")
        return st
    entry = repo / "transfer_learning" / "generate_external_grn.py"
    if not entry.exists():
        st["missing"].append(f"找不到入口 {entry}")
    # GRNPredictionDataset(root="dataset") 是**相对当前工作目录**解析的(dataset.py:12-14),
    # 不是相对仓库根;上游 README 从仓库根启动脚本,故两处都查。
    ds_candidates = [repo / "dataset" / "processed_dataset.pt",
                     repo / "transfer_learning" / "dataset" / "processed_dataset.pt"]
    if not any(p.exists() for p in ds_candidates):
        st["missing"].append(
            "缺 PECA2 预处理数据集 processed_dataset.pt(GRNPredictionDataset 的 root='dataset',"
            f" 按启动目录解析,已查 {[str(p) for p in ds_candidates]};"
            " 需从 https://www.scidb.cn/en/s/VNvY3e 下载原始 PECA2 数据,仓库不含)")
    ckpts = list((repo / "transfer_learning" / "result").rglob("*.pt")) if repo.exists() else []
    if not ckpts:
        st["missing"].append(
            "缺训练好的迁移模型权重 result/multi_to_multi_generalization/**/fold_*.pt"
            "(仓库不含,需自行训练:script_train_multi_to_multi.py)")
    try:
        import torch_geometric  # noqa: F401
    except Exception as e:
        st["missing"].append(f"缺 torch_geometric(上游 environment.yml 锁 2.3.1):{type(e).__name__}")
    st["available"] = not st["missing"]
    return st


# =============================================================================
# 6. 出图(全部走 pubstyle;不用条形图)
# =============================================================================
def fig_pgm_fit(losses_by_state, net_by_state, k_by_state, out):
    set_pub_style(base_size=9)
    fig, axes = plt.subplots(1, 2, figsize=(NATURE_W2, 2.9))
    cols = pal(len(losses_by_state), "npg")

    for c, (st, ls) in zip(cols, losses_by_state.items()):
        axes[0].plot(np.arange(len(ls)), ls, lw=1.6, color=c, label=st)
    axes[0].set_xlabel("SGD iteration")
    axes[0].set_ylabel("Negative log-likelihood")
    axes[0].set_title("PGM optimisation")
    axes[0].legend(title="Cell state")

    for c, st in zip(cols, net_by_state):
        axes[1].scatter(net_by_state[st]["Cov"], k_by_state[st], s=14, alpha=0.75,
                        color=c, edgecolor="none", label=st)
    axes[1].axhline(0, color="grey", lw=0.8, ls="--")
    axes[1].axvline(0, color="grey", lw=0.8, ls="--")
    axes[1].set_xlabel("Edge covariance  Cov(TF, TG)")
    axes[1].set_ylabel("Learned coupling  k")
    axes[1].set_title("Fitted edge coupling")
    axes[1].legend(title="Cell state")

    fig.tight_layout()
    save_fig(fig, out)
    plt.close(fig)


def fig_deltax_heatmap(dxw: pd.DataFrame, state: str, out, top_n=22):
    set_pub_style(base_size=9)
    order = dxw.abs().max(axis=1).sort_values(ascending=False).index[:top_n]
    m = dxw.loc[order]
    v = float(np.abs(m.to_numpy()).max()) or 1.0

    fig, ax = plt.subplots(figsize=(NATURE_W2 * 0.72, 0.19 * len(order) + 1.5))
    im = ax.imshow(m.to_numpy(), aspect="auto", cmap=CMAP_DIVERGE, vmin=-v, vmax=v)
    ax.set_xticks(range(m.shape[1]))
    ax.set_xticklabels(m.columns, rotation=45, ha="right")
    ax.set_yticks(range(m.shape[0]))
    ax.set_yticklabels(m.index)
    ax.set_xlabel("Knocked-out TF")
    ax.set_ylabel("Target gene")
    ax.set_title(f"Simulated KO response  $\\Delta X$  ({state})")
    cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cb.set_label("$\\Delta X$")
    fig.tight_layout()
    save_fig(fig, out)
    plt.close(fig)


def fig_master_tf(scores: pd.DataFrame, out):
    set_pub_style(base_size=9)
    trans = list(dict.fromkeys(scores["transition"]))
    fig, axes = plt.subplots(1, len(trans), figsize=(NATURE_W1 * len(trans), 3.2),
                             squeeze=False)
    cols = pal(2, "npg")
    smax = float(scores["obs_shift_norm"].max()) or 1.0
    for ax, tr in zip(axes[0], trans):
        s = scores[scores["transition"] == tr].sort_values("cosine")
        y = np.arange(len(s))
        c = [cols[0] if v < 0 else cols[1] for v in s["cosine"]]
        sz = 20 + 200 * (s["obs_shift_norm"] / smax)
        ax.hlines(y, 0, s["cosine"], color=c, lw=1.4, alpha=0.85)
        ax.scatter(s["cosine"], y, s=sz, color=c, zorder=3,
                   edgecolor="white", linewidth=0.6)
        ax.set_yticks(y)
        ax.set_yticklabels(s["TF"])
        ax.axvline(0, color="grey", lw=0.8, ls="--")
        ax.set_xlabel("Cosine(observed $\\Delta$Expr, simulated $\\Delta X$)")
        ax.set_title(f"Transition {tr}")
    axes[0][0].set_ylabel("Knocked-out TF")
    # 点大小 = 靶基因上真实表达变化幅度 ||ΔExpr||,用来识破"方向对但没动"的假阳性
    for lab, frac in (("low", 0.15), ("high", 1.0)):
        axes[0][-1].scatter([], [], s=20 + 200 * frac, color="0.45",
                            edgecolor="white", linewidth=0.6,
                            label=f"{lab}  ({frac * smax:.0f})")
    axes[0][-1].legend(title="||observed $\\Delta$Expr||", loc="upper left",
                       bbox_to_anchor=(1.02, 1.0), labelspacing=1.4, borderpad=0.8)
    fig.tight_layout()
    save_fig(fig, out)
    plt.close(fig)


# =============================================================================
# 7. 主流程
# =============================================================================
def main() -> int:
    ap = argparse.ArgumentParser(description="584 · CellPolaris GRN 迁移 / PGM 敲除 / 主控 TF")
    ap.add_argument("--grn", default=None, help="GRN 文件 TF/TG/Score 制表符")
    ap.add_argument("--datadir", default=None, help="example_data 目录")
    ap.add_argument("--outdir", default=None, help="输出目录,默认 results/")
    ap.add_argument("--edge-num", type=int, default=3000, help="保留的最强边数(上游 run.py 用 3000)")
    ap.add_argument("--columns", type=int, default=50, help="用几列 metacell(上游写死 50)")
    ap.add_argument("--min-std", type=float, default=1.0, help="基因跨 metacell 标准差下限")
    ap.add_argument("--iters", type=int, default=49, help="SGD 迭代数(上游 range(1,50))")
    ap.add_argument("--lr-k", type=float, default=0.1)
    ap.add_argument("--lr-mu", type=float, default=0.01)
    ap.add_argument("--lr-sigma", type=float, default=0.01)
    ap.add_argument("--clip-grad", type=float, default=1.0, help="0 表示关闭")
    ap.add_argument("--ko-factor", type=float, default=5.0, help="deltaX 里上游写死的系数 5")
    ap.add_argument("--cellpolaris-repo", default=None, help="上游仓库路径(启用 GRN 迁移段检查)")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    np.random.seed(args.seed)
    torch.manual_seed(args.seed)

    datadir = Path(args.datadir) if args.datadir else ROOT / "example_data"
    outdir = Path(args.outdir) if args.outdir else ROOT / "results"
    outdir.mkdir(parents=True, exist_ok=True)
    assets = ROOT / "assets"
    ensure_example_data(datadir, seed=args.seed)

    grn_file = Path(args.grn) if args.grn else datadir / "grn.txt"
    trajectories = [ln.split("\t") for ln in
                    (datadir / "trajectory.txt").read_text(encoding="utf-8").splitlines()
                    if ln.strip() and not ln.startswith("#")]
    states = list(dict.fromkeys([s for t in trajectories for s in t]))
    print(f"Step 0 · GRN = {grn_file.name};状态 = {states};轨迹 = {trajectories}")

    # --- Step 1-3:逐状态跑 PGM 与 deltaX -------------------------------------
    losses_by_state, net_by_state, k_by_state, dx_by_state = {}, {}, {}, {}
    pseudobulk = {}
    for st in states:
        mc = datadir / f"expr_metacell_{st}.csv"
        pb = datadir / f"expr_pseudobulk_{st}.txt"
        if not mc.exists():
            print(f"  [跳过] {st}:缺 {mc.name}")
            continue
        print(f"Step 1 · [{st}] 准备数据 ...")
        roots, tgs, fnum, net, Xr, Xt = prepare(
            grn_file, mc, args.edge_num, args.columns, args.min_std)
        print(f"    边 {len(net)} · root TF {len(roots)} · TG {len(tgs)} · metacell {Xr.shape[0]}")

        print(f"Step 2 · [{st}] 拟合高斯 PGM ...")
        mu, sigma, k, losses, gidx = fit_pgm(
            roots, tgs, fnum, net, Xr, Xt, n_iter=args.iters, lr_mu=args.lr_mu,
            lr_sigma=args.lr_sigma, lr_k=args.lr_k, clip_grad=args.clip_grad, seed=args.seed)

        print(f"Step 3 · [{st}] 模拟 TF 敲除,算 deltaX ...")
        dxw = compute_deltax(net, mu, sigma, k, gidx, ko_factor=args.ko_factor)
        dxw.to_csv(outdir / f"deltax_{st}.csv")
        pd.DataFrame({"edge": net["TF"] + "->" + net["TG"], "Score": net["Score"],
                      "Alpha": net["Alpha"], "Cov": net["Cov"], "k": k}
                     ).to_csv(outdir / f"pgm_params_{st}.csv", index=False)

        losses_by_state[st], net_by_state[st], k_by_state[st], dx_by_state[st] = \
            losses, net, k, dxw
        pseudobulk[st] = pd.read_table(pb, header=None, index_col=0,
                                       names=["Gene", "Exp"])["Exp"]

    if not dx_by_state:
        print("没有任何状态跑通,退出。")
        return 1

    # --- Step 4:主控 TF 打分 -------------------------------------------------
    print("Step 4 · 沿轨迹给 TF 打分(余弦相似度)...")
    scores = master_tf_scores(dx_by_state, pseudobulk, trajectories)
    scores = scores.sort_values(["transition", "cosine"]).reset_index(drop=True)
    scores.to_csv(outdir / "master_tf_scores.csv", index=False)
    for tr, g in scores.groupby("transition"):
        top = g.nsmallest(3, "cosine")
        print(f"    {tr}  最负(推动该转换): " +
              ", ".join(f"{r.TF}(cos{r.cosine:+.3f}, |dExpr|={r.obs_shift_norm:.1f})"
                        for r in top.itertuples()))

    # --- Step 5:GRN 迁移段守卫检查 -------------------------------------------
    print("Step 5 · 检查上游 GRN 迁移(transfer_learning)可用性 ...")
    repo = Path(args.cellpolaris_repo) if args.cellpolaris_repo else None
    tl = check_transfer_learning(repo)
    if tl["available"]:
        print("    环境齐备。生成 GRN 请直接跑上游命令:")
    else:
        print("    不可用,缺:")
        for m in tl["missing"]:
            print(f"      - {m}")
        print("    齐备后的上游命令(读自 generate_external_grn.py 的 argparse):")
    print("      python transfer_learning/generate_external_grn.py \\\n"
          "        --species sc_mouse --fold 0 --transfer_loss_type graph_mixup \\\n"
          "        --model ncf --seed 0 --device cuda:0 --top_ratio 0.2 \\\n"
          "        --rna_seq_path <pseudobulk1.txt> <pseudobulk2.txt>")

    # --- Step 6:出图 ---------------------------------------------------------
    print("Step 6 · 出图 ...")
    fig_pgm_fit(losses_by_state, net_by_state, k_by_state, assets / "584_fig1_pgm_fit")
    first = list(dx_by_state)[0]
    fig_deltax_heatmap(dx_by_state[first], first, assets / "584_fig2_deltax_heatmap")
    fig_master_tf(scores, assets / "584_fig3_master_tf")

    summary = {
        "module": "584_cellpolaris_grn_transfer",
        "upstream": {"paper": "Adv Sci 2026;13(12):e08697", "pmid": "41498638",
                     "doi": "10.1002/advs.202508697",
                     "repo": "https://github.com/xCompass-AI/CellPolaris"},
        "seed": args.seed,
        "session": {"python": sys.version.split()[0], "numpy": np.__version__,
                    "pandas": pd.__version__, "torch": torch.__version__},
        "params": vars(args),
        "states": list(dx_by_state),
        "edges_per_state": {s: int(len(net_by_state[s])) for s in net_by_state},
        "final_nll": {s: float(losses_by_state[s][-1]) for s in losses_by_state},
        "transfer_learning_stage": tl,
        "top_master_tf": {tr: g.nsmallest(3, "cosine")[["TF", "cosine", "obs_shift_norm"]]
                          .to_dict("records") for tr, g in scores.groupby("transition")},
    }
    (outdir / "584_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"完成。结果 → {outdir}  图 → {assets}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
