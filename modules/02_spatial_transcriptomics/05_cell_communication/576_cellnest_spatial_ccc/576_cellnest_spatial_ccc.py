"""576 · CellNEST — 空间转录组细胞通讯(图注意力)+ 可跑朴素基线.

CellNEST (Zohora et al., Nature Methods 2025, 22:1505-1519; PMID 40481363,
doi:10.1038/s41592-025-02721-3) 用 GATv2 图注意力在空间转录组上直接推断
配体-受体通讯,并给出注意力导出的置信度与 relay(接力)网络。
仓库 https://github.com/schwartzlab-methods/CellNEST

本模块两条路:
  · 基线(默认,本机 CPU 即跑):空间受限的 LR 共表达乘积打分 —— "无注意力"下限。
    这是 CellNEST 这类模型必须打败的朴素对照,不给对照就不报告模型结果。
  · CellNEST 路径:CellNEST 是 Linux+GPU+singularity 的命令行流水线,本机不可能跑。
    --run-cellnest 只做守卫式检查并打印从官方 vignette 抄录的真实命令;
    --cellnest-csv 可读入真实 CellNEST postprocess 输出,复用同一套出图并与基线比对。

API 来源 —— 全部对照本地克隆的上游源码逐条核实(2026-07-21):
  CLI 子命令表   cellnest (bash 分发脚本) L3-L57
  preprocess 参数 data_preprocess_CellNEST.py L28-L64 (argparse)
  run 参数        run_CellNEST.py L22-L44
  postprocess 参数 output_postprocess_CellNEST.py L42-L55
  visualize 参数  output_visualization_CellNEST.py L88-L120
  relay_extract   extract_relay_cellnest.py L61-L68
  relay_confidence relay_confidence.py L183-L186
  输出 9 列列名   output_postprocess_CellNEST.py L275 / L333
  模型            GATv2Conv_CellNEST.py L19 class GATv2Conv → CCC_gat.py L78-L79
  LR 数据库列名   database/CellNEST_database.csv 表头 = Ligand,Receptor,Annotation,Reference
  上游许可证      LICENSE = GNU GPL v3

本模块**不 import 任何 CellNEST 代码**(上游无 setup.py/pyproject,不是可安装的 Python 包),
只复用其 CLI 命令字符串与输出 CSV schema,两者均已在上述源码位置逐字核对。
"""
from __future__ import annotations

import argparse
import os
import shutil
import sys
import warnings

warnings.filterwarnings("ignore")

# Windows 控制台默认 GBK,中文/✓✗ 会 UnicodeEncodeError;强制 UTF-8 输出
for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

HERE = os.path.dirname(os.path.abspath(__file__))
EXAMPLE = os.path.join(HERE, "example_data")
RESULTS = os.path.join(HERE, "results")
ASSETS = os.path.join(HERE, "assets")
FRAMEWORK = os.path.abspath(os.path.join(HERE, "..", "..", "..", "_framework"))
sys.path.insert(0, FRAMEWORK)

# CellNEST postprocess 写出的真实列顺序(headerless CSV)。
# 来源:output_postprocess_CellNEST.py L275 / L333 —— csv_record.append([...])
# 注意:该文件里的 'component' 列被上游硬编码为 -1(L284 / L342 `label = -1`),
#       真正的连通分量是在**可视化步骤**里算的(output_visualization_CellNEST.py L291,
#       connected_components(directed=True, connection='weak'))。读入真实输出时别指望
#       component 有意义 —— 本模块基线自己算 component,做法同上游可视化步骤但**不是**
#       postprocess CSV 里那一列。
CELLNEST_COLUMNS = ["from_cell", "to_cell", "ligand", "receptor", "edge_rank",
                    "component", "from_id", "to_id", "attention_score"]

# 官方 vignette 里的真实命令(逐字抄录自 vignette/workflow.md L24/L47/L64/L75/L160/L184
# 与 README.md L9 的 singularity 行,占位符用 <> 标出),守卫路径打印给用户。
# 每个参数名都已在对应脚本的 argparse 里核对过(见模块顶部 docstring 的行号表)。
CELLNEST_CLI = """\
# 0) 环境(Linux + NVIDIA GPU;本机 Windows 无法运行)
singularity pull cellnest_image.sif library://fatema/collection/cellnest_image.sif:latest
#    或克隆仓库后:  sudo bash setup.sh

# ★ 上游已知缺陷(2026-07-21 对照源码,原样报告,不要被绊住):
#   a) setup.sh 里写的是  chmod +x nest; cp nest $HOME/.local/bin/  —— 但仓库里的脚本
#      叫 `cellnest`,没有 `nest` 文件,直接跑会失败。手动改成 cellnest 或自己 cp。
#   b) cellnest 分发脚本 L44 把 relay_extract 指向 extract_relay_nest.py,
#      而仓库里的文件是 extract_relay_cellnest.py —— relay_extract 子命令会静默无输出。
#      绕过办法:直接 python -u extract_relay_cellnest.py <同样的参数>

# 1) 预处理(建图 + Quantile 归一化;输入必须是原始 counts)
cellnest preprocess --data_name='<SAMPLE>' --data_from='data/<SAMPLE>/' --data_type=visium
#    单细胞分辨率数据加 --distance_measure='knn';忽略自分泌加 --block_autocrine=1
#    已归一化的数据必须加 --skip_normalize=1
#    归一化用的是 qnorm.quantile_normalize(分位数归一化,data_preprocess_CellNEST.py L113),
#      不是 sklearn.QuantileTransformer —— 本模块基线用后者,是近似不是复刻。
#    邻域: vignette/workflow.md 写「默认 --neighborhood_threshold = spot_diameter*4
#      (--spot_diameter=89.43)」,但源码 L246-251 实际是「第一个 spot 到最近邻的距离 ×4」,
#      且 argparse 里**根本没有 --spot_diameter 这个参数**(文档与实现不一致,以源码为准)。
#    --filter_min_cell: vignette 说默认 5,argparse L34 实际 default=1(同上,文档≠实现)。
#    --threshold_gene_exp 默认 98 (L35) —— 文档与实现一致。

# 2) 训练(推荐 5 个种子后做 ensemble;单次约 13 小时 @ V100)
nohup cellnest run --data_name='<SAMPLE>' --num_epoch 80000 --manual_seed='yes' \\
      --seed=1 --model_name='CellNEST_<SAMPLE>' --run_id=1 > run1.log &

# 3) 后处理(按 rank 乘积做 ensemble,取 top 20%)
cellnest postprocess --data_name='<SAMPLE>' --model_name='CellNEST_<SAMPLE>' --total_runs=5

# 4) 可视化 / relay 网络 / 置信度
cellnest visualize --data_name='<SAMPLE>' --model_name='CellNEST_<SAMPLE>' --top_edge_count=40000
cellnest relay_extract --data_name='<SAMPLE>' --metadata='metadata/' \\
      --top_ccc_file='output/<SAMPLE>/CellNEST_<SAMPLE>_ccc_list_top3000.csv' \\
      --output_path='output/<SAMPLE>/'
cellnest relay_confidence --input_path='output/<SAMPLE>/CellNEST_<SAMPLE>_relay_pattern_count.csv' \\
      --output_path='output/<SAMPLE>/relay_confidence_score_for_top3kCCC.csv' \\
      --organism='human' --database_dir='database/'

# 产物: output/<SAMPLE>/CellNEST_<SAMPLE>_top20percent.csv  (9 列)
#        from_cell,to_cell,ligand,receptor,edge_rank,component,from_id,to_id,attention_score
#   注意:上游用 to_csv(..., header=False) 写出,但被写的 DataFrame 第 0 行本身就是列名
#   (output_postprocess_CellNEST.py L292/L396 把 csv_record[0]=表头行一起放进 df),
#   所以文件第一行实际上是这串列名。本模块读入时自动判别有/无表头两种情况。
# 把该文件喂回本模块:  python 576_cellnest_spatial_ccc.py --cellnest-csv <该文件>
"""


# ---------------------------------------------------------------- 合成示例数据
def make_example(seed: int = 0) -> None:
    """生成小型合成空间转录组(synthetic, for demo only)。

    设计:两块空间区域(region A / region B),各自富集一组配体-受体对,
    使得"真信号"是空间局部的 —— 基线若有效,应把边集中在区域内部。
    """
    import numpy as np
    import pandas as pd

    rng = np.random.default_rng(seed)
    os.makedirs(EXAMPLE, exist_ok=True)

    # --- 蜂窝状 spot 网格(近似 Visium 布局),两个区域 ---
    coords, region = [], []
    n_row, n_col = 14, 14
    pitch = 100.0                      # spot 中心间距(与 spot_diameter=89.43 同量级)
    for r in range(n_row):
        for c in range(n_col):
            x = c * pitch + (pitch / 2 if r % 2 else 0)
            y = r * pitch * 0.866
            coords.append((x, y))
            region.append("A" if (x + y) < (n_col * pitch * 0.75) else "B")
    coords = np.asarray(coords)
    region = np.asarray(region)
    n_spot = len(coords)
    barcodes = [f"SPOT_{i:04d}" for i in range(n_spot)]

    # --- 基因集合:两组 LR + 背景基因 ---
    lr_A = [("TGFB1", "TGFBR1"), ("TGFB1", "TGFBR2"), ("CCL19", "CCR7"), ("CXCL12", "CXCR4")]
    lr_B = [("VEGFA", "KDR"), ("PDGFB", "PDGFRB"), ("IL6", "IL6R"), ("WNT5A", "FZD4")]
    lr_pairs = lr_A + lr_B
    sig_genes = sorted({g for pair in lr_pairs for g in pair})
    bg_genes = [f"BG{i:03d}" for i in range(60)]
    genes = sig_genes + bg_genes

    # 背景:负二项样的计数
    counts = rng.negative_binomial(4, 0.35, size=(n_spot, len(genes))).astype(float)

    # 区域特异地把对应 LR 基因抬高
    gidx = {g: k for k, g in enumerate(genes)}
    for pairs, reg in ((lr_A, "A"), (lr_B, "B")):
        mask = region == reg
        for lig, rec in pairs:
            for g in (lig, rec):
                boost = rng.negative_binomial(18, 0.35, size=mask.sum())
                counts[mask, gidx[g]] += boost

    pd.DataFrame(counts.astype(int), index=barcodes, columns=genes).to_csv(
        os.path.join(EXAMPLE, "spatial_counts.csv"), index_label="barcode")
    pd.DataFrame({"barcode": barcodes, "x": coords[:, 0], "y": coords[:, 1],
                  "region": region}).to_csv(
        os.path.join(EXAMPLE, "spatial_coordinates.csv"), index=False)
    # 列名与 CellNEST 官方数据库一致:Ligand,Receptor,Annotation,Reference
    pd.DataFrame({
        "Ligand":    [p[0] for p in lr_pairs],
        "Receptor":  [p[1] for p in lr_pairs],
        "Annotation": ["Secreted Signaling"] * len(lr_pairs),
        "Reference": ["synthetic, for demo only"] * len(lr_pairs),
    }).to_csv(os.path.join(EXAMPLE, "lr_pairs.csv"), index=False)
    print(f"[example_data] 合成数据已写出: {n_spot} spots × {len(genes)} genes, "
          f"{len(lr_pairs)} LR pairs (synthetic, for demo only)")


# ---------------------------------------------------------------- 基线实现
def baseline_ccc(counts, coords_df, lr_df, *, neighborhood, thr_gene_exp,
                 filter_min_cell, top_percent, block_autocrine, seed):
    """朴素基线:空间受限的配体-受体共表达乘积(无注意力、无图神经网络)。

    对应 CellNEST 论文里被比较的一类"共表达乘积"方法,是本模块的可跑下限。
    预处理规则参照 CellNEST vignette 的描述(分位数归一化 + 每个 spot 取高分位活跃基因
    + spot_diameter*4 邻域),但**打分本身是共表达乘积,不是 CellNEST 的注意力**。
    三处刻意的近似(不要当成复刻 CellNEST):
      · 归一化用 sklearn.QuantileTransformer,上游用 qnorm.quantile_normalize;
      · 邻域半径用 spot_diameter*4(vignette 的说法),上游源码实际用最近邻距离*4;
      · component 由本函数自算,上游 postprocess CSV 里那一列是 -1。
    """
    import numpy as np
    import pandas as pd
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    from scipy.spatial import cKDTree
    from sklearn.preprocessing import QuantileTransformer

    rng = np.random.default_rng(seed)

    # 1) 基因过滤:至少在 filter_min_cell 个 spot 中表达(CellNEST --filter_min_cell)
    keep = (counts.values > 0).sum(axis=0) >= filter_min_cell
    X = counts.loc[:, keep]
    print(f"  [1] 基因过滤 --filter_min_cell={filter_min_cell}: "
          f"{counts.shape[1]} -> {X.shape[1]} genes")

    # 2) 分位数归一化(CellNEST 默认对表达矩阵做 Quantile transform)
    qt = QuantileTransformer(output_distribution="uniform", random_state=seed,
                             n_quantiles=min(1000, X.shape[0]))
    Xn = pd.DataFrame(qt.fit_transform(X.values), index=X.index, columns=X.columns)

    # 3) 每个 spot 内取高分位基因为"活跃"(CellNEST --threshold_gene_exp)
    cutoff = np.percentile(Xn.values, thr_gene_exp, axis=1, keepdims=True)
    active = Xn.values >= cutoff
    print(f"  [2] 分位数归一化 + 活跃基因阈值 --threshold_gene_exp={thr_gene_exp}: "
          f"平均每 spot {active.sum(1).mean():.1f} 个活跃基因")

    # 4) 空间邻域图(CellNEST 默认 --neighborhood_threshold = spot_diameter*4)
    XY = coords_df[["x", "y"]].values
    tree = cKDTree(XY)
    pairs = tree.query_pairs(r=neighborhood, output_type="ndarray")
    # 通讯有方向:i->j 与 j->i 都要
    directed = np.vstack([pairs, pairs[:, ::-1]])
    if not block_autocrine:                       # 自分泌 = 自环
        self_loops = np.repeat(np.arange(len(XY))[:, None], 2, axis=1)
        directed = np.vstack([directed, self_loops])
    print(f"  [3] 邻域图 r={neighborhood:.1f}: {len(directed)} 条有向 spot-spot 连接"
          f"{' (含自分泌)' if not block_autocrine else ' (--block_autocrine)'}")

    # 5) LR 打分:score = normed(ligand @ sender) * normed(receptor @ receiver)
    col = {g: k for k, g in enumerate(X.columns)}
    lr_use = [(l, r) for l, r in zip(lr_df["Ligand"], lr_df["Receptor"])
              if l in col and r in col]
    if not lr_use:
        raise SystemExit("配体-受体对与表达矩阵无交集,检查基因命名(需同一 symbol 体系)")

    src, dst = directed[:, 0], directed[:, 1]
    recs = []
    for lig, rec in lr_use:
        li, ri = col[lig], col[rec]
        ok = active[src, li] & active[dst, ri]     # 两端都必须活跃
        if not ok.any():
            continue
        s = Xn.values[src[ok], li] * Xn.values[dst[ok], ri]
        recs.append(pd.DataFrame({"from_id": src[ok], "to_id": dst[ok],
                                  "ligand": lig, "receptor": rec, "score": s}))
    if not recs:
        raise SystemExit(
            "没有任何 LR 对在两端同时通过活跃基因阈值 —— 请调低 --threshold-gene-exp "
            "或调大 --neighborhood 后重试")
    edges = pd.concat(recs, ignore_index=True)
    print(f"  [4] LR 打分: {len(lr_use)} 个 LR 对 -> {len(edges)} 条候选通讯边")

    # 6) 取 top_percent(CellNEST postprocess --top_percent=20)
    edges = edges.sort_values("score", ascending=False).reset_index(drop=True)
    n_top = max(1, int(round(len(edges) * top_percent / 100)))
    top = edges.iloc[:n_top].copy()
    top["edge_rank"] = np.arange(1, len(top) + 1)
    print(f"  [5] 取 top {top_percent}%: {len(top)} 条边 "
          f"(score {top['score'].min():.3f} – {top['score'].max():.3f})")

    # 7) 连通分量标签。做法对齐上游**可视化步骤**(output_visualization_CellNEST.py L291
    #    用 scipy connected_components 给 spot 分配 component);上游 postprocess CSV 的
    #    component 列本身是硬编码 -1(L284/L342),不要拿来比对。
    n = len(XY)
    adj = coo_matrix((np.ones(len(top)), (top["from_id"], top["to_id"])), shape=(n, n))
    _, labels = connected_components(adj, directed=False)
    top["component"] = labels[top["from_id"].values]

    bc = coords_df["barcode"].values
    out = pd.DataFrame({
        "from_cell": bc[top["from_id"].values],
        "to_cell":   bc[top["to_id"].values],
        "ligand":    top["ligand"].values,
        "receptor":  top["receptor"].values,
        "edge_rank": top["edge_rank"].values,
        "component": top["component"].values,
        "from_id":   top["from_id"].values,
        "to_id":     top["to_id"].values,
        "attention_score": top["score"].values,   # 基线里这是共表达乘积,不是注意力
    })[CELLNEST_COLUMNS]

    # 8) 空间置换零模型:打乱 spot 的空间位置(表达谱整体搬家),邻域图结构不变,
    #    只破坏"表达与位置的对应"。注意:分位数归一化让分数的边际分布几乎不受置换影响,
    #    真正的空间信号体现在**通过双端活跃门的候选边数量**上,所以按 LR 对比较边数。
    n_perm = 20
    obs_cnt = edges.groupby(["ligand", "receptor"]).size()
    null_cnt = {k: [] for k in obs_cnt.index}
    for _ in range(n_perm):
        perm = rng.permutation(len(XY))
        for lig, rec in lr_use:
            li, ri = col[lig], col[rec]
            k = (lig, rec)
            if k in null_cnt:
                null_cnt[k].append(int((active[perm[src], li] & active[perm[dst], ri]).sum()))
    ctrl = pd.DataFrame({
        "ligand":   [k[0] for k in obs_cnt.index],
        "receptor": [k[1] for k in obs_cnt.index],
        "n_observed": obs_cnt.values,
        "n_permuted_mean": [float(np.mean(null_cnt[k])) for k in obs_cnt.index],
        "n_permuted_sd":   [float(np.std(null_cnt[k])) for k in obs_cnt.index],
    })
    ctrl["enrichment"] = ctrl["n_observed"] / ctrl["n_permuted_mean"].replace(0, np.nan)
    # 经验 p:置换中达到或超过观测边数的比例(加一平滑)
    ctrl["p_perm"] = [ (1 + sum(np.asarray(null_cnt[k]) >= obs_cnt[k])) / (n_perm + 1)
                       for k in obs_cnt.index ]
    ctrl = ctrl.sort_values("enrichment", ascending=False).reset_index(drop=True)
    print(f"  [6] 空间置换零模型 ({n_perm} 次): 候选边 {len(edges)} (观测) vs "
          f"{ctrl['n_permuted_mean'].sum():.0f} (置换均值), "
          f"中位富集 {ctrl['enrichment'].median():.2f}x")

    return out, ctrl, coords_df


# ---------------------------------------------------------------- 出图
def make_figures(ccc, coords_df, ctrl, outdir, label):
    """4 张顶刊风图。规矩:不用条形图 —— 用散点/线段/lollipop/heatmap/violin。"""
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    from matplotlib.collections import LineCollection
    from pubstyle import CMAP_CONT, NATURE_W1, pal, save_fig, set_pub_style

    set_pub_style(base_size=10)
    os.makedirs(outdir, exist_ok=True)
    XY = coords_df[["x", "y"]].values
    made = []

    # --- Fig 1 空间通讯图:spot 散点 + top 边线段(viridis 按分数) ---
    fig, ax = plt.subplots(figsize=(NATURE_W1 * 1.6, NATURE_W1 * 1.5))
    ax.scatter(XY[:, 0], XY[:, 1], s=16, c="#DCDCDC", edgecolors="none", zorder=1)
    sub = ccc.nlargest(min(1200, len(ccc)), "attention_score")
    seg = np.stack([XY[sub["from_id"].values], XY[sub["to_id"].values]], axis=1)
    lc = LineCollection(seg, array=sub["attention_score"].values, cmap=CMAP_CONT,
                        linewidths=0.7, alpha=0.85, zorder=2)
    ax.add_collection(lc)
    cb = fig.colorbar(lc, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label("Communication score", fontsize=9)
    ax.set(xlabel="Spatial x", ylabel="Spatial y",
           title=f"Spatially resolved cell-cell communication ({label})")
    ax.set_aspect("equal")
    save_fig(fig, os.path.join(outdir, "fig1_spatial_ccc_map"))
    plt.close(fig)
    made.append("fig1_spatial_ccc_map.png")

    # --- Fig 2 lollipop:各 LR 对保留的边数(明确不用条形图) ---
    cnt = (ccc.assign(lr=ccc["ligand"] + " → " + ccc["receptor"])
              .groupby("lr").agg(n=("edge_rank", "size"),
                                 mean_score=("attention_score", "mean"))
              .sort_values("n").tail(15))
    fig, ax = plt.subplots(figsize=(NATURE_W1 * 1.6, 0.30 * len(cnt) + 1.4))
    y = np.arange(len(cnt))
    ax.hlines(y, 0, cnt["n"], color="#C8C8C8", linewidth=1.4, zorder=1)
    sc = ax.scatter(cnt["n"], y, c=cnt["mean_score"], cmap=CMAP_CONT,
                    s=70, zorder=2, edgecolors="white", linewidths=0.6)
    cb = fig.colorbar(sc, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label("Mean score", fontsize=9)
    ax.set_yticks(y, cnt.index, fontsize=8)
    ax.set(xlabel="Retained communication edges", ylabel="",
           title=f"Top ligand–receptor pairs ({label})")
    ax.set_xlim(left=0)
    save_fig(fig, os.path.join(outdir, "fig2_lr_lollipop"))
    plt.close(fig)
    made.append("fig2_lr_lollipop.png")

    # --- Fig 3 heatmap:LR 对 × 空间连通分量的平均分数 ---
    keep_comp = ccc["component"].value_counts().head(10).index
    m = (ccc[ccc["component"].isin(keep_comp)]
         .assign(lr=ccc["ligand"] + " → " + ccc["receptor"])
         .pivot_table(index="lr", columns="component",
                      values="attention_score", aggfunc="mean"))
    m = m.loc[m.mean(axis=1).sort_values(ascending=False).index]
    fig, ax = plt.subplots(figsize=(NATURE_W1 * 1.5, 0.28 * len(m) + 1.6))
    im = ax.imshow(m.values, cmap=CMAP_CONT, aspect="auto")
    ax.set_xticks(range(m.shape[1]), [str(c) for c in m.columns], fontsize=8)
    ax.set_yticks(range(m.shape[0]), m.index, fontsize=8)
    cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
    cb.set_label("Mean score", fontsize=9)
    ax.set(xlabel="Spatial component", ylabel="",
           title=f"LR activity across spatial components ({label})")
    save_fig(fig, os.path.join(outdir, "fig3_lr_component_heatmap"))
    plt.close(fig)
    made.append("fig3_lr_component_heatmap.png")

    # --- Fig 4 dumbbell:每个 LR 对 观测 vs 空间置换的候选边数(空间特异性对照) ---
    #     分数的边际分布对置换几乎不敏感(分位数归一化所致),真正的空间信号在边数上,
    #     所以这里比边数而不是比分数。
    if ctrl is not None and len(ctrl):
        d = ctrl.assign(lr=ctrl["ligand"] + " → " + ctrl["receptor"]) \
                .sort_values("enrichment").tail(15)
        y = np.arange(len(d))
        c_obs, c_null = pal(2)
        fig, ax = plt.subplots(figsize=(NATURE_W1 * 1.7, 0.30 * len(d) + 1.6))
        ax.hlines(y, d["n_permuted_mean"], d["n_observed"],
                  color="#C8C8C8", linewidth=1.6, zorder=1)
        ax.errorbar(d["n_permuted_mean"], y, xerr=d["n_permuted_sd"], fmt="none",
                    ecolor="#9E9E9E", elinewidth=1.0, capsize=2, zorder=2)
        ax.scatter(d["n_permuted_mean"], y, s=62, color=c_null, zorder=3,
                   edgecolors="white", linewidths=0.6, label="Spatial permutation")
        ax.scatter(d["n_observed"], y, s=62, color=c_obs, zorder=3,
                   edgecolors="white", linewidths=0.6, label="Observed")
        ax.set_yticks(y, d["lr"], fontsize=8)
        ax.set(xlabel="Candidate communication edges", ylabel="",
               title=f"Spatial specificity control ({label})")
        ax.set_xlim(left=0)
        ax.legend(loc="lower left", fontsize=8)   # 右侧被观测点占满,图例放左下
        save_fig(fig, os.path.join(outdir, "fig4_permutation_control"))
        plt.close(fig)
        made.append("fig4_permutation_control.png")

    return made


# ---------------------------------------------------------------- CellNEST 守卫
def cellnest_guard(repo: str | None) -> None:
    """守卫式引用封装:CellNEST 是 Linux+GPU CLI 流水线,本机不运行,只做检查+指路。"""
    print("\n=== CellNEST 真实流水线检查 ===")
    exe = shutil.which("cellnest")
    ok = True
    if exe is None:
        print("  ✗ 未找到 `cellnest` 命令(需在 Linux 上 `sudo bash setup.sh` 或用 singularity 镜像)")
        ok = False
    else:
        print(f"  ✓ 找到 cellnest: {exe}")
        # 注意:上游 `cellnest` 只是个 if/elif 分发脚本(cellnest L3-L57),**没有 --help 分支**,
        # 传 --help 会静默落空、退出码 0 —— 所以这里不调用它,也不拿它当可用性证据。
        # 各子命令的帮助要向下游脚本要,例如 `python -u data_preprocess_CellNEST.py --help`。
        print("    (上游 cellnest 无 --help 分支,子命令帮助见 "
              "`python -u data_preprocess_CellNEST.py --help`)")
    if repo and not os.path.isdir(repo):
        print(f"  ✗ --cellnest-repo 路径不存在: {repo}")
        ok = False
    try:
        import torch
        print(f"  {'✓' if torch.cuda.is_available() else '✗'} torch CUDA available = "
              f"{torch.cuda.is_available()}")
        if not torch.cuda.is_available():
            ok = False
    except Exception:
        print("  ✗ 未安装 torch")
        ok = False
    if sys.platform.startswith("win"):
        print("  ✗ 当前是 Windows;CellNEST 官方仅在 Linux (CentOS 7) + NVIDIA GPU 上测试")
        ok = False

    if not ok:
        print("\n本机无法运行 CellNEST。以下是官方 vignette 的真实命令,请在 Linux GPU 机上执行:\n")
        print(CELLNEST_CLI)
        print("上述参数名均已对照上游 argparse 逐条核实(见模块顶部行号表);"
              "跨版本使用请以你那份源码的 argparse 为准。")
    print("=== 基线结果不受影响,已照常产出 ===\n")


def read_cellnest_csv(path: str, coords_df):
    """读入真实 CellNEST postprocess 输出(*_top20percent.csv,9 列)。

    上游写盘用 header=False,但表头行被当成数据行写在第一行,所以文件通常仍带列名;
    这里对两种情况都做判别(见 CELLNEST_CLI 里的溯源注释)。
    """
    import pandas as pd
    head = pd.read_csv(path, nrows=1, header=None)
    has_header = str(head.iloc[0, 0]).strip() == "from_cell"
    df = pd.read_csv(path, header=0 if has_header else None)
    if df.shape[1] != len(CELLNEST_COLUMNS):
        raise SystemExit(
            f"{path} 有 {df.shape[1]} 列,期望 {len(CELLNEST_COLUMNS)} 列: {CELLNEST_COLUMNS}")
    df.columns = CELLNEST_COLUMNS
    n = len(coords_df)
    if df["from_id"].max() >= n or df["to_id"].max() >= n:
        raise SystemExit("CellNEST CSV 的 from_id/to_id 超出坐标文件行数,坐标与输出不匹配")
    print(f"[cellnest] 读入 {len(df)} 条真实 CellNEST 通讯边: {path}")
    return df


# ---------------------------------------------------------------- 依赖快照
def save_session(outdir: str, args) -> None:
    """落盘依赖版本 + 本次参数(铁律6:可复现)。"""
    import importlib
    import platform
    lines = [f"python: {sys.version.split()[0]}  ({platform.platform()})",
             f"seed:   {args.seed}", "", "packages:"]
    for m in ("numpy", "pandas", "scipy", "sklearn", "matplotlib"):
        try:
            lines.append(f"  {m}: {importlib.import_module(m).__version__}")
        except Exception as e:
            lines.append(f"  {m}: <unavailable: {type(e).__name__}>")
    lines += ["", "args:"] + [f"  {k}: {v}" for k, v in sorted(vars(args).items())]
    with open(os.path.join(outdir, "session_info.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


# ---------------------------------------------------------------- main
def main() -> None:
    p = argparse.ArgumentParser(description="576 · CellNEST 空间通讯 + 可跑基线")
    p.add_argument("--counts", default=os.path.join(EXAMPLE, "spatial_counts.csv"))
    p.add_argument("--coords", default=os.path.join(EXAMPLE, "spatial_coordinates.csv"))
    p.add_argument("--lr-db", default=os.path.join(EXAMPLE, "lr_pairs.csv"))
    p.add_argument("--outdir", default=RESULTS)
    p.add_argument("--spot-diameter", type=float, default=89.43,
                   help="Visium spot 直径;89.43 取自 CellNEST vignette 的说明文字"
                        "(注:上游 argparse 里并没有 --spot_diameter 这个参数)")
    p.add_argument("--neighborhood", type=float, default=None,
                   help="邻域半径,默认 spot_diameter*4(vignette 的说法;上游源码实际是"
                        "最近邻距离*4)")
    p.add_argument("--threshold-gene-exp", type=float, default=80.0,
                   help="活跃基因分位阈;CellNEST 默认 98(源码已核实),合成小数据用 80")
    p.add_argument("--filter-min-cell", type=int, default=5,
                   help="基因至少在几个 spot 中表达;CellNEST vignette 说默认 5,"
                        "但其 argparse 实际 default=1")
    p.add_argument("--top-percent", type=float, default=20.0)
    p.add_argument("--block-autocrine", action="store_true", help="忽略自分泌(自环)")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--run-cellnest", action="store_true",
                   help="检查真实 CellNEST 流水线可用性并打印官方命令")
    p.add_argument("--cellnest-repo", default=None)
    p.add_argument("--cellnest-csv", default=None,
                   help="真实 CellNEST *_top20percent.csv,读入后复用同一套出图")
    p.add_argument("--regen-example", action="store_true", help="重新生成合成示例数据")
    p.add_argument("--save-assets", action="store_true", help="同时把图写进 assets/")
    args = p.parse_args()

    import numpy as np
    import pandas as pd

    np.random.seed(args.seed)
    os.makedirs(args.outdir, exist_ok=True)

    if args.regen_example or not os.path.exists(args.counts):
        make_example(seed=0)

    print("=== 576 · CellNEST 空间细胞通讯 ===")
    counts = pd.read_csv(args.counts, index_col=0)
    coords_df = pd.read_csv(args.coords)
    lr_df = pd.read_csv(args.lr_db)
    print(f"输入: {counts.shape[0]} spots × {counts.shape[1]} genes · "
          f"{len(lr_df)} LR pairs")

    nb = args.neighborhood if args.neighborhood else args.spot_diameter * 4
    print("\n[基线] 空间受限 LR 共表达乘积(无注意力下限)")
    ccc, ctrl, coords_df = baseline_ccc(
        counts, coords_df, lr_df, neighborhood=nb,
        thr_gene_exp=args.threshold_gene_exp, filter_min_cell=args.filter_min_cell,
        top_percent=args.top_percent, block_autocrine=args.block_autocrine,
        seed=args.seed)

    # 与 CellNEST 一致的 headerless 输出 + 带表头的易读版
    ccc.to_csv(os.path.join(args.outdir, "baseline_top_ccc_cellnest_schema.csv"),
               index=False, header=False)
    ccc.to_csv(os.path.join(args.outdir, "baseline_top_ccc.csv"), index=False)

    summary = (ccc.assign(lr=ccc["ligand"] + "→" + ccc["receptor"])
                  .groupby("lr").agg(n_edges=("edge_rank", "size"),
                                     mean_score=("attention_score", "mean"),
                                     best_rank=("edge_rank", "min"))
                  .sort_values("n_edges", ascending=False))
    summary.to_csv(os.path.join(args.outdir, "baseline_lr_summary.csv"))
    ctrl.to_csv(os.path.join(args.outdir, "baseline_permutation_control.csv"), index=False)
    print(f"\n[输出] {args.outdir}/baseline_top_ccc.csv "
          f"(+ CellNEST 9 列 schema 版) · baseline_lr_summary.csv")

    print("\n[出图]")
    figs = make_figures(ccc, coords_df, ctrl, args.outdir, "baseline")
    for f in figs:
        print(f"  {args.outdir}/{f}")

    # 真实 CellNEST 输出:同一套图 + 与基线的秩一致性
    if args.cellnest_csv:
        real = read_cellnest_csv(args.cellnest_csv, coords_df)
        make_figures(real, coords_df, None, args.outdir, "CellNEST")
        key = ["from_id", "to_id", "ligand", "receptor"]
        merged = ccc.merge(real, on=key, suffixes=("_base", "_nest"))
        if len(merged) > 2:
            from scipy.stats import spearmanr
            rho, pv = spearmanr(merged["edge_rank_base"], merged["edge_rank_nest"])
            print(f"[比对] 基线 vs CellNEST 共享 {len(merged)} 条边, "
                  f"edge_rank Spearman rho={rho:.3f} (p={pv:.2g})")
            pd.DataFrame({"n_shared": [len(merged)], "spearman_rho": [rho],
                          "p_value": [pv]}).to_csv(
                os.path.join(args.outdir, "baseline_vs_cellnest_concordance.csv"),
                index=False)
        else:
            print("[比对] 共享边过少,跳过秩一致性")

    if args.run_cellnest:
        cellnest_guard(args.cellnest_repo)

    save_session(args.outdir, args)

    if args.save_assets:
        import shutil as sh
        os.makedirs(ASSETS, exist_ok=True)
        for f in figs:
            sh.copy(os.path.join(args.outdir, f), os.path.join(ASSETS, f))
        print(f"[assets] 展示图已复制到 {ASSETS}")

    print("\n完成。")


if __name__ == "__main__":
    main()
