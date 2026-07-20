"""生成 579 模块的合成示例数据(synthetic, for demo only)。

写出 example_data/ 下 5 个 CSV:一张空间转录组切片(spots)+ 两个单细胞模态
(RNA、第二模态如甲基化 gene-activity),都带 ground-truth 层标签,
用来对"细胞→spot 的空间分配"做可量化打分。

只需跑一次;主脚本 579_simo_spatial_multiomics.py 会直接读这些 CSV。
"""
from __future__ import annotations
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "example_data")
SEED = 2026

LAYERS = ["L2_3", "L4", "L5", "L6"]      # 4 个皮层样"层",沿 y 轴分带
N_GENES = 60
N_SPOTS_PER_ROW = 14
N_ROWS = 16
N_RNA_CELLS = 320
N_MOD2_CELLS = 180


def _programs(rng):
    """每层一套基因程序:层特异 marker + 相邻层共享梯度 + 共享背景。

    相邻层程序刻意重叠(真实皮层就是连续梯度而非离散块),
    否则任务过于简单——单纯按相关性贪心就能 100% 命中,基线对比失去意义。
    """
    prog = rng.gamma(1.2, 1.0, size=(len(LAYERS), N_GENES)) * 0.6
    block = N_GENES // len(LAYERS)
    for i in range(len(LAYERS)):                       # 每层的 marker 段
        prog[i, i * block:(i + 1) * block] += rng.uniform(1.2, 2.0, block)
    # 相邻层互相渗透:marker 不是层独占的,而是逐层衰减
    smooth = prog.copy()
    for i in range(len(LAYERS)):
        for j in range(len(LAYERS)):
            if i != j:
                smooth[i] += prog[j] * 0.35 ** abs(i - j)
    return smooth


def main():
    os.makedirs(OUT, exist_ok=True)
    rng = np.random.default_rng(SEED)
    prog = _programs(rng)
    genes = [f"Gene{i:03d}" for i in range(N_GENES)]

    # ---- 空间转录组 spots:规则网格,y 决定层归属 ----
    xs, ys, spot_layer = [], [], []
    for r in range(N_ROWS):
        for c in range(N_SPOTS_PER_ROW):
            xs.append(c + 0.5 * (r % 2))               # 交错网格,像 Visium
            ys.append(r)
            spot_layer.append(LAYERS[min(r * len(LAYERS) // N_ROWS, len(LAYERS) - 1)])
    n_spots = len(xs)
    spot_ids = [f"spot{i:04d}" for i in range(n_spots)]
    li = np.array([LAYERS.index(l) for l in spot_layer])
    # spot 是多细胞混合:主导层 + 邻层污染,更接近真实 Visium
    mix = np.zeros((n_spots, len(LAYERS)))
    mix[np.arange(n_spots), li] = 0.7
    for k in range(n_spots):
        for nb in (li[k] - 1, li[k] + 1):
            if 0 <= nb < len(LAYERS):
                mix[k, nb] += 0.15
    mix /= mix.sum(1, keepdims=True)
    st_expr = rng.poisson(np.clip(mix @ prog, 0.05, None) * 6.0).astype(float)

    pd.DataFrame(st_expr, index=spot_ids, columns=genes).to_csv(
        os.path.join(OUT, "st_expression.csv"))
    pd.DataFrame({"spot": spot_ids, "x": xs, "y": ys, "layer": spot_layer}).to_csv(
        os.path.join(OUT, "st_coords.csv"), index=False)

    # ---- 模态 1:scRNA-seq,与 ST 共享基因,但存在平台效应与 dropout ----
    cell_layer = rng.choice(LAYERS, N_RNA_CELLS)
    ci = np.array([LAYERS.index(l) for l in cell_layer])
    depth = rng.uniform(0.5, 1.8, N_RNA_CELLS)[:, None]     # 细胞间测序深度差异
    # 平台效应:同一基因在 ST 与 scRNA 上捕获效率不同(不同建库化学),
    # 这正是"跨模态对齐"要解决的问题,不加就等于假设两平台可直接比
    platform = rng.lognormal(0, 0.35, N_GENES)[None, :]
    lam = prog[ci] * depth * platform * 4.0
    rna = rng.poisson(np.clip(lam, 0.02, None)).astype(float)
    rna *= rng.random(rna.shape) > 0.25                     # dropout:单细胞稀疏性
    rna_ids = [f"rna{i:04d}" for i in range(N_RNA_CELLS)]
    pd.DataFrame(rna, index=rna_ids, columns=genes).to_csv(
        os.path.join(OUT, "sc_rna_expression.csv"))
    pd.DataFrame({"cell": rna_ids, "layer": cell_layer}).to_csv(
        os.path.join(OUT, "sc_rna_meta.csv"), index=False)

    # ---- 模态 2:非转录组(如甲基化 gene-activity),与表达负相关 ----
    m2_layer = rng.choice(LAYERS, N_MOD2_CELLS)
    mi = np.array([LAYERS.index(l) for l in m2_layer])
    # 甲基化在高表达基因上偏低 -> 取负号再压到 [0,1],模拟 mCH/mCG gene body 信号
    raw = -prog[mi] + rng.normal(0, 1.8, (N_MOD2_CELLS, N_GENES))
    m2 = (raw - raw.min()) / (raw.max() - raw.min())
    m2_ids = [f"met{i:04d}" for i in range(N_MOD2_CELLS)]
    pd.DataFrame(m2, index=m2_ids, columns=genes).to_csv(
        os.path.join(OUT, "sc_mod2_gene_activity.csv"))
    pd.DataFrame({"cell": m2_ids, "layer": m2_layer}).to_csv(
        os.path.join(OUT, "sc_mod2_meta.csv"), index=False)

    print(f"[579] synthetic data -> {OUT}")
    print(f"       spots={n_spots}  rna_cells={N_RNA_CELLS}  mod2_cells={N_MOD2_CELLS}  genes={N_GENES}")


if __name__ == "__main__":
    main()
