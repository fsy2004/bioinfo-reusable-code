# -*- coding: utf-8 -*-
"""生成 590 模块的合成示例数据(synthetic, for demo only)。

只需跑一次,产物已提交到 example_data/。重跑会得到同样的文件(固定种子)。
"""
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).parent
OUT = HERE / "example_data"
OUT.mkdir(exist_ok=True)

SEED = 42
np.random.seed(SEED)          # 固定全局种子(可复现)
rng = np.random.default_rng(SEED)

N_GENES = 60
N_CTRL = 300
N_STIM = 400

genes = [f"G{i:03d}" for i in range(N_GENES)]

# 真实扰动效应:只有前 20 个基因真正响应,其余为噪声基因
true_delta = np.zeros(N_GENES)
true_delta[:20] = rng.normal(0, 1.5, 20)

base = rng.uniform(1.0, 4.0, N_GENES)          # 每个基因的基线表达
noise_sd = rng.uniform(0.25, 0.6, N_GENES)     # 每个基因的细胞间噪声

ctrl = base[None, :] + rng.normal(0, noise_sd[None, :], (N_CTRL, N_GENES))
stim = base[None, :] + true_delta[None, :] + rng.normal(0, noise_sd[None, :], (N_STIM, N_GENES))

obs = pd.DataFrame(np.vstack([ctrl, stim]), columns=genes)
obs.insert(0, "group", ["control"] * N_CTRL + ["stimulated"] * N_STIM)
obs.insert(0, "cell_id", [f"ctrl_{i:04d}" for i in range(N_CTRL)] +
                         [f"stim_{i:04d}" for i in range(N_STIM)])
obs.to_csv(OUT / "observed.csv", index=False)

# ---- 三个待评估的"候选预测" ---------------------------------------------
# 评估集只用 stimulated 的后一半(见主脚本的 train/eval 切分),预测细胞数与之对齐即可
N_PRED = 200
preds = []

# goodModel:方向对、幅度略缩(0.85),噪声接近真实
d = true_delta * 0.85
X = base[None, :] + d[None, :] + rng.normal(0, noise_sd[None, :], (N_PRED, N_GENES))
preds.append(("goodModel", X))

# meanShiftModel:只给出常数均值位移,细胞间无变异(典型"只学到 pseudobulk"的模型)
X = np.tile(base + true_delta, (N_PRED, 1))
preds.append(("meanShiftModel", X))

# shuffledModel:把真实效应打乱到错误的基因上(方向/靶基因都错)
d = rng.permutation(true_delta)
X = base[None, :] + d[None, :] + rng.normal(0, noise_sd[None, :], (N_PRED, N_GENES))
preds.append(("shuffledModel", X))

frames = []
for name, X in preds:
    df = pd.DataFrame(X, columns=genes)
    df.insert(0, "method", name)
    df.insert(0, "cell_id", [f"{name}_{i:04d}" for i in range(N_PRED)])
    frames.append(df)
pd.concat(frames).to_csv(OUT / "predicted.csv", index=False)

print("wrote", OUT / "observed.csv", OUT / "predicted.csv")
