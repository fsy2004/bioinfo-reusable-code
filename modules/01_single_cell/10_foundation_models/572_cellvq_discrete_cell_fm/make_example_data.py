"""生成 example_data/ 的小型合成 scRNA 计数矩阵(synthetic, for demo only)。

只需跑一次,产物已提交进库;主脚本默认直接读该 CSV,不依赖本文件。
"""
import numpy as np, pandas as pd, os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "example_data")

def main(n_per_type=110, n_genes=160, seed=572):
    rng = np.random.default_rng(seed)
    types = ["Tcell", "Bcell", "Mono", "NK", "DC"]
    n_types = len(types)
    # 每种细胞有一段特异高表达基因块,其余基因共享背景 → 可分但有重叠
    base = rng.gamma(2.0, 1.2, size=n_genes)          # 基础表达水平
    blocks = np.array_split(np.arange(n_genes), n_types)

    X, labels = [], []
    for ti, t in enumerate(types):
        mu = np.tile(base, (n_per_type, 1))
        # marker 上调幅度刻意压小 → 类型间有真实重叠,连续嵌入不会顶到 100%,
        # 否则天花板效应会让"离散化损失多少信息"这个比较失去分辨力
        mu[:, blocks[ti]] *= rng.uniform(1.5, 1.9)
        # 相邻类型共享一部分程序(如 T/NK),制造连续过渡而非干净团块
        mu[:, blocks[(ti + 1) % n_types]] *= 1.15
        # 细胞内测序深度差异(真实数据的主要技术噪声源)
        depth = rng.lognormal(0.0, 0.55, size=(n_per_type, 1))
        counts = rng.poisson(np.clip(mu * depth, 0.02, None))
        X.append(counts)
        labels += [t] * n_per_type

    X = np.vstack(X).astype(int)
    genes = [f"G{i:03d}" for i in range(n_genes)]
    cells = [f"C{i:04d}" for i in range(X.shape[0])]
    df = pd.DataFrame(X, index=cells, columns=genes)
    df.insert(0, "cell_type", labels)
    df.index.name = "cell_id"

    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, "synthetic_counts.csv")
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write("# synthetic, for demo only -- 5 cell types x 110 cells x 160 genes, seed=572\n")
        df.to_csv(fh)
    print("wrote", path, df.shape)

if __name__ == "__main__":
    main()
