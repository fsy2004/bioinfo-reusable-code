"""生成 571 模块的小型合成 CITE-seq 示例数据(仅演示用,非真实数据)。

只需运行一次,产物已提交到 example_data/。主脚本默认直接读这些 CSV。
合成设计:5 个细胞类型 → 每类特异表达一组基因;每个表面蛋白由其"同源基因"
(cognate gene)驱动 + 加性噪声,其中 4 个蛋白刻意设为弱耦合(转录本预测不了蛋白),
用来考察不同方法在"RNA-蛋白脱耦"位点上的差异。
"""
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).parent
OUT = HERE / "example_data"
OUT.mkdir(exist_ok=True)

rng = np.random.default_rng(20260720)

n_cells, n_genes, n_prot = 900, 80, 14
cell_types = ["CD4T", "CD8T", "Bcell", "Mono", "NK"]
labels = rng.choice(cell_types, n_cells, p=[0.3, 0.2, 0.15, 0.25, 0.10])

gene_names = [f"GENE{i:03d}" for i in range(n_genes)]
prot_names = ["CD3", "CD4", "CD8A", "CD19", "CD14", "CD56", "CD16", "CD25",
              "CD27", "CD45RA", "CD127", "HLA-DR", "PD1", "CD38"]
# 每个蛋白指定一个同源基因(前 14 个基因),其余基因为背景/细胞类型标记
cognate = {p: gene_names[i] for i, p in enumerate(prot_names)}
# 弱耦合蛋白:RNA 几乎不携带其蛋白信息(模拟 CD127/PD1 这类翻译后调控强的位点)
weak = {"CD127", "PD1", "CD25", "CD38"}

# --- 细胞类型潜在程序 ---------------------------------------------------------
type_idx = {t: i for i, t in enumerate(cell_types)}
prog = rng.normal(0, 1, (len(cell_types), n_genes)) * 1.2
latent = prog[[type_idx[t] for t in labels], :]
depth = rng.lognormal(0, 0.25, n_cells)[:, None]           # 测序深度差异
mu = np.exp(1.2 + latent) * depth
rna = rng.poisson(np.clip(mu, 0.02, 500)).astype(int)

# --- 蛋白:同源基因 + 类型效应 + 噪声 ----------------------------------------
rna_log = np.log1p(rna / rna.sum(1, keepdims=True) * 1e4)
prot = np.zeros((n_cells, n_prot))
for j, p in enumerate(prot_names):
    g = gene_names.index(cognate[p])
    w = 0.15 if p in weak else 1.0
    signal = w * rna_log[:, g] + (0.8 if p not in weak else 0.15) * latent[:, (g + 7) % n_genes]
    prot[:, j] = signal + rng.normal(0, 0.6, n_cells)
# 转回类 ADT 计数尺度(负二项式风格),保持非负整数
prot_counts = rng.poisson(np.clip(np.exp(prot + 2.0), 0.05, 5000)).astype(int)

NOTE = "# synthetic, for demo only -- not real CITE-seq data (571_captain_rna_protein_fm)\n"


def _dump(df, name, index=True):
    """写 CSV,首行加 `# synthetic` 注释(主脚本以 comment='#' 读取)。"""
    with open(OUT / name, "w", newline="") as fh:
        fh.write(NOTE)
        df.to_csv(fh, index=index)


_dump(pd.DataFrame(rna, index=[f"CELL{i:04d}" for i in range(n_cells)],
                   columns=gene_names), "citeseq_rna_counts.csv")
_dump(pd.DataFrame(prot_counts, index=[f"CELL{i:04d}" for i in range(n_cells)],
                   columns=prot_names), "citeseq_adt_counts.csv")
_dump(pd.DataFrame({"cell": [f"CELL{i:04d}" for i in range(n_cells)],
                    "cell_type": labels}), "cell_meta.csv", index=False)
_dump(pd.DataFrame({"protein": prot_names,
                    "cognate_gene": [cognate[p] for p in prot_names]}),
      "protein_gene_map.csv", index=False)
print("wrote synthetic CITE-seq demo to", OUT)
