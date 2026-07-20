"""生成 583 模块的合成示例数据 —— synthetic, for demo only.

产出的三个文件严格对齐 KEGNI 上游仓库的真实输入格式(逐一核对自
https://github.com/Lipxiao/KEGNI):

1. expression.csv          行=基因, 列=细胞, 第一列为基因名索引
                           (对齐 dataset/MAEDataset.py: pd.read_csv(input, index_col=0, header=0))
2. knowledge_graph.tsv     三列 TSV, 无表头: head \t relation \t tail
                           (对齐 dataloader/kge_dataloader.py: pd.read_csv(path, sep='\t', header=None))
                           不出现在表达矩阵里的节点被 KEGNI 视为 "kgg"(纯知识图节点)
3. ground_truth_network.csv  表头 Gene1,Gene2
                           (对齐 eval.py 与 data/GroundTruth/*/*-ChIP-network.csv)

设计意图(为了让基线对比有意义,而不是让某一路碾压):
- 一半真实边在表达上强相关(相关性基线抓得到);
- 另一半真实边被噪声淹没(相关性抓不到,只能靠知识先验);
- 知识图覆盖 3/4 真实边,但同时混入大量假边(先验单独用precision有限);
=> 表达 + 知识融合才应该最好。这正是 KEGNI 论文主张的动机。
"""
from __future__ import annotations
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SEED = 2024

N_TF, N_PER_TF, N_CELLS = 8, 4, 300
N_DECOY = 24                      # 不受 TF 调控、但被共同混杂因子带出假相关的基因


def build(seed: int = SEED):
    np.random.seed(seed)                 # 固定种子:合成数据必须逐字节可复现
    rng = np.random.default_rng(seed)
    tfs = [f"TF{i + 1:02d}" for i in range(N_TF)]
    targets = [f"TG{i + 1:03d}" for i in range(N_TF * N_PER_TF)]
    decoys = [f"DC{i + 1:03d}" for i in range(N_DECOY)]
    genes = tfs + targets + decoys

    # ---- 潜在 TF 活性(细胞维度的平滑信号)-------------------------------
    t = np.linspace(0, 1, N_CELLS)
    act = {}
    for k, tf in enumerate(tfs):
        phase = rng.uniform(0, 2 * np.pi)
        act[tf] = np.sin(2 * np.pi * (1 + k % 3) * t + phase) + 0.25 * rng.normal(size=N_CELLS)

    expr = {}
    edges = []                      # (TF, target, 是否表达可检测)
    for k, tf in enumerate(tfs):
        expr[tf] = 3.0 + 1.2 * act[tf] + 0.3 * rng.normal(size=N_CELLS)
        for j in range(N_PER_TF):
            tg = targets[k * N_PER_TF + j]
            strong = j < 2                       # 前两个靶点强、后两个被噪声淹没
            s = 1.10 if strong else 0.16
            noise = 0.35 if strong else 1.10
            expr[tg] = 3.0 + s * act[tf] + noise * rng.normal(size=N_CELLS)
            edges.append((tf, tg, strong))

    # ---- 混杂因子:让一批非靶基因也跟某些 TF 相关(制造假阳性)----------
    conf = np.cos(2 * np.pi * 2 * t) + 0.2 * rng.normal(size=N_CELLS)
    for i, dc in enumerate(decoys):
        lead = 0.9 if i < N_DECOY // 2 else 0.3
        expr[dc] = 3.0 + lead * conf + 0.5 * rng.normal(size=N_CELLS)
    for k in (0, 3, 5):                          # 三个 TF 也载荷混杂因子
        expr[tfs[k]] = expr[tfs[k]] + 0.8 * conf

    mat = pd.DataFrame({g: expr[g] for g in genes}).T
    mat.columns = [f"Cell{i + 1:04d}" for i in range(N_CELLS)]
    mat = mat.round(4)
    mat.index.name = None

    # ---- 知识图:TF -in_pathway-> PW <-in_pathway- target ----------------
    # 覆盖 3/4 真实边(含全部弱边),再掺入假边,先验单独用不够准。
    triples = []
    covered = [e for i, e in enumerate(edges) if i % 4 != 0]
    pw_of = {}
    for n, (tf, tg, _) in enumerate(covered):
        pw = f"PATHWAY{(n % 6) + 1:02d}"
        pw_of.setdefault(tf, pw)
        triples.append((tf, "in_pathway", pw_of[tf]))
        triples.append((tg, "in_pathway", pw_of[tf]))
    # 少量直接的基因-基因先验(真边)
    for tf, tg, _ in covered[:6]:
        triples.append((tf, "interacts_with", tg))
    # 假先验:decoy 与随机 target 也被塞进同样的 pathway
    for dc in decoys:
        triples.append((dc, "in_pathway", f"PATHWAY{rng.integers(1, 7):02d}"))
    for _ in range(40):
        a, b = rng.choice(targets, 2, replace=False)
        triples.append((str(a), "interacts_with", str(b)))
    # 知识图里还有一些与表达矩阵无关的纯 KG 节点(KEGNI 的 kgg-kgg 三元组)
    for i in range(6):
        triples.append((f"PATHWAY{i + 1:02d}", "part_of", "PATHWAY_ROOT"))

    kg = pd.DataFrame(sorted(set(triples)))
    gt = pd.DataFrame([(a, b) for a, b, _ in edges], columns=["Gene1", "Gene2"])
    return mat, kg, gt


if __name__ == "__main__":
    mat, kg, gt = build()
    mat.to_csv(os.path.join(HERE, "expression.csv"))
    kg.to_csv(os.path.join(HERE, "knowledge_graph.tsv"), sep="\t",
              header=False, index=False)
    gt.to_csv(os.path.join(HERE, "ground_truth_network.csv"), index=False)
    print(f"expression {mat.shape}  kg {kg.shape}  ground-truth edges {len(gt)}")
