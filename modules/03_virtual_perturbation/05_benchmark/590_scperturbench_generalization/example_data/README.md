# example_data — **synthetic, for demo only**

这两个 CSV 由同目录上一层的 `make_example_data.py` 生成(seed=42),**不是真实实验数据**,
仅用于本模块的冒烟测试与 README 展示图。任何结论都不具生物学意义。

| 文件 | 内容 |
|---|---|
| `observed.csv` | 真实观测:300 control + 400 stimulated 细胞 × 60 基因 |
| `predicted.csv` | 三个人造"候选模型"的预测:`goodModel` / `meanShiftModel` / `shuffledModel`,各 200 细胞 |

设计意图:`meanShiftModel` 只给出正确的**均值位移**、细胞间变异被压扁 —— 它在
`mse` / `pearson_distance` / `common_deg` 上都排第 1,却在 `sym_kldiv` 上垫底(第 5)。
单看一个指标会得出完全相反的结论,这正是 scPerturBench 想让人看见的东西。

另注:按上游 `trainMean.py` 口径实现的朴素基线 `trainMean`,在这份示例上**平均排名第一
(1.75)**,压过两个"模型" —— 与论文结论同向。示例是合成的,该现象只作演示,不构成证据。
