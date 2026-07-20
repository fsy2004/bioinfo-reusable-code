# example_data — synthetic, for demo only

本目录三个文件均为**合成数据**,由固定种子(`set.seed(562)`)生成,只用于冒烟测试与生成
README 展示图。**不是真实实验数据,不可用于任何生物学结论。**

| 文件 | 内容 |
|---|---|
| `counts.csv` | 700 基因 × 235 细胞的原始 count 矩阵(首列 `Gene` 为基因名) |
| `cell_metadata.csv` | 每个细胞的真实类型标签(T / B / Mono / Rare) |
| `ground_truth_hvg.csv` | 90 个"真实高变基因"及其表达档位(HI / MID / LO) |

## 生成逻辑(为什么这样造)

数据刻意做成**不同 HVG 方法各有盲区**的样子,因为这正是 mixhvg 论文里"混合优于单一"
的前提。若数据太容易,所有方法都满分,比较就没有意义。

- **HI / MID / LO 三档 marker**(各 30 个):高 / 中 / 低平均表达。vst 在高表达端把
  marker 压得很低,dispersion 类反之 —— 盲区互补。
- **Rare 稀有群**(仅 25 细胞):它的 marker 有 2/3 落在 LO 档。漏掉低表达 marker 的方法
  在下游聚类里就分不出这个群,silhouette 会掉下来。
- **TRAP 基因**(90 个):高表达 + 轻度过离散(负二项 size=20),但**没有任何分组结构**,
  是均值-方差类方法的经典假阳性来源。
- **BG 背景基因**(520 个)+ 细胞测序深度差异(lognormal),模拟真实数据的噪声底。

fold change 被刻意压小(HI 1.9× / MID 2.3× / LO 6.0×),让 recall 不饱和。
