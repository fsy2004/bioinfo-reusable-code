# 19_多组学整合_分型模板

本模块用于多组学潜变量整合、监督多组学分类、共浸润 pattern 和无监督分型。

## 脚本

| 脚本 | 作用 |
|---|---|
| `083_MOFA_DIABLO_多组学潜变量整合.R` | 用 MOFA2 做无监督多组学因子分析，或用 mixOmics/DIABLO 做监督多组学整合。 |
| `084_NMF_ConsensusClusterPlus_共浸润分型.R` | 对表达、免疫评分或空间生态位矩阵做 NMF 和共识聚类。 |

## 推荐输入

- 多组学矩阵：每个矩阵行为 feature、列为 sample。
- metadata：样本分组、疾病状态、临床表型。
- 免疫浸润、空间生态位或通路活性矩阵。

## 推荐输出

- MOFA latent factors、feature weights。
- DIABLO sample variates。
- NMF metagene/pattern、样本亚型。
- consensus clustering 稳定性结果。
