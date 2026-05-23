# 13_转录因子调控_基因组圈图

本模块用于候选基因集的 TF/motif 调控解释、regulon 活性推断和染色体/圈图展示。

## 脚本

| 脚本 | 作用 |
|---|---|
| `047_RcisTarget_转录因子Motif网络.R` | 对候选基因集做 motif/TF 富集并构建调控网络。 |
| `053_circlize_基因染色体圈图.R` | 绘制目标基因在染色体上的圈图。 |
| `081_pySCENIC_Regulon_TF活性.py` | 调用 pySCENIC 的 GRN、ctx、AUCell 三步流程。 |

## 推荐链条

WGCNA/DEG/marker genes -> RcisTarget motif -> pySCENIC regulon -> AUCell/decoupler 活性 -> 轨迹、空间或扰动结果解释。
