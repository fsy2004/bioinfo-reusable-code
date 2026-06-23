# 13 · 转录因子调控 / 基因组圈图

候选基因集的 TF/motif 调控解释、regulon 活性推断与染色体圈图展示。

| 模块 | 用途 | 语言 | 输出图 | 状态 |
|------|------|------|--------|:---:|
| [053 染色体圈图](053_circlize染色体圈图/) | 基因染色体位置 circos | R | 圈图 | ✅ turnkey |
| 047 RcisTarget motif-TF 网络 | motif/TF 富集 + 调控网络 | R | 网络 · Sankey | ⏭️ 需 cisTarget DB |
| 081 pySCENIC regulon | GRN+ctx+AUCell,TF 活性 | Python | UMAP · 热图 | ⏭️ 重型环境 |

> **053** turnkey:坐标表即出顶刊圈图(遵循 [统一框架规范](../_framework/CONVENTIONS.md))。
> **047**:需 RcisTarget motif 排名数据库(GB 级);**081**:pySCENIC(Python + GRNBoost,重型)。两者保留原脚本作参考。

## 推荐链条
WGCNA/DEG/marker genes → RcisTarget motif → pySCENIC regulon → AUCell/decoupler 活性 → 轨迹/空间/扰动解释。
