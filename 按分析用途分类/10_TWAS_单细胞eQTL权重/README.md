# 10 · TWAS / 单细胞 eQTL 权重

转录组关联分析(TWAS)的权重训练与模型拟合。

| 模块 | 用途 | 语言 | 状态 |
|------|------|------|:---:|
| 036-039 OneK1K TWAS 权重 | sc-eQTL 同/异质成分拟合、权重预处理与生成 | R | ⏭️ 重型 |
| 040-042 FUSION TWAS | FUSION targetC / S+targetC / S+allC 模型 | R | ⏭️ 重型 |

> ⏭️ **状态**:TWAS 流程依赖 FUSION、plink、大型 LD 参考面板与 sc-eQTL 权重文件,属计算密集 + 外部工具链,**本地未渲染**,保留原脚本作参考。
> 上游可接 09 类(GWAS 处理);权重生成后用 FUSION 官方流程做关联检验。配图规范见 [统一框架](../_framework/CONVENTIONS.md)。
