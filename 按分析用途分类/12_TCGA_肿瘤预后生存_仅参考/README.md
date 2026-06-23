# 12 · TCGA 肿瘤预后生存(仅参考)

肿瘤队列预后分析(单基因 / 风险签名)。非肿瘤研究若有随访数据亦可套用。

| 模块 | 用途 | 语言 | 输出图 | 状态 |
|------|------|------|--------|:---:|
| [048 单基因多终点生存](048_TCGA单基因生存曲线/) | OS/DSS/DFI/PFI KM | R | 4 终点 KM | ✅ |
| [057 预后风险模型](057_TCGA预后风险模型/) | 风险签名五件套 | R | 风险分布·状态·热图·KM·时间ROC | ✅ |
| [060 免疫双蝴蝶图](060_TCGA免疫双蝴蝶图/) | 基因-免疫相关蝴蝶图 | R | butterfly | ✅ |
| 497 scSurvival | 单细胞-队列联合生存(外部包) | Python | 队列生存 | ⏭️ 外部包 |

```bash
Rscript 048_TCGA单基因生存曲线/048_single_gene_survival.R --gene TP53   # 单基因
Rscript 057_TCGA预后风险模型/057_prognostic_risk_model.R                  # 风险签名
```

> 048/057 turnkey(遵循 [统一框架规范](../_framework/CONVENTIONS.md))。497_scSurvival 为外部 Python 包,保留作参考。
