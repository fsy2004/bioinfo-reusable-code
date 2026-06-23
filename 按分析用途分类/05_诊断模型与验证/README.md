# 05 · 诊断模型与验证

把筛选出的特征基因(04 类)落地成 logistic 诊断模型,并做内部 + 外部验证。

| 模块 | 用途 | 语言 | 输出图 |
|------|------|------|--------|
| [016 诊断模型](016_诊断模型_ROC校准DCA/) | 建模 + 内部评价 | R | 列线图 · 校准 · DCA · ROC · OR森林 · 箱线 |
| [063 外部验证](063_GEO诊断模型验证/) | 独立队列验证 | R | 训练/验证 ROC · 校准 |

```bash
Rscript 016_诊断模型_ROC校准DCA/016_diagnostic_model.R                      # 内部:建模+评价
Rscript 063_GEO诊断模型验证/063_diagnostic_validation.R                     # 外部:独立队列验证
```

> 遵循 [统一框架规范](../_framework/CONVENTIONS.md)。上游接 04 类特征基因。
