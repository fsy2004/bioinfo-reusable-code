# 11 · WGCNA 共表达网络

| 模块 | 用途 | 语言 | 输出图 |
|------|------|------|--------|
| [054 WGCNA 共表达网络](054_WGCNA共表达网络/) | 软阈值+模块+模块-性状 | R | 无标度拟合 · 模块树状图 · 模块-性状热图 |

```bash
Rscript 054_WGCNA共表达网络/054_WGCNA_coexpression.R --input data/expr.csv --traits data/traits.csv
```

> turnkey;遵循 [统一框架规范](../_framework/CONVENTIONS.md)。关键模块基因可接 007 富集 / 047 TF 调控。
