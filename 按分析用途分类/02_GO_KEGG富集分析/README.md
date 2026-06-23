# 02 · GO / KEGG 富集分析

把一组候选基因(差异基因、模块基因、机器学习特征、靶点交集等)注释到 GO 功能与 KEGG 通路。

| 模块 | 用途 | 语言 | 输出图 |
|------|------|------|--------|
| [007 GO/KEGG 富集](007_GO_KEGG富集分析/) | GO(BP/CC/MF)+ KEGG 通路富集 | R | GO 分面点图 · KEGG 棒棒糖 · 基因–通路网络 |

```bash
Rscript 007_GO_KEGG富集分析/007_GO_KEGG_enrichment.R                 # 跑示例
Rscript 007_GO_KEGG富集分析/007_GO_KEGG_enrichment.R --input data/genes.csv
```

> 遵循 [统一框架规范](../_framework/CONVENTIONS.md)。上游常接 03 类(差异基因)或 04 类(特征基因)。
