# CNV 分析模板：GISTIC2 / CNVkit / inferCNV

本文件是 CNV 分析入口说明，当前不直接绑定某一个外部软件环境。

## 适用数据

| 数据类型 | 推荐工具 | 输出 |
|---|---|---|
| bulk WES/WGS CNV segment | GISTIC2 | peak、amplification/deletion genes |
| panel/WES target coverage | CNVkit | segment、gene-level CNV、scatter/diagram |
| scRNA-seq 推断 CNV | inferCNV | cell-level CNV heatmap |

## 标准输入

- sample metadata：`sample, group, batch`
- segment 文件：`sample, chrom, start, end, log2/copy_number`
- gene annotation：可选，用于 gene-level 汇总

## 推荐输出

- `cnv_segments_standardized.tsv`
- `cnv_gene_level.tsv`
- `cnv_recurrent_peaks.tsv`
- `cnv_group_comparison.tsv`
- `figures/cnv_heatmap.pdf`
- `figures/cnv_oncoprint.pdf`

## 文章链条位置

候选基因或分型 -> CNV 支持 -> 表达/CNV 相关 -> 生存、免疫或药敏关联。

## 后续封装建议

1. 对 GISTIC2：写 shell/PowerShell wrapper，固定 marker、seg、refgene 输入。
2. 对 CNVkit：写 batch wrapper，输出统一 segment 表。
3. 对 inferCNV：接 Seurat 注释，固定 reference cell type 和 tumor/target cell type。
