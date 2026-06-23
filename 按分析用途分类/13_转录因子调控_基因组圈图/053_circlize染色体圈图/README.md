# 053 · circlize 基因染色体圈图

> 基因坐标表 → 一条命令 → 顶刊级染色体圈图(ideogram + 基因标签)。

| | |
|---|---|
| **语言 / 主依赖** | R · `circlize` |
| **一句话用途** | 把目标基因画到染色体上展示分布 |
| **输入** | `example_data/gene_positions.csv` |
| **输出** | `results/` · 展示图见 `assets/` |

---

## ① 输入数据

CSV,列:`Gene`, `Chr`(如 `chr17`), `Start`, `End`(基因组坐标,可从 NCBI Gene/UCSC 获取)。

## ② 方法 / 原理

`circos.initializeWithIdeogram` 初始化染色体框架 → 外圈彩色染色体轨道 + 中圈 cytoband ideogram → `circos.genomicLabels` 内圈基因标签引线。

## ③ 用途

把候选基因集(差异/特征/靶点基因)的染色体分布做成一区风格补充图,直观展示基因组定位与聚集。

## ④ 特点 / 亮点

- **Turnkey**:单坐标表即跑;期刊配色染色体。
- **顶刊图**:三层圈图(染色体 + cytoband + 基因标签)。
- 支持 `--genome hg38/hg19/mm10` 等(首次联网取 cytoband)。

## ⑤ 输出结果图

| 文件 | 图型 | 说明 |
|------|------|------|
| `assets/Chromosome_circos.png` | circos 圈图 | 基因染色体分布 |

![circos](assets/Chromosome_circos.png)

---

## 运行

```bash
Rscript 053_chromosome_circos.R                              # 示例(hg38)
Rscript 053_chromosome_circos.R --input data/gene_positions.csv --genome hg38
```

## 依赖安装

```r
install.packages("circlize")
```
