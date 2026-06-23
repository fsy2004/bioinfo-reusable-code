# 015 · 多方法特征交集 Venn / UpSet

> 一个目录里的多份基因列表 → 一条命令 → 全局/两两交集 + Venn(≤3集)+ UpSet 图。

| | |
|---|---|
| **语言 / 主依赖** | R · `theme_pub`(零依赖 Venn)+ `UpSetR` |
| **一句话用途** | 多 ML 方法特征求交、可视化一致性 |
| **输入** | `example_data/gene_sets/`(3 份基因列表) |
| **输出** | `results/` 交集表+图 · 展示图见 `assets/` |

---

## ① 输入数据

**输入是一个目录**(`--input`),内含 ≥2 份基因列表(csv/txt,首列=基因名;文件名 = 集合名)。

## ② 方法 / 原理

读取各列表 → `Reduce(intersect)` 求全局交集 + 两两交集 → Venn(≤3 集用零依赖 `venn_pub`)与 UpSet(任意集合数)。

## ③ 用途

LASSO/RF/SVM-RFE 等多法筛选后取**一致特征基因**(交集 = 最稳健候选),供诊断/预后建模。

## ④ 特点 / 亮点

- **Turnkey + 零依赖 Venn**:不需 ggvenn/VennDiagram,`venn_pub` 直接出期刊级 Venn。
- 全局交集 + 两两交集表 + UpSet,一次到位。

## ⑤ 输出结果图

| 文件 | 图型 | 说明 |
|------|------|------|
| `assets/Feature_Venn.png` | Venn | ≤3 集合时的交集 Venn |
| `assets/Feature_UpSet.png` | UpSet | 任意集合数的交集条形 |
| `results/global_intersection.csv` · `pairwise_intersection.csv` | 表 | 交集基因 |

![Venn](assets/Feature_Venn.png)

---

## 运行

```bash
Rscript 015_feature_intersection.R                          # 示例(3 集)
Rscript 015_feature_intersection.R --input data/gene_sets   # 你的目录
```

## 依赖安装

```r
install.packages("UpSetR")   # Venn 由共享 theme_pub.R 提供,无需额外包
```
