# 003 · CTD ∩ SwissTarget 化合物靶点 Venn

> 多份化合物靶点列表(CTD / SwissTargetPrediction) → 一条命令 → 交集/并集 + Venn + 集合柱状图。

| | |
|---|---|
| **语言 / 主依赖** | R · `theme_pub`(零依赖 Venn)+ `UpSetR` |
| **输入** | `example_data/`(CTD + SwissTarget 靶点 csv) |
| **输出** | `results/` 交并集表 + `assets/` 图 |

## ① 输入数据
`--input` 目录,内含 ≥2 份靶点列表(csv;自动识别 `Gene`/`Gene Symbol`/首列)。文件名=集合名。

## ② 方法 / 原理
读取各列表 → `Reduce(union/intersect)` 求并集/交集 → `venn_pub`(≤3集)+ 集合大小柱状图 + UpSet(≥3集)。

## ③ 用途
网络药理学中合并多数据库的化合物靶点,得到更全的候选靶点集(并集)或高置信靶点(交集)。

## ④ 特点 / 亮点
Turnkey;零依赖期刊级 Venn;自动识别基因列。

## ⑤ 输出结果图
| 文件 | 图型 |
|------|------|
| `assets/Target_Venn.png` | Venn |
| `assets/Set_size_bar.png` | 集合大小柱状 |

![Venn](assets/Target_Venn.png)

## 运行
```bash
Rscript 003_target_intersection_venn.R                       # 示例
Rscript 003_target_intersection_venn.R --input data/lists    # 你的目录
```
依赖:`install.packages("UpSetR")`(Venn 由 theme_pub.R 提供)
