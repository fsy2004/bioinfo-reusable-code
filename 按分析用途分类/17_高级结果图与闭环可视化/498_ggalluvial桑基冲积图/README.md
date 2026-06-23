# 498 · ggalluvial 桑基/冲积图

> 多层流向表(药物→hub→通路 等) → 一条命令 → 顶刊级冲积/桑基图。

| | |
|---|---|
| **语言 / 主依赖** | R · `ggalluvial` `ggplot2` |
| **一句话用途** | 展示多层级实体之间的流向关系 |
| **输入** | `example_data/flow_table.csv` |
| **输出** | `results/` · 展示图见 `assets/` |

---

## ① 输入数据

CSV 长表:前若干列 = 各层(2-4 层,如 `Drug`/`Hub`/`Pathway`),最后一列 `Freq`(流宽;缺省每行计 1)。

## ② 方法 / 原理

`ggalluvial::to_lodes_form` 转换 → `geom_flow`(sigmoid 曲带)+ `geom_stratum`(层块)+ 标签。

## ③ 用途

网络药理学"药物→hub→通路"、细胞通讯"配体→受体→细胞"等多层关系的经典展示(论文常见 Fig)。

## ④ 特点 / 亮点

- **Turnkey**:自动识别层数与频次列;2-4 层通用。
- **顶刊图**:期刊配色 + sigmoid 流带,信息密度高。

## ⑤ 输出结果图

| 文件 | 图型 | 说明 |
|------|------|------|
| `assets/Alluvial.png` | 冲积/桑基 | 多层流向 |

![alluvial](assets/Alluvial.png)

---

## 运行

```bash
Rscript 498_alluvial_sankey.R                              # 示例
Rscript 498_alluvial_sankey.R --input data/flow_table.csv
```

## 依赖安装

```r
install.packages(c("ggalluvial","ggplot2"))
```
