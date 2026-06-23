# 022 · 分子对接结合能可视化

> 化合物×靶点结合能矩阵 → 一条命令 → 结合能热图(蓝=强结合)+ 最强结合排序气泡图。

| | |
|---|---|
| **语言 / 主依赖** | R · `ComplexHeatmap` `ggplot2` |
| **一句话用途** | 对接结合能的矩阵热图与排序可视化 |
| **输入** | `example_data/binding_energy.csv` |
| **输出** | `results/` 表+图 · 展示图见 `assets/` |

---

## ① 输入数据

CSV:首列 `Target`(靶点名),其余列 = 各化合物的对接结合能(kcal/mol,越负结合越强)。单化合物时只需一列。

## ② 方法 / 原理

结合能矩阵 → ComplexHeatmap 聚类热图(蓝→白渐变,标注数值)→ 每靶点取最强结合化合物,排序气泡图。

## ③ 用途

网络药理学/对接验证:直观比较多个化合物对多个靶点的亲和力,锁定强结合的化合物-靶点对(通常 < −7 kcal/mol 视为良好结合)。

## ④ 特点 / 亮点

- **Turnkey**:矩阵即跑;自动聚类、数值标注、强结合配色。
- **顶刊图**:对接结合能热图 + 最强结合气泡图。

## ⑤ 输出结果图

| 文件 | 图型 | 说明 |
|------|------|------|
| `assets/Binding_heatmap.png` | 热图 | 化合物×靶点结合能 |
| `assets/Binding_bubble.png` | 气泡/棒棒糖 | 每靶点最强结合 |

![heatmap](assets/Binding_heatmap.png)

---

## 运行

```bash
Rscript 022_docking_binding_energy.R                              # 示例
Rscript 022_docking_binding_energy.R --input data/binding_energy.csv
```

## 依赖安装

```r
install.packages(c("ggplot2","circlize")); BiocManager::install("ComplexHeatmap")
```
