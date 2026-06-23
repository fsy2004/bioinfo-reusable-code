# 032 · 孟德尔随机化(MR)分析

> harmonized 工具变量数据 → 一条命令 → IVW/Egger/加权中位数 + 散点/森林/漏斗/留一图。

| | |
|---|---|
| **语言 / 主依赖** | R · `ggplot2`(MR 核心自包含,**不依赖 TwoSampleMR**) |
| **一句话用途** | 因果推断:暴露对结局的 MR 估计与敏感性分析 |
| **输入** | `example_data/harmonized_data.csv` |
| **输出** | `results/` 估计表+图 · 展示图见 `assets/` |

---

## ① 输入数据

harmonized CSV,需含:`SNP`, `beta.exposure`, `se.exposure`, `beta.outcome`, `se.outcome`(即 `TwoSampleMR::harmonise_data` 的输出格式;可由 028-031 的 GWAS 处理流程产生)。

## ② 方法 / 原理

- **IVW**(固定效应):过原点的逆方差加权回归 → 主因果估计。
- **MR-Egger**:带截距加权回归 → 截距检验水平多效性。
- **加权中位数**:对 Wald 比的加权中位数 → 对部分无效工具稳健。
- **敏感性**:留一法、漏斗图。

> 方法引用:Burgess *et al.*, *Eur J Epidemiol* 2017。核心 MR 自包含实现,便携、可离线运行。

## ③ 用途

用遗传工具变量推断暴露对结局的因果效应(规避混杂/反向因果),并以多方法 + 敏感性分析检验稳健性。

## ④ 特点 / 亮点

- **Turnkey + 便携**:harmonized 表即跑,不需安装/调用重型 TwoSampleMR、LD 服务。
- **顶刊图**:MR 散点(三方法线)· 森林(单SNP+合并)· 漏斗 · 留一法,全 theme_pub。

## ⑤ 输出结果图

| 文件 | 图型 | 说明 |
|------|------|------|
| `assets/MR_scatter.png` | 散点 | SNP 效应 + 三方法因果斜率 |
| `assets/MR_forest.png` | 森林 | 单 SNP Wald + 合并估计 |
| `assets/MR_funnel.png` · `MR_leaveoneout.png` | 漏斗/留一 | 多效性/稳健性 |

![scatter](assets/MR_scatter.png)

---

## 运行

```bash
Rscript 032_MR_analysis.R                                       # 示例
Rscript 032_MR_analysis.R --input data/harmonized.csv
```

## 依赖安装

```r
install.packages("ggplot2")   # 核心 MR 为自包含实现,无需 TwoSampleMR
```
