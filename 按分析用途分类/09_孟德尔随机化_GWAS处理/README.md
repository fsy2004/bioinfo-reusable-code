# 09 · 孟德尔随机化 / GWAS 处理

GWAS 工具变量准备 → MR 因果推断 → 敏感性/方向性/中介分析。**[032](032_MR_TwoSampleMR分析/)** 是 turnkey 旗舰(自包含 MR,不依赖 TwoSampleMR)。

## ✅ turnkey 旗舰

| 模块 | 用途 | 语言 | 输出图 |
|------|------|------|--------|
| [032 MR 分析](032_MR_TwoSampleMR分析/) | IVW/Egger/WM + 敏感性 | R | 散点 · 森林 · 漏斗 · 留一 |

## 📦 GWAS 工具变量处理流程(上游 helper)

| 模块 | 作用 |
|------|------|
| 028 VCF 显著 SNP 筛选 · 029 LD clumping · 030 添加 EAF · 031 弱工具变量(F)过滤 | 从原始 GWAS 产出 harmonized 工具变量(→ 032 输入) |

## ⏭️ MR 进阶变体(保留参考)

| 模块 | 方法 | 依赖 |
|------|------|------|
| 033 MR 备用模板 | 基础 MR | 同 032 |
| 043 MR + 方向性 · 055 免疫细胞双向 MR | Steiger 方向性 | TwoSampleMR |
| 075 MR + coloc 因果证据链 | colocalization | coloc · LocusZoom |
| 079 pQTL MVMR 蛋白中介 | 多变量 MR | MVMR |
| 497 lavaan SEM 中介路径 | 结构方程 | lavaan |

> 旗舰 032 遵循 [统一框架规范](../_framework/CONVENTIONS.md);其核心 IVW/Egger/加权中位数为自包含实现,可直接套用到 033/043/055 的数据。进阶变体(coloc/MVMR/SEM)依赖专用包,保留原脚本作参考。
