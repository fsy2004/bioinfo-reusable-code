# 08 · 单细胞 / 空间转录组 / 细胞轨迹

本类别体量最大。**[046](046_单细胞发表级图/)** 是 turnkey 旗舰(Seurat 标准流程 + 发表级标配图,零改动即跑出真图);其余多为重型/专用流程(CellChat / monocle / Scissor / 空间 / 深度学习),保留原脚本作引擎与参考。

## ✅ turnkey 旗舰

| 模块 | 用途 | 语言 | 输出图 |
|------|------|------|--------|
| [046 单细胞发表级图](046_单细胞发表级图/) | Seurat 全流程 + 标配图 | R | UMAP · marker点图 · marker热图 · FeaturePlot · violin |

## 📦 引擎 / 数据前处理

| 模块 | 作用 |
|------|------|
| 023 RDS 对象结构检查 · 024/025 单细胞 QC 整理 | 数据读取/质控/整理(供 046 等上游) |

## ⏭️ 重型 / 专用(保留参考)

| 模块 | 方法 | 为何未本地渲染 |
|------|------|----------------|
| 026 Seurat 全流程 + scTenifoldKnk 敲除 | 虚拟敲除 | scTenifoldKnk 重型 |
| 049 手工注释 + CellChat + monocle | 细胞通讯+轨迹 | CellChat/monocle 重型 |
| 051 CellChat 细胞通讯 | 圆图/chord/气泡 | CellChat 重型 |
| 058 Scissor 疾病相关细胞 | 表型关联细胞 | Scissor + 队列 |
| 044 AD 单细胞 + monocle | 拟时序 | monocle 重型 |
| 027 / 050 空间转录组 | Visium 空间分析 | 需空间数据 |
| 062 scTour · 082 Palantir/Slingshot | 轨迹/向量场 | Python 深度学习 / 重型 |
| 061 scFOCAL/CellOracle GUI 准备 · 491 scTour 环境 | 扰动入口/环境 | 外部 GUI/环境 |

> 旗舰 046 遵循 [统一框架规范](../_framework/CONVENTIONS.md)。重型模块如需运行,见各脚本头部依赖说明。
