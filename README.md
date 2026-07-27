# Reusable Bioinformatics Code Library · 可复用生物信息学代码库

> **A domain-organised library of reusable R and Python analyses.**
> Find an analysis by the question you want to answer, open the matching module,
> and use its README for the exact input, method, dependencies, and output figures.
>
> **按分析问题组织的 R / Python 生信模块库。** 先按「我要做什么分析」定位域和子类,
> 再打开具体模块的 README,查看输入格式、分析方法、依赖与输出图。

<p>
<img alt="modules" src="https://img.shields.io/badge/modules-196-blue">
<img alt="domains" src="https://img.shields.io/badge/domains-10-informational">
<img alt="languages" src="https://img.shields.io/badge/R%204.4%20%7C%20Python%203.12-informational">
<img alt="figures" src="https://img.shields.io/badge/figures-vector%20PDF%20%2B%20300dpi%20PNG-success">
</p>

---

## 这套库能做什么？ · What analyses are included?

仓库目前包含 **196 个模块**,按 **10 个分析域、53 个子类**组织。下面是首页的
「分析导航」；每个域的链接会进入该域 README,完整的逐模块索引见
[`modules/CATALOG.md`](modules/CATALOG.md)。

| 域 | 主要分析内容 | 进入目录 |
|---|---|---|
| **01 单细胞分析** | 数据读入与质控、批次整合、细胞注释与分型、组成差异、差异表达/pseudobulk、轨迹与 RNA velocity、CNV/克隆、通路与 TF 活性、单细胞与 bulk 表型关联、基础模型 | [`modules/01_single_cell`](modules/01_single_cell) |
| **02 空间转录组** | 细胞分割、空间域与 SVG、空间统计、解卷积与单细胞映射、切片配准与三维重建、空间细胞通讯、空间多组学、空间基础模型 | [`modules/02_spatial_transcriptomics`](modules/02_spatial_transcriptomics) |
| **03 虚拟扰动** | 虚拟敲除、基因调控网络推断、因果表示与反事实、药物扰动与响应、扰动预测基准 | [`modules/03_virtual_perturbation`](modules/03_virtual_perturbation) |
| **04 因果推断与遗传** | 工具变量准备、两样本 MR、cis-MR 与药靶、中介与 MVMR、共定位、TWAS/sc-eQTL、稳健 MR 估计量 | [`modules/04_causal_inference_genetics`](modules/04_causal_inference_genetics) |
| **05 机器学习** | 特征筛选、分类模型、生存机器学习、模型解释、不确定性量化、泛化与外部验证 | [`modules/05_machine_learning`](modules/05_machine_learning) |
| **06 Bulk 组学** | 差异表达、GO/KEGG 与 GSEA、WGCNA/共表达网络、多组学整合与分型、突变/甲基化/蛋白质组/代谢组 | [`modules/06_bulk_omics`](modules/06_bulk_omics) |
| **07 临床与转化** | 诊断模型、预后与生存、免疫浸润与解卷积、药物警戒、疾病负担与人群队列 | [`modules/07_clinical_translational`](modules/07_clinical_translational) |
| **08 结构与药物设计** | 分子对接、分子动力学、虚拟筛选与打分 | [`modules/08_structure_drug_design`](modules/08_structure_drug_design) |
| **09 网络药理学** | 靶点数据库提取、靶点交集与集合关系、成药性评分 | [`modules/09_network_pharmacology`](modules/09_network_pharmacology) |
| **10 可视化** | 高级科研图型、出版级图形模板、外部绘图资源 | [`modules/10_visualization`](modules/10_visualization) |

### 按研究问题找入口

- 想回答「哪些细胞、状态或亚群发生变化」：从 **01 单细胞** 或 **02 空间转录组**开始。
- 想回答「哪些基因、通路或调控因子发生变化」：查看 **01 单细胞**、**06 Bulk 组学** 的差异表达、富集与活性分析。
- 想回答「细胞之间如何通讯、空间上是否形成生态位」：查看 **02 空间转录组** 的空间统计、细胞通讯与空间多组学。
- 想回答「某个基因/通路/药物靶点是否可能具有因果作用」：查看 **03 虚拟扰动**、**04 因果推断与遗传**、**09 网络药理学**。
- 想建立「诊断、预后或风险预测模型」：查看 **05 机器学习** 与 **07 临床与转化**。
- 想回答「候选小分子如何结合、是否稳定、如何筛选」：查看 **08 结构与药物设计**。
- 想把结果整理成论文图：查看 **10 可视化**,或直接复用各分析模块的 `assets/` 与统一绘图框架。

---

## 模块怎么用？ · How a module works

每个模块是一个相对独立的分析目录,通常包含：

```text
<domain>/<subcategory>/<NNN_module>/
├── <NNN_module>.R | .py     # 主脚本
├── README.md                # 输入、方法、依赖、输出与限制
├── example_data/            # 小型示例或合成数据（如适用）
└── assets/                  # 已渲染的 PNG/PDF 预览图（如适用）
```

推荐顺序：

1. 先读模块 README,确认输入对象、样本/细胞层级、依赖和适用范围。
2. 先用 `example_data/` 跑通示例,再用 `--input` / `--outdir` 指向自己的数据。
3. 对基础模型、服务器方法或外部工具链,同时查看状态标记和
   [`SERVER_DEPENDENCIES.md`](modules/_framework/SERVER_DEPENDENCIES.md)。

---

## Quick start · 快速开始

```bash
git clone https://github.com/fsy2004/bioinfo-reusable-code.git
cd bioinfo-reusable-code

# 进入对应域的目录,再按模块 README 运行
cd modules
Rscript 06_bulk_omics/01_differential_expression/010_geo_deg_volcano_heatmap_pca/010_*.R
python 02_spatial_transcriptomics/02_domains_svg_stats/543_squidpy_spatial_statistics/543_*.py
```

模块的具体命令和参数以各自 README 为准。默认目标环境为 **R 4.4** 和
**Python 3.12**；不同模块的依赖可能需要本地、服务器或 GPU 环境。

---

## 状态标记 · Status legend

| 标记 | 含义 |
|---|---|
| ✅ | 自带/合成示例可本机运行,无需修改 |
| 🟡 | 真实基线或核心流程可本机运行,完整方法需安装服务器依赖 |
| 🔴 | 重型、GPU 或外部工具链；保留守卫式封装或参考实现 |
| 📄 | 模板或上游脚本；需要自行准备数据和依赖 |
| 📦 | Vendored/外部工具清单或本地来源,不是完整 turnkey 模块 |
| 🗃️ | 仅本地保留,由 Git 忽略,不会发布到仓库 |

“可运行”不等于“适合所有数据”。请优先检查每个模块 README 中的输入限制、样本层级、
基线、随机种子和输出解释。

---

## 统一框架 · Shared framework

[`modules/_framework/`](modules/_framework/) 提供跨模块的约定与工具：

- `CATALOG.md`：按域 → 子类的完整模块索引，包含用途、输入→输出、依赖、语言、图型与状态。
- `TOOL_SELECTION_GUIDE.md`：按分析任务选择方法和模块。
- `CONVENTIONS.md`：目录结构、复用优先、路径、随机种子和脚本约定。
- `theme_pub.R` / `pubstyle.py`：统一配色、主题和 PDF + 300 dpi PNG 导出。
- `ANALYSIS_TEMPLATE/`：带配置、检查点和环境快照的新项目脚手架。
- `QUALITY_CHECKLIST.md` / `qc_lint.py`：分析质量清单与静态检查。
- `SERVER_DEPENDENCIES.md`：服务器/GPU/外部工具链依赖说明。
- `BENCHMARKS_AND_CRITIQUES.md`：复杂模型、基础模型和扰动预测的比较尺度与常见陷阱。

---

## 可复现约定 · Reproducibility

- 优先复用真实工具和现有模块,不凭记忆重写分析 API。
- 使用相对路径、固定随机种子和模块级环境说明；避免硬编码本机路径。
- 运行输出默认写入 `results/`，本地运行结果、缓存和上游源码不作为仓库内容发布。
- 图形优先输出矢量 PDF 与 300 dpi PNG；复杂模型尽量同时报告透明、可比较的基线。
- 结论应以模块实际输出和数据层级为准，不把示例数据结果当作真实研究结论。

---

## License · 许可

Each module follows the license of the tools and methods it uses. Vendored third-party code
keeps its original license; see the relevant module README and upstream repository.
每个模块遵循所用工具与方法的许可。第三方代码保留原始许可证,详见对应模块 README 与上游仓库。
