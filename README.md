# Reusable Bioinformatics Code Library

按研究问题组织的 R/Python 生物信息学代码库：先选领域，再进入可独立运行的模块，用自己的数据替换合成示例。

<!-- readme-metrics:start -->
![modules](https://img.shields.io/badge/modules-188-2563EB)
![catalog entries](https://img.shields.io/badge/catalog_entries-202-334155)
![domains](https://img.shields.io/badge/domains-15-0F766E)
![subcategories](https://img.shields.io/badge/subcategories-59-7C3AED)
![outputs](https://img.shields.io/badge/figures-PDF%20%2B%20300dpi%20PNG-15803D)
<!-- readme-metrics:end -->

## 这是什么

本库保存可复用的分析入口、输入样例、来源说明和投稿级展示图。189 个三位编号模块构成主要分析单元；完整目录另含 14 个模板、索引与本地资源入口，因此总目录为 203 条。数字由文件树自动生成，不手工维护。

每个新模块必须满足三条底线：从真实工具或已发表方法起步；默认示例可以零改动运行；示例成功只证明接口和审计路径可工作，不代替真实数据分析或科学验证。

## Quick start

先在 [`modules/CATALOG.md`](modules/CATALOG.md) 按领域、用途或图型定位模块，再阅读该模块的 README。

```powershell
git clone https://github.com/fsy2004/bioinfo-reusable-code.git
cd bioinfo-reusable-code

# 示例：审计 ClairS 兼容 VCF；默认读取合成数据
python modules\12_genome_annotation_variants\01_somatic_variants\598_clairs_longread_somatic_audit\598_clairs_longread_somatic_audit.py
```

默认结果写入模块自己的 `results/`；展示用 PNG 保存在 `assets/`。R 模块通常以 `Rscript <module>.R` 运行，Python 模块以 `python <module>.py` 运行。完整参数、依赖和真实输入规格以模块 README 为准。

## 按问题选领域

- 细胞状态、亚群、轨迹或空间生态位：01–03。
- 遗传因果、QTL、预测模型和 bulk 组学：04–06。
- 临床转化、药物、结构与网络：07–09。
- 投稿图与证据合成：10–11。
- 基因组变异/注释、比较基因组、微生物组、表观组和蛋白质组：12–16。

<!-- readme-domains:start -->
| Domain | Entries | Entry |
|---|---:|---|
| 01 · 单细胞分析 — Single-cell analysis | 33 | [`modules/01_single_cell`](modules/01_single_cell) |
| 02 · 空间转录组 — Spatial transcriptomics | 24 | [`modules/02_spatial_transcriptomics`](modules/02_spatial_transcriptomics) |
| 03 · 虚拟扰动技术 — Virtual perturbation | 26 | [`modules/03_virtual_perturbation`](modules/03_virtual_perturbation) |
| 04 · 因果推断与遗传流行病 — Causal inference & genetics | 29 | [`modules/04_causal_inference_genetics`](modules/04_causal_inference_genetics) |
| 05 · 机器学习 — Machine learning | 18 | [`modules/05_machine_learning`](modules/05_machine_learning) |
| 06 · Bulk 组学 — Bulk omics | 19 | [`modules/06_bulk_omics`](modules/06_bulk_omics) |
| 07 · 临床与转化 — Clinical & translational | 20 | [`modules/07_clinical_translational`](modules/07_clinical_translational) |
| 08 · 结构与药物设计 — Structure & drug design | 6 | [`modules/08_structure_drug_design`](modules/08_structure_drug_design) |
| 09 · 网络药理学 — Network pharmacology | 8 | [`modules/09_network_pharmacology`](modules/09_network_pharmacology) |
| 10 · 可视化 — Visualization | 13 | [`modules/10_visualization`](modules/10_visualization) |
| 12 · 基因组注释与变异 — Genome annotation & variants | 2 | [`modules/12_genome_annotation_variants`](modules/12_genome_annotation_variants) |
| 13 · 比较基因组与系统发育 — Comparative genomics & phylogenetics | 1 | [`modules/13_comparative_genomics_phylogenetics`](modules/13_comparative_genomics_phylogenetics) |
| 14 · 微生物组与宏基因组 — Microbiome & metagenomics | 1 | [`modules/14_microbiome_metagenomics`](modules/14_microbiome_metagenomics) |
| 15 · 表观基因组与染色质 — Epigenomics & chromatin | 1 | [`modules/15_epigenomics_chromatin`](modules/15_epigenomics_chromatin) |
| 16 · 蛋白质组与代谢组 — Proteomics & metabolomics | 1 | [`modules/16_proteomics_metabolomics`](modules/16_proteomics_metabolomics) |
<!-- readme-domains:end -->

## 模块契约

```text
<domain>/<subcategory>/<NNN_module>/
├── <NNN_module>.R | .py     # 主入口
├── README.md                # 输入、方法、用途、特点、输出与边界
├── example_data/            # 小型合成或公开样例
├── assets/                  # README 展示图
└── results/                 # 运行生成；Git 忽略
```

状态、服务器依赖和选型建议分别见：

- [`modules/CATALOG.md`](modules/CATALOG.md)：全部条目、状态、I/O、依赖和图型。
- [`modules/_framework/TOOL_SELECTION_GUIDE.md`](modules/_framework/TOOL_SELECTION_GUIDE.md)：按分析问题选工具。
- [`modules/_framework/SERVER_DEPENDENCIES.md`](modules/_framework/SERVER_DEPENDENCIES.md)：重型/GPU/上游环境。
- [`modules/_framework/CONVENTIONS.md`](modules/_framework/CONVENTIONS.md)：模块工程与科研质量边界。

## 2026 方法更新

2026-08-20 冻结扫描把原有 11 个领域扩展到 16 个，新增 ClairS、EviAnn、OrthoFinder v3、CoverM、damidBind 和 DiaReport 的轻量输出审计模块。完整上游方法仍在其官方环境运行；本地模块负责稳定 I/O、错误拦截、透明摘要和展示图。

来源、许可、观察区和下一轮范围见 [`METHOD_SCAN_2026_08_20.md`](modules/_framework/METHOD_SCAN_2026_08_20.md)。该登记表覆盖当前全领域框架，但不是系统综述，也不声称穷尽所有 2026 方法。

## 科学边界

- 合成数据用于接口回归和图形演示，不产生生物学结论。
- 供体/样本层推断、协变量、批次、随机种子和真实数据来源必须在具体项目中重新确认。
- 基础模型、虚拟扰动和药物预测用于假设生成；关键结论需要简单基线、外部数据或实验支持。
- 第三方代码、数据库和方法各自遵循上游许可；本仓库根目录当前没有统一 LICENSE 文件。

## 开发与持续更新

```powershell
python modules\_framework\gen_catalog.py
python scripts\update_readme.py
python scripts\update_readme.py --check
python -m unittest tests.test_readme_update tests.test_2026_method_modules -v
python modules\_framework\qc_lint.py modules\12_genome_annotation_variants
git diff --check
```

README 的动态指标和领域表由 `scripts/update_readme.py` 维护；GitHub Actions 在 push、pull request 和每周计划任务中检查漂移。写作范式、触发条件和人工维护边界见 [`docs/README_MAINTENANCE.md`](docs/README_MAINTENANCE.md)。


## 贡献

提交新模块前，请先搜索目录和现有 Issue，核实论文、官方代码、许可、真实 CLI/API 与输出格式；随后按测试→实现→目录生成→质量检查的顺序提交。问题与建议可通过 [GitHub Issues](https://github.com/fsy2004/bioinfo-reusable-code/issues) 反馈。
