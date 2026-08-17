# Reusable Bioinformatics Code Library

A domain-organized library of reusable R and Python bioinformatics
analyses. Pick the question you want to answer, open the matching
module, and follow its README for the exact input, method,
dependencies, and output figures.

Module READMEs are written in Chinese unless noted otherwise; the
[module catalog](modules/CATALOG.md) gives a full English-indexed
overview.

<p>
<img alt="modules" src="https://img.shields.io/badge/modules-196-blue">
<img alt="domains" src="https://img.shields.io/badge/domains-10-informational">
<img alt="languages" src="https://img.shields.io/badge/R%204.4%20%7C%20Python%203.12-informational">
<img alt="figures" src="https://img.shields.io/badge/figures-vector%20PDF%20%2B%20300dpi%20PNG-success">
</p>

## What's inside

196 modules organized into 10 domains and 53 subcategories. Each row
below links to the domain directory; the complete per-module index is
[`modules/CATALOG.md`](modules/CATALOG.md).

| Domain | Main analyses | Directory |
|---|---|---|
| 01 Single-cell | loading & QC, batch integration, annotation & cell typing, composition, differential expression / pseudobulk, trajectory & RNA velocity, CNV / clonality, pathway & TF activity, sc-to-bulk phenotype association, foundation models | [`modules/01_single_cell`](modules/01_single_cell) |
| 02 Spatial transcriptomics | cell segmentation, spatial domains & SVG, spatial statistics, deconvolution & single-cell mapping, slice registration & 3D reconstruction, spatial cell-cell communication, spatial multi-omics, spatial foundation models | [`modules/02_spatial_transcriptomics`](modules/02_spatial_transcriptomics) |
| 03 Virtual perturbation | virtual knockout, GRN inference, causal representation & counterfactuals, drug perturbation & response, perturbation prediction benchmarks | [`modules/03_virtual_perturbation`](modules/03_virtual_perturbation) |
| 04 Causal inference & genetics | instrument preparation, two-sample MR, cis-MR & drug targets, mediation & MVMR, colocalization, TWAS / sc-eQTL, robust MR estimators | [`modules/04_causal_inference_genetics`](modules/04_causal_inference_genetics) |
| 05 Machine learning | feature selection, classification, survival ML, model explanation, uncertainty quantification, generalization & external validation | [`modules/05_machine_learning`](modules/05_machine_learning) |
| 06 Bulk omics | differential expression, GO/KEGG & GSEA, WGCNA / co-expression networks, multi-omics integration & subtyping, mutation / methylation / proteomics / metabolomics | [`modules/06_bulk_omics`](modules/06_bulk_omics) |
| 07 Clinical & translational | diagnostic models, prognosis & survival, immune infiltration & deconvolution, pharmacovigilance, disease burden & population cohorts | [`modules/07_clinical_translational`](modules/07_clinical_translational) |
| 08 Structure & drug design | molecular docking, molecular dynamics, virtual screening & scoring | [`modules/08_structure_drug_design`](modules/08_structure_drug_design) |
| 09 Network pharmacology | target database extraction, target intersection & set relations, druggability scoring | [`modules/09_network_pharmacology`](modules/09_network_pharmacology) |
| 10 Visualization | advanced scientific plot types, publication-grade figure templates, external plotting resources | [`modules/10_visualization`](modules/10_visualization) |

### Start from your research question

- "Which cells, states, or subpopulations change?" → start at **01 Single-cell** or **02 Spatial transcriptomics**.
- "Which genes, pathways, or regulators change?" → differential expression, enrichment, and activity analyses in **01 Single-cell** and **06 Bulk omics**.
- "How do cells communicate, and do spatial niches form?" → spatial statistics, cell-cell communication, and spatial multi-omics in **02 Spatial transcriptomics**.
- "Could a gene / pathway / drug target be causal?" → **03 Virtual perturbation**, **04 Causal inference & genetics**, **09 Network pharmacology**.
- "Build a diagnostic, prognostic, or risk model" → **05 Machine learning** and **07 Clinical & translational**.
- "How does a candidate small molecule bind, and is it stable?" → **08 Structure & drug design**.
- "Turn results into paper figures" → **10 Visualization**, or reuse each module's `assets/` and the shared plotting framework.

## MetaWingman skill (systematic review & meta-analysis)

[`MetaWingman`](https://github.com/fsy2004/MetaWingman) is an
end-to-end skill for systematic review and meta-analysis. Skills are
portable across LLM hosts, and MetaWingman is no exception: it runs on
your existing host model and requires no extra model API. It is not a
"meta-analysis app" that only wraps statistical commands — it organizes
methodological decisions, an auditable project state, literature search
and lawful retrieval, statistical code, writing conventions, and an AI
review loop into one workflow. A Gitee mirror lives at
[`fsy2004/meta-wingman`](https://gitee.com/fsy2004/meta-wingman).

### Stage coverage

| Stage | What the skill does | Human / audit checkpoints that remain |
|---|---|---|
| Topic & feasibility | frames the question, estimand, novelty and combinability; live-checks existing reviews, registrations, key studies | research question, scope, go / no-go decision |
| Protocol & registration | protocol structure, eligibility criteria, search plan, pre-specified analysis and bias control | protocol freeze, justification and approval of major deviations |
| Search & retrieval | auditable PubMed / Europe PMC / ClinicalTrials.gov search; Crossref / DOI checks; Unpaywall and open full text; imports authorized database exports | commercial database logins, CAPTCHAs, institutional licenses, final search confirmation |
| Dedup & screening | deterministic dedup, dual title/abstract and full-text screening, conflict arbitration, per-record exclusion reasons, PRISMA counts | independent judgment, conflict arbitration, full-text conclusions that cannot be inferred automatically |
| Extraction & appraisal | `record_id -> report_id -> study/trial_id -> result_id` lineage, multi-report / multi-timepoint / multi-arm relations, design-matched tools (RoB 2, ROBINS-I, QUADAS-2, ...) | source-locating checks, double-checking key numbers, final domain judgment |
| Quantitative / narrative synthesis | pairwise meta, network meta, diagnostic, prognostic, prevalence/proportions, dose-response, Bayesian, multilevel/RVE, IPD, sequential analysis, sensitivity & influence; SWiM when pooling is inappropriate | model and effect-measure choice, clinical-heterogeneity interpretation, the pooling decision |
| Certainty, writing & review | GRADE, absolute effects, PRISMA reporting, unified numbers / terminology / figure style, AI multi-perspective review–revision–recheck loop, living-review updates | strength of conclusions, author responsibility, pre-submission item-by-item verification |

### Safeguards

- Live verification of time-sensitive facts, guidelines,
  registrations, retractions, and references before anything is
  finalized.
- No fabricated database access, hit counts, PDFs, screening
  decisions, extracted values, analysis outputs, confidence intervals,
  GRADE ratings, or completion status; missing data stays missing and
  is queued for follow-up.
- Full texts come only from public APIs, open-access links,
  user-provided files, or authorized institutional routes.
- API keys are read from environment variables or user-approved secret
  stores, never written to project files or git history.
- Raw responses, queries, timestamps, file hashes, and decision logs
  are kept at every step.

### What it ships & how to invoke

The MetaWingman repository includes 12 methodological references, 8
Python workflow scripts, 61 R analysis manifests, 15 task-oriented R
adapters (plus one shared adapter), and 26 lower-level R tool modules.
All 61 manifests pass migration checks with their declared example
inputs; these checks prove the interfaces run, not that any specific
project's data, methods, or interpretation are correct.

```text
$metawingman

Research question: ...
Current stage: topic / protocol / search / screening / extraction /
appraisal / analysis / writing / review / update
Existing materials: protocol, search string, RIS/CSV, PDFs, extraction
tables, or analysis data (if any)
Expected output: decision records, reproducible project, figures,
GRADE tables, manuscript, or review report
```

The skill first identifies the review type and current stage, then
loads only the relevant methodology notes and tools. Statistical
analyses follow the protocol, estimand, study design, effect measure,
dependency structure, and information available — having every tool is
never a reason to run every analysis.

## How a module works

Each module is a self-contained analysis directory:

```text
<domain>/<subcategory>/<NNN_module>/
├── <NNN_module>.R | .py     # main script
├── README.md                # input, method, dependencies, output, limitations
├── example_data/            # small example or synthetic data (if applicable)
└── assets/                  # rendered PNG/PDF preview figures (if applicable)
```

Recommended order:

1. Read the module README first: input objects, sample/cell level,
   dependencies, scope.
2. Run the example in `example_data/`, then point `--input` /
   `--outdir` at your own data.
3. For foundation models, server methods, or external toolchains, also
   check the status marks and
   [`SERVER_DEPENDENCIES.md`](modules/_framework/SERVER_DEPENDENCIES.md).

## Quick start

```bash
git clone https://github.com/fsy2004/bioinfo-reusable-code.git
cd bioinfo-reusable-code

# enter the domain directory, then run per the module README
cd modules
Rscript 06_bulk_omics/01_differential_expression/010_geo_deg_volcano_heatmap_pca/010_*.R
python 02_spatial_transcriptomics/02_domains_svg_stats/543_squidpy_spatial_statistics/543_*.py
```

Exact commands and arguments are per module README. The default target
environments are **R 4.4** and **Python 3.12**; individual modules may
need local, server, or GPU environments.

## Status legend

| Mark | Meaning |
|---|---|
| ✅ | bundled/synthetic example runs locally as-is |
| 🟡 | real baseline or core pipeline runs locally; the full method needs server dependencies |
| 🔴 | heavy, GPU, or external toolchain; guarded wrapper or reference implementation |
| 📄 | template or upstream script; prepare your own data and dependencies |
| 📦 | vendored/external tool list or local source, not a complete turnkey module |
| 🗃️ | kept locally only, git-ignored, never published |

"Runs" does not mean "fits all data". Check each module README for
input limits, sample level, baselines, random seeds, and output
interpretation.

## Shared framework

[`modules/_framework/`](modules/_framework/) provides cross-module
conventions and tooling:

- `CATALOG.md` — full domain → subcategory module index: purpose,
  input → output, dependencies, language, figures, status.
- `TOOL_SELECTION_GUIDE.md` — pick methods and modules by analysis task.
- `CONVENTIONS.md` — directory layout, reuse-first, paths, random
  seeds, scripting conventions.
- `theme_pub.R` / `pubstyle.py` — unified palettes, themes, PDF + 300
  dpi PNG export.
- `ANALYSIS_TEMPLATE/` — scaffold for new projects with config,
  checkpoints, and environment snapshots.
- `QUALITY_CHECKLIST.md` / `qc_lint.py` — analysis quality checklist
  and static checks.
- `SERVER_DEPENDENCIES.md` — server / GPU / external toolchain notes.
- `BENCHMARKS_AND_CRITIQUES.md` — comparison scales and common pitfalls
  for complex models, foundation models, and perturbation prediction.

## Reproducibility

- Reuse real tools and existing modules; never rewrite analysis APIs
  from memory.
- Use relative paths, fixed random seeds, and per-module environment
  notes; no hardcoded machine paths.
- Outputs default to `results/`; local run outputs, caches, and
  upstream sources are not published as repository content.
- Figures prefer vector PDF plus 300 dpi PNG; complex models should
  also report transparent, comparable baselines.
- Conclusions follow the module's actual outputs and data level;
  example-data results are never presented as real findings.

## License

Each module follows the license of the tools and methods it uses.
Vendored third-party code keeps its original license; see the relevant
module README and upstream repository.

---

## 中文说明

# 可复用生物信息学代码库

按分析问题组织的 R / Python 生信模块库：先按「我要做什么分析」定位域和
子类，再打开具体模块的 README，查看输入格式、分析方法、依赖与输出图。

模块 README 以中文为主（少数为英文）；[模块索引](modules/CATALOG.md)提供
全量逐模块总览。

<p>
<img alt="modules" src="https://img.shields.io/badge/modules-196-blue">
<img alt="domains" src="https://img.shields.io/badge/domains-10-informational">
<img alt="languages" src="https://img.shields.io/badge/R%204.4%20%7C%20Python%203.12-informational">
<img alt="figures" src="https://img.shields.io/badge/figures-vector%20PDF%20%2B%20300dpi%20PNG-success">
</p>

## 这套库能做什么

仓库包含 **196 个模块**，按 **10 个分析域、53 个子类**组织。下面每个域
链接进入该域目录；完整的逐模块索引见 [`modules/CATALOG.md`](modules/CATALOG.md)。

| 域 | 主要分析内容 | 进入目录 |
|---|---|---|
| 01 单细胞分析 | 数据读入与质控、批次整合、细胞注释与分型、组成差异、差异表达/pseudobulk、轨迹与 RNA velocity、CNV/克隆、通路与 TF 活性、单细胞与 bulk 表型关联、基础模型 | [`modules/01_single_cell`](modules/01_single_cell) |
| 02 空间转录组 | 细胞分割、空间域与 SVG、空间统计、解卷积与单细胞映射、切片配准与三维重建、空间细胞通讯、空间多组学、空间基础模型 | [`modules/02_spatial_transcriptomics`](modules/02_spatial_transcriptomics) |
| 03 虚拟扰动 | 虚拟敲除、基因调控网络推断、因果表示与反事实、药物扰动与响应、扰动预测基准 | [`modules/03_virtual_perturbation`](modules/03_virtual_perturbation) |
| 04 因果推断与遗传 | 工具变量准备、两样本 MR、cis-MR 与药靶、中介与 MVMR、共定位、TWAS/sc-eQTL、稳健 MR 估计量 | [`modules/04_causal_inference_genetics`](modules/04_causal_inference_genetics) |
| 05 机器学习 | 特征筛选、分类模型、生存机器学习、模型解释、不确定性量化、泛化与外部验证 | [`modules/05_machine_learning`](modules/05_machine_learning) |
| 06 Bulk 组学 | 差异表达、GO/KEGG 与 GSEA、WGCNA/共表达网络、多组学整合与分型、突变/甲基化/蛋白质组/代谢组 | [`modules/06_bulk_omics`](modules/06_bulk_omics) |
| 07 临床与转化 | 诊断模型、预后与生存、免疫浸润与解卷积、药物警戒、疾病负担与人群队列 | [`modules/07_clinical_translational`](modules/07_clinical_translational) |
| 08 结构与药物设计 | 分子对接、分子动力学、虚拟筛选与打分 | [`modules/08_structure_drug_design`](modules/08_structure_drug_design) |
| 09 网络药理学 | 靶点数据库提取、靶点交集与集合关系、成药性评分 | [`modules/09_network_pharmacology`](modules/09_network_pharmacology) |
| 10 可视化 | 高级科研图型、出版级图形模板、外部绘图资源 | [`modules/10_visualization`](modules/10_visualization) |

### 按研究问题找入口

- 想回答「哪些细胞、状态或亚群发生变化」：从 **01 单细胞** 或 **02 空间转录组**开始。
- 想回答「哪些基因、通路或调控因子发生变化」：查看 **01 单细胞**、**06 Bulk 组学** 的差异表达、富集与活性分析。
- 想回答「细胞之间如何通讯、空间上是否形成生态位」：查看 **02 空间转录组** 的空间统计、细胞通讯与空间多组学。
- 想回答「某个基因/通路/药物靶点是否可能具有因果作用」：查看 **03 虚拟扰动**、**04 因果推断与遗传**、**09 网络药理学**。
- 想建立「诊断、预后或风险预测模型」：查看 **05 机器学习** 与 **07 临床与转化**。
- 想回答「候选小分子如何结合、是否稳定、如何筛选」：查看 **08 结构与药物设计**。
- 想把结果整理成论文图：查看 **10 可视化**，或直接复用各分析模块的 `assets/` 与统一绘图框架。

## 系统综述与 Meta 分析 Skill（MetaWingman）

[`MetaWingman`](https://github.com/fsy2004/MetaWingman) 是系统综述与
Meta 分析全流程 skill。skill 天然跨 LLM 通用，MetaWingman 同样如此：跑在
你现有的宿主模型上，不需要额外模型 API。它不是只封装统计命令的「Meta
软件」，而是把方法学决策、可审计项目状态、文献检索与合法获取、统计代码、
写作规范和 AI 审稿闭环组织在一起。Gitee 镜像见
[`fsy2004/meta-wingman`](https://gitee.com/fsy2004/meta-wingman)。

### 能覆盖什么

| 阶段 | Skill 提供的辅助 | 必须保留的人工/审计关口 |
|---|---|---|
| 选题与可行性 | 明确问题框架、估计目标、创新性与可合并性，联网核查既有综述、注册与关键研究 | 研究问题、范围和继续/停止决定 |
| 方案与注册 | 协议结构、资格标准、检索计划、分析和偏倚控制预设 | 协议冻结、重要偏离的说明与批准 |
| 检索与文献获取 | PubMed、Europe PMC、ClinicalTrials.gov 的可审计检索；Crossref/DOI 核验；Unpaywall 与开放全文获取；导入授权数据库导出 | 商业数据库登录、验证码、机构许可和最终检索确认 |
| 去重与纳排 | 确定性去重、双人题录/摘要与全文筛选、冲突仲裁、逐篇排除理由、PRISMA 计数 | 独立判断、冲突仲裁和不可自动推断的全文结论 |
| 提取与质量评价 | `record_id -> report_id -> study/trial_id -> result_id` 谱系，多报告/多时间点/多臂关系，RoB 2、ROBINS-I、QUADAS-2 等设计匹配工具 | 原文定位复核、关键数值双查与最终领域判断 |
| 定量/叙述综合 | 双臂 Meta、网络 Meta、诊断、预后、患病率/比例、剂量反应、Bayesian、多层/RVE、IPD、序贯分析、敏感性与影响分析；不宜合并时转 SWiM | 模型与效应量选择、临床异质性解释、是否合并的决定 |
| 确定性、写作与审稿 | GRADE、绝对效应、PRISMA 报告、统一数字/术语/图表样式、AI 多视角审稿—修订—复核循环、living review 更新 | 结论强度、作者责任、投稿前逐项核验 |

### 行为边界

- 时效性事实、指南、注册、撤稿和参考文献在定稿前联网核验。
- 不虚构数据库访问、检索条数、PDF、筛选决定、提取值、分析输出、置信区间、
  GRADE 等级或完成状态；缺失数据保留为缺失并进入追索队列。
- 文献下载只走公开 API、开放获取链接、用户提供文件或用户已获授权的机构路径。
- API 密钥从环境变量或经用户批准的密钥存储读取，不写入项目文件或 Git 历史。
- 每一步保留原始响应、查询式、时间戳、文件哈希与决策日志。

### 内置方法与调用

MetaWingman 仓库包括 **12 份方法学参考、8 个 Python 工作流脚本、61 个 R 分析清单、15 个面向任务的 R 适配器（另有 1 个公共适配器）和 26 个底层 R 工具模块**。61 个清单均已用声明的示例输入完成迁移验证；这些验证证明接口能够运行，不代替具体课题的数据核查、方法选择和科学解释。

```text
$metawingman

研究问题：……
当前阶段：选题 / 协议 / 检索 / 纳排 / 提取 / 评价 / 分析 / 写作 / 审稿 / 更新
已有材料：协议、检索式、RIS/CSV、PDF、提取表或分析数据（如有）
期望输出：决策记录、可复现项目、图表、GRADE 表、稿件或审稿报告
```

Skill 会先识别综述类型和当前阶段，再只加载相关方法说明与工具；统计分析服从协议、估计目标、研究设计、效应量、依赖结构和信息量，不以「工具齐全」为理由堆叠不必要分析。

## 模块怎么用

每个模块是一个相对独立的分析目录，通常包含：

```text
<domain>/<subcategory>/<NNN_module>/
├── <NNN_module>.R | .py     # 主脚本
├── README.md                # 输入、方法、依赖、输出与限制
├── example_data/            # 小型示例或合成数据（如适用）
└── assets/                  # 已渲染的 PNG/PDF 预览图（如适用）
```

推荐顺序：

1. 先读模块 README，确认输入对象、样本/细胞层级、依赖和适用范围。
2. 先用 `example_data/` 跑通示例，再用 `--input` / `--outdir` 指向自己的数据。
3. 对基础模型、服务器方法或外部工具链，同时查看状态标记和
   [`SERVER_DEPENDENCIES.md`](modules/_framework/SERVER_DEPENDENCIES.md)。

## Quick start

```bash
git clone https://github.com/fsy2004/bioinfo-reusable-code.git
cd bioinfo-reusable-code

# 进入对应域的目录，再按模块 README 运行
cd modules
Rscript 06_bulk_omics/01_differential_expression/010_geo_deg_volcano_heatmap_pca/010_*.R
python 02_spatial_transcriptomics/02_domains_svg_stats/543_squidpy_spatial_statistics/543_*.py
```

模块的具体命令和参数以各自 README 为准。默认目标环境为 **R 4.4** 和
**Python 3.12**；不同模块的依赖可能需要本地、服务器或 GPU 环境。

## 状态标记

| 标记 | 含义 |
|---|---|
| ✅ | 自带/合成示例可本机运行，无需修改 |
| 🟡 | 真实基线或核心流程可本机运行，完整方法需安装服务器依赖 |
| 🔴 | 重型、GPU 或外部工具链；保留守卫式封装或参考实现 |
| 📄 | 模板或上游脚本；需要自行准备数据和依赖 |
| 📦 | Vendored/外部工具清单或本地来源，不是完整 turnkey 模块 |
| 🗃️ | 仅本地保留，由 Git 忽略，不会发布到仓库 |

「可运行」不等于「适合所有数据」。请优先检查每个模块 README 中的输入限制、样本层级、基线、随机种子和输出解释。

## 统一框架

[`modules/_framework/`](modules/_framework/) 提供跨模块的约定与工具：

- `CATALOG.md`：按域 → 子类的完整模块索引，包含用途、输入→输出、依赖、语言、图型与状态。
- `TOOL_SELECTION_GUIDE.md`：按分析任务选择方法和模块。
- `CONVENTIONS.md`：目录结构、复用优先、路径、随机种子和脚本约定。
- `theme_pub.R` / `pubstyle.py`：统一配色、主题和 PDF + 300 dpi PNG 导出。
- `ANALYSIS_TEMPLATE/`：带配置、检查点和环境快照的新项目脚手架。
- `QUALITY_CHECKLIST.md` / `qc_lint.py`：分析质量清单与静态检查。
- `SERVER_DEPENDENCIES.md`：服务器/GPU/外部工具链依赖说明。
- `BENCHMARKS_AND_CRITIQUES.md`：复杂模型、基础模型和扰动预测的比较尺度与常见陷阱。

## 可复现约定

- 优先复用真实工具和现有模块，不凭记忆重写分析 API。
- 使用相对路径、固定随机种子和模块级环境说明；避免硬编码本机路径。
- 运行输出默认写入 `results/`，本地运行结果、缓存和上游源码不作为仓库内容发布。
- 图形优先输出矢量 PDF 与 300 dpi PNG；复杂模型尽量同时报告透明、可比较的基线。
- 结论应以模块实际输出和数据层级为准，不把示例数据结果当作真实研究结论。

## License

每个模块遵循所用工具与方法的许可。第三方代码保留原始许可证，详见对应模块 README 与上游仓库。
