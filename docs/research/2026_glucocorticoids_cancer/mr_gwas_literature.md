# 2026 GWAS/MR 方法准备：糖皮质激素与肿瘤方向

检索日期：2026-10-19。范围：2026 年正式研究/方法与设计指导，另列预印本和历史基础文献。检索是定向方法搜集，不是系统综述，也不声称穷尽全年论文或所有分析方法。逐项记录见 `mr_gwas_papers.json`；仓库与 commit/许可核验见 `mr_gwas_sources.json`。当前未收集课题数据或选择最终数据库版本。

## 直接相关与可迁移论文

| 日期 | 文献与证据读取层级 | 可学习方法 | 必须保留的边界 |
|---|---|---|---|
| 2026-04-18 | Hong et al. [Glucocorticoid signaling / NSCLC](https://doi.org/10.1038/s41698-026-01421-1)，已读出版商全文 | cis-eQTL MR、SMR/HEIDI、coloc、独立 outcome 复制和组织语境 | 暴露是广泛 GRG 表达，不能解释为外源激素用量或滥用；正文未找到已验证作者代码库 |
| 2026-01-14 | Went et al. [B-cell malignancy proteome MR](https://doi.org/10.1016/j.bneo.2026.100191)，PubMed 摘要 | 蛋白组 MR 与共定位/细胞/药靶证据整合 | 精确估计器和代码尚未从全文核实 |
| 2026-05-01 | Zhi et al. [Gastric cancer multi-ancestry targets](https://doi.org/10.3389/fonc.2026.1821512)，出版商方法片段与 PubMed 元数据 | 祖源分层 pQTL MR、SMR/HEIDI、coloc | 各祖源功效和 instrument coverage 不同；未验证作者代码库 |
| 2026-05-08 | Forde et al. [MR-SimSS](https://doi.org/10.1371/journal.pgen.1011949)，全文 | summary 层面模拟拆分；winner's curse、overlap 与弱工具敏感性 | 需足够工具和方法假设；不修复暴露构念错误；[真实作者包](https://github.com/amandaforde/mr.simss) |
| 2026-02-03 | Zou et al. [mvSuSiE](https://doi.org/10.1038/s41588-025-02486-7)，出版商 PDF 方法/代码声明 | 多性状联合 fine-mapping、PIP、可信集 | 多性状与多祖源设计不同；[作者包](https://github.com/stephenslab/mvsusieR) |
| 2026-03-30 | Wu et al. [GWFM](https://doi.org/10.1038/s41588-026-02549-3)，全文，另有 5 月勘误 | SBayesRC/GCTB 全基因组 fine-mapping、局部/全局可信集 | 需完整 summary 与匹配 LD、MCMC 和较高计算成本；[GCTB](https://github.com/jianzeng/GCTB) |
| 2026-07-21 | Shen et al. [Multi-ancestry colocalization](https://doi.org/10.1371/journal.pgen.1012221)，全文 | SuSiEx/MsCAVIAR 与 coloc/eCAVIAR 组合 | 需祖源分层 summary 与各自 LD；[MACA](https://github.com/cthyshn/MACA) 许可仍是占位文本，只登记参考 |
| 2026-05-08 | Reales et al. [eQTL-GWAS coloc design](https://doi.org/10.1371/journal.pgen.1012141)，全文 | shared-variant 覆盖审计、prior 敏感性、Bayesian FDR、细胞类型及共调控 | 多基因 coloc 不能单独选出功能介导基因；[harmonization code](https://github.com/GRealesM/GWAS_tools) |
| 2026-04-15 | Sanderson et al. [MR future directions](https://doi.org/10.1038/s41588-026-02546-6)，摘要和参考表 | 明确 estimand、实证评估 MR 假设、triangulation | Perspective，不是新估计器或药物安全结果；另有 4 月勘误 |
| 2026-10-15 | Zhang et al. [MRZ](https://doi.org/10.1038/s41467-026-78442-8)，accepted early version 摘要 | IV-outcome z mixture 区分有效/相关/非相关多效性工具 | 只读到摘要；未找到经核验代码，列方法观察清单，不据摘要复刻 |

另列 [2026-05-04 cross-ancestry breast-cancer protein MR 预印本](https://doi.org/10.64898/2026.05.04.26352055)，用于独立蛋白组、亚型与跨祖源复制设计参考，明确未经同行评审。[2021 cortisol/cancer MR](https://doi.org/10.1038/s41416-021-01505-8) 是历史背景，不混入 2026 新文献计数。

## 暴露构念：三个遗传路线不能互换

1. **内源性 cortisol**：遗传预测晨间皮质醇或相关生物标志物变化，估计终身生物学差异。它不能编码外源制剂、剂量、给药途径或超适应证使用。
2. **用药倾向 GWAS**：关联可来自疾病适应证、就医、开药及其他行为。需查原始 phenotype；“使用过/当前使用”不能自动解释为药效，更不能解释为滥用。
3. **NR3C1/cis-eQTL 或 pQTL 药靶 MR**：遗传代理特定组织中的表达或蛋白丰度；需目标、组织、方向及 on-target 验证。表达增加不必然等于受体激活，终身微小扰动不能直接换算临床剂量。

因此当前适合建设可替换暴露的基础方法，而将“滥用”因果主张保留到有可验证的剂量/疗程/适应证测量或明确 misuse GWAS 时再判定。观察性用药问题与遗传分子问题可构成证据互补，不能假定它们估计同一个效应。

## GitHub 方法资源与去重路线

| 方法 | 核验的来源 | 在现有库中的处理 |
|---|---|---|
| 基础 two-sample MR、harmonization、IVW/Egger/median、LOO | [TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR) | 复用 032/519 等现成模块，未在本轮重验旧科学结果 |
| OpenGWAS metadata/query/local LD 接口 | [ieugwasr](https://github.com/MRCIEU/ieugwasr) | 复用 instrument preparation；本轮不下载数据、不假定账户权限 |
| MVMR 与条件强度 | [MVMR](https://github.com/WSpiller/MVMR) | 复用 079/534；610 只审计条件强度是否有真实来源 |
| global/outlier/distortion pleiotropy tests | [MR-PRESSO](https://github.com/rondolab/MR-PRESSO) | 作为现成 MR 敏感性依赖，不再复制主流程 |
| correlated-IV/其他 MR 估计器 | [MendelianRandomization](https://github.com/cran/MendelianRandomization)（作者网站指向 CRAN 镜像） | 共用估计器依赖；不是所有稳健 MR 都适用于很少 SNP 的 cis 区域 |
| 弱工具稳健 MR | [mr.raps](https://github.com/qingyuanzhao/mr.raps) | 登记 canonical 作者仓，和 MR-SimSS 联合路线分别标记 |
| coloc / SuSiE fine-mapping | [coloc](https://github.com/chr1swallace/coloc)、[susieR](https://github.com/stephenslab/susieR) | 新 609：稠密区域真包 API、signed LD 顺序守卫、H0-H4/prior grid |
| 多性状 fine-mapping | [mvsusieR](https://github.com/stephenslab/mvsusieR) | 先登记高级扩展，避免与单性状 609 混淆 |
| 全基因组 fine-mapping | [GCTB](https://github.com/jianzeng/GCTB)、[SBayesRC](https://github.com/zhilizheng/SBayesRC) | 软件角色不同；2026 GWFM 正文指向 GCTB，不把 R 版 PRS 教程误当 GWFM 全部实现 |
| 跨祖源 fine-mapping/coloc | [SuSiEx](https://github.com/getian107/SuSiEx)、[MsCAVIAR](https://github.com/nlapier2/MsCAVIAR)、[MACA](https://github.com/cthyshn/MACA)、[CAVIAR/eCAVIAR](https://github.com/fhormoz/caviar) | 独立扩展；MACA 是论文链接代码，许可未明，不 vendor |
| GWAS summary 处理与 QC | [GWASLab](https://github.com/Cloufield/gwaslab)、[GWAS_tools](https://github.com/GRealesM/GWAS_tools) | 新 608 只提供轻量输入审计；完整规范化/参考校验仍由正式工具执行 |
| GWAS 工具衔接 | [gwasglue](https://github.com/MRCIEU/gwasglue) | README 标实验/开发中且部分方法 TODO；保留可选桥接，不能承诺全方法已实现 |
| LD/GWAS 基础与 PGS scoring | [PLINK2](https://github.com/chrchang/plink-ng)、[--score 官方说明](https://www.cog-genomics.org/plink/2.0/score) | scorer 资源核验，不代表 2026 breast 预印本作者的确切评分实现已核验 |
| LDSC 遗传相关 | [经典 bulik/ldsc](https://github.com/bulik/ldsc)、[Python3 CBIIT/ldsc](https://github.com/CBIIT/ldsc) | 经典 README 在 2026-01-16 指向新版；祖源匹配 LD scores；遗传相关不等于因果 |
| 独立队列效应 meta-analysis | [meta](https://github.com/guido-s/meta)、[metafor](https://github.com/wviechtb/metafor) | 参考实现；先核实样本独立与尺度，不把重叠联盟当独立复制 |

许可保留 API 原字段与 DESCRIPTION 实际声明两个层级：TwoSampleMR/ieugwasr 的 API `NOASSERTION` 不等于“无许可”，其 DESCRIPTION 声明 `MIT + file LICENSE`；MVMR/MR-PRESSO 声明 `GPL-3`；MendelianRandomization 声明 `GPL-2 | GPL-3`；coloc 声明 `GPL`；susieR/mvsusieR 声明 `BSD_3_clause + file LICENSE`。不擅自把不精确的 GPL 声明映射为某个 SPDX 版本。MACA 和 GWAS_tools 未建立有效源码许可，保持 reference-only。

GitHub 公共 API 触发共享未认证 rate limit 后，改用公开 `git ls-remote` 与 commit-pinned raw 文件；不读取或暴露凭据。个别仓库的 commit 未核得，明确留 null，不能把“链接可访问”写成“复现通过”。第三方 README 全文是工作证据，不应整体公开到用户库；公开索引只保留来源链接、短摘要、许可与 commit。

## 优先补齐的 3–5 个模块

1. **608 gwas_summary_qc**：build/祖源/单位、重复和数值、signed LD SNP/allele 顺序审计。不会假装执行 clumping 或全 GWAS。
2. **609 susie_coloc**：真 SuSiE RSS、多信号 coloc、prior sensitivity、PIP 与 H0–H4。拒绝非空输出目录，保留原始结果。
3. **610 mr_exposure_contract**：暴露构念、overlap 与 marginal F 审计；MVMR 条件 F 缺失明确为 null/UNAVAILABLE。
4. **MR-SimSS 适配器**：后续可围绕已公开真实包补 winner's curse/overlap 敏感性；不与现有 533 CARE/RIVW 简单重复。
5. **跨祖源 coloc 适配层**：只有祖源分层数据和对应 LD 齐备再实现；当前登记 MACA/SuSiEx 及许可边界。

当前前三项已在 workspace `work/mr_modules/` 准备，主代理负责入库。608/610 用合成数据验证输入与语义守卫；609 在本机真实 susieR 0.14.2、coloc 5.2.3 验证共享信号支持 H4、分离信号支持 H3、翻转 allele 拒绝及非空输出拒绝。它们是方法接口与错误处理验证，不是课题科学分析完成证据。
