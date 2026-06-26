# 21_疾病负担_共病_GBD_NHANES_CHARLS

本模块用于补齐疾病负担与共病数据库分析能力，优先服务 GBD、NHANES、CHARLS 三类公开/半公开数据库的共病研究。每个子目录现各含一个**接地于真实工具、开箱即跑的 turnkey 模块**(合成数据 + 出图 + README)。

## 模块

| 子目录 | 模块 | 用途 |
|--------|------|------|
| `01_GBD` | [527 GBD burden trend](01_GBD/527_gbd_burden_trend/) | ASR 趋势 + EAPC + Das Gupta 分解 + SDI |
| `02_NHANES` | [528 NHANES survey-weighted](02_NHANES/528_nhanes_survey_weighted/) | svydesign → svymean/svyby/svyglm 加权估计 |
| `03_CHARLS` | [529 CHARLS longitudinal](03_CHARLS/529_charls_longitudinal_cohort/) | 波次描述 + 等百分位等值化 + LMM + Cox/KM |
| `04_comorbidity_network` | [530 comorbidity network](04_comorbidity_network/530_comorbidity_network/) | 2×2 关联 → igraph → Louvain 社区 + hub |

这 4 个模块均**接地于 `99_external_sources/` 下已克隆的真实工具仓库**(R-script-for-GBD/GBD2021、nhanes、charls_memory_equating、comorbidity_networks),用本机已装包实现,外部/重型步骤(Joinpoint/BAPC/equate)在各模块 README 诚实标注替代方案。`99_external_sources/` 与选题/文献草稿仍本地保留(不入库)。

## 护栏

1. 选题前必须查宽：题名核心变量 + 数据库 + 方法组合，包含 PubMed、bioRxiv、Google Scholar 和 GitHub。
2. 通用 pipeline 不是卖点；卖点应来自具体疾病问题、数据互证和真实证据锚。
3. NHANES 机器学习必须低维、流行病学约束、报告 calibration 与 decision curve，并使用真外部验证。
4. CHARLS 是纵向证据锚，但数据获取、变量 harmonization 和失访处理要先确认。
5. GBD 是疾病负担锚，不足以单独证明个体共病因果；需要 NHANES/CHARLS/MR 或机制数据补强。
