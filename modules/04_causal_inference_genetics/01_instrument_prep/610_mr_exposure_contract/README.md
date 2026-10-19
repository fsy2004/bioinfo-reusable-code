# 610 · MR 暴露构念与工具元数据契约

用途：方法准备时审计暴露到底代表什么，避免把“开药倾向”“内源性 cortisol”“NR3C1 表达”写成“糖皮质激素滥用”。仅计算每 SNP 近似 marginal F=`(beta/se)^2`，不拟合 MR 或生成条件 F。

```powershell
python mr_exposure_contract.py --contract contract.json --instruments instruments.tsv --out contract_audit.json
```

输入 `instruments.tsv`：`exposure_id,variant_id,beta,se`（制表符）。`contract.json`：顶层需要 `research_question,intended_estimand,design,outcome,sample_overlap,exposures`。design 为 `UVMR/MVMR/cis_MR`。每个 exposure 有 `id,phenotype,proxy_kind,source,genome_build,ancestry,beta_unit,n`。outcome 有 `id,phenotype,source,genome_build,ancestry,effect_scale,n`；若 effect_scale=`log_odds` 还需 `n_cases,n_controls`。

## 暴露种类

| proxy_kind | 代表的构念与限制 |
|---|---|
| prescribed_use_propensity | 遗传预测用药倾向；含适应证、就医、开药通路，不能直接解释为剂量效应或滥用 |
| endogenous_cortisol | 内源性 cortisol 水平；不能代表外源激素的剂量、途径、疗程 |
| gene_expression | 特定组织/细胞中目标基因表达；NR3C1 表达不自动等同受体活性、药理激动或药物滥用 |
| protein_abundance | 蛋白丰度；仍需方向、靶点特异性、平台和组织的验证 |
| misuse_or_abuse | 需有明确测量的滥用定义及验证；仅相关代理不能填成这一类 |

分子种类另需 `target_gene,tissue_or_cell_type,cis_window,molecular_qtl_type`。misuse_or_abuse 另需 `measured_definition` 和明确 `validated_misuse_phenotype=true`；脚本仅审计声明，并不能代替原始 phenotype 核验。顶层若 `claim="glucocorticoid_abuse_causal_effect"`，其他四类代理均以错误返回，提示重新界定问题。

`sample_overlap.status` 取 `unknown/none_documented/overlap_documented`；有文档状态须 `evidence` 来源，未知不等于不重叠。可选 fraction 取 [0,1]。祖源不同会报警；build 不同会报错，要求先显式 harmonization。

MVMR 需要 `conditional_strength={method,covariance_source,F_by_exposure:{exposure_id:value}}`，必须是其他真实 MVMR 工具已计算的值。未提供时输出 `UNAVAILABLE` 和 null，绝不以 marginal F 代替。已提供值标为 `USER_SUPPLIED_UNVERIFIED`。conditional F 低于 10 只作常见弱工具提示；F>10 不能证明工具有效。

输出 `contract_audit.json` 及退出码：2 有错误；0 无错误但仍可 `REVIEW_REQUIRED`。通过仅代表元数据契约通过，不能说明因果成立。正式分析还需工具独立性、LD、等位基因、样本重叠、多效性、共定位、多重检验及效应尺度审查。
既有输出路径会明确拒绝；复跑应使用新路径。病例数、对照数和总样本量均须有限且为正数。

## 来源与当前证据

- [2026 glucocorticoid-related genes/NSCLC](https://doi.org/10.1038/s41698-026-01421-1)：cis-eQTL 分子暴露实例，不能转写为用药滥用。
- [cortisol/cancer MR (2021)](https://doi.org/10.1038/s41416-021-01505-8)：内源晨间 cortisol 暴露，保留为历史基础文献。
- [MVMR 官方实现](https://github.com/WSpiller/MVMR)：条件工具强度需要联合资料与协方差假设。
- [2026 MR assumptions perspective](https://doi.org/10.1038/s41588-026-02546-6)：对 MR 假设与证据整合的边界。

2026-10-19：完成构念不等价、marginal F、样本重叠及 MVMR 缺失条件强度的合成行为验证；未收集真实课题数据、未执行科学 MR。
