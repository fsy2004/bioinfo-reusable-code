# 2026 生物信息学方法扫描（冻结于 2026-08-20）

## 范围与判定

本清单覆盖 2026-01-01 至 2026-08-20，重点复核上一版清单冻结日 2026-06-26 之后的方法。检索入口包括 Nature Methods、Bioinformatics、Genome Biology、Nature Biotechnology 等期刊官网，以及论文链接的官方代码仓库。纳入模块需要：论文或正式方法说明、官方代码、可核实许可、稳定输入/输出契约、可在本库诚实实现的本地审计边界。没有明确许可或稳定运行契约的候选保留在观察区。

本文件是方法候选登记，不是系统综述；“已覆盖”仅指本库存在相应分析入口，不代表对该领域所有方法进行了穷尽性检索。

## 全领域覆盖矩阵

| 领域 | 本库状态 | 2026-08-20 决策 |
|---|---|---|
| 单细胞、空间转录组、虚拟扰动 | 已有独立域 | 保持现有模块；SpaMTP 进入观察区 |
| 遗传因果、QTL、群体遗传 | 已有 MR/QTL 主干 | ARGscape 进入观察区；不在许可未明确时复制代码 |
| 机器学习与基础模型 | 已有独立域 | 继续要求简单基线、外部验证与不确定性 |
| bulk RNA、转录组装 | 已有 bulk 域 | StringTie3 进入观察区，待形成稳定增量模块 |
| 临床转化、药物、结构、网络药理 | 已有独立域 | 不因“新”而重复已有能力 |
| 证据合成与统计可视化 | 已有独立域 | 保持 MetaWingman 与 597 工具包入口 |
| 基因组变异与注释 | 原目录缺少顶层域 | 新增 598 ClairS 与 599 EviAnn |
| 比较基因组与系统发育 | 原目录缺少顶层域 | 新增 600 OrthoFinder v3 |
| 微生物组与宏基因组 | 原目录缺少顶层域 | 新增 601 CoverM |
| 表观基因组与染色质 | 原目录缺少顶层域 | 新增 602 damidBind |
| 蛋白质组与代谢组 | 原仅混在 bulk 子类 | 新增 603 DiaReport，建立独立顶层域 |
| 多组学、免疫组学、临床组学 | 已在 bulk/临床/单细胞交叉覆盖 | 后续按稳定输出契约拆域，不建空壳 |
| 图像组学、病理、成像质谱 | 尚无稳定通用模块 | 登记为下一轮检索范围 |
| 宏生态、病毒组、真菌组 | 可由宏基因组域扩展 | 按真实项目需求增设子类 |

## 本轮落地方法

| 方法 | 论文/官方源 | 代码与许可 | 本地模块边界 |
|---|---|---|---|
| ClairS | doi:10.1038/s41592-026-03152-4 | https://github.com/HKU-BAL/ClairS · BSD-3-Clause | VCF 输出审计 + 已核实 CLI |
| EviAnn | doi:10.1038/s41592-026-03156-0 | https://github.com/alekseyzimin/EviAnn_release · GPL-3.0 | GFF3 注释审计 + 已核实 CLI |
| OrthoFinder v3 | https://orthofinder.github.io/OrthoFinder/ | https://github.com/OrthoFinder/OrthoFinder · GPL-3.0 | Orthogroups 表审计 + 已核实 CLI |
| CoverM | doi:10.1093/bioinformatics/btaf147 | https://github.com/wwood/CoverM · GPL-3.0 | 丰度长表审计 + 已核实 CLI |
| damidBind | doi:10.1093/bioinformatics/btag512 | https://github.com/marshall-lab/damidBind · GPL-3.0 | 差异结合导出表审计 |
| DiaReport | doi:10.1093/bioinformatics/btag435 | https://github.com/Gevaert-Lab/diareport · MIT | 差异蛋白导出表审计 |

## 观察区

| 候选 | 来源 | 暂缓原因 |
|---|---|---|
| SpaMTP | doi:10.1038/s41592-026-03140-8 · https://github.com/GenomicsMachineLearning/SpaMTP | 官方仓库许可元数据未明确；依赖模型与空间数据，先不复制代码 |
| ARGscape | doi:10.1093/bioinformatics/btag420 · https://github.com/chris-a-talbot/argscape | 官方仓库许可元数据未明确，先保留方法入口 |
| StringTie3 | https://github.com/gpertea/stringtie | 已有转录组装能力且本轮重点是缺失域；下次按真实长读长项目需求增补 |

## 持续更新规则

1. 每季度或真实项目提出新数据类型时，更新冻结日期和检索日志。
2. 新候选先核实论文、官方仓库、许可、实际 CLI/API 和输出文件，再决定模块化。
3. 先写失败测试，后写输出审计；重型上游只生成命令或守卫式调用。
4. 更新 `new_modules.json` 与 `gen_catalog.py` 后重新生成目录，运行 `qc_lint.py`、模块测试与 `git diff --check`。
5. 未通过上述闸门的条目留在观察区，不以空目录占位。
