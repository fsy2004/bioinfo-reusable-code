# MR与共定位方法库维护入口

权威方法元数据在 [模块611](../../../modules/04_causal_inference_genetics/09_method_library/611_mr_coloc_library/README.md) 的 `methods.json`；按问题读一条方法及关联来源，避免扫描全部源码。机器索引为 [methods JSONL](../../../catalog/mr_coloc_methods.jsonl)、[sources JSONL](../../../catalog/mr_coloc_sources.jsonl)。本次固定2026-10-21来源审计，今后仍需核新版与新原始论文。

## 覆盖与状态

分别统计：43方法族概览、具体算法/模式/参数变体、多个软件实现、设计/诊断、共定位和纯fine-mapping。总记录数不能叫作独立MR算法数。重名算法保留一个稳定ID、合并软件实现与别名。无公开代码/正式输入未核的方法保留文献参考；不据摘要编造估计器。当前分类与源码状态以 `python scripts/mr_coloc.py check` 输出和机器索引为准。

已补入传统/稳健/弱工具MR、相关多效性、选择与重叠校正、多变量/多结局、cis区域、个体/非线性方法，以及ABF、SuSiE、条件共定位、多性状/多组织和多祖源共定位。高级方法包括MRMix、MR-APSS、MRAID、BWMR、LHC-MR、MR-CUE、PMR/moPMR、MR2/MrDAG、MR-Corr2、MR-LDP、MAPLE、cisMR-cML、METEOR等；新方法的出版/预印本状态、代码和运行证据分别标记。

本机接口证据见 [runtime inventory](runtime_inventory_20261021.json) 与 [真实原生行为检查](native_validation_20261021.md)。后者使用合成数据调用真实包，不是统计性能比较或课题科学分析。RAPS源码可获取，但本机依赖未安装，原生调用明确failed。源码取得、函数导出、实际返回、科学有效性是四种不同状态。

## 按问题择优

| 数据与风险 | 候选路线 | 决定前必须核对 |
|---|---|---|
| 仅1个合格IV | Wald ratio | 强度、功能与排除限制；不能做Egger/PRESSO作为无多效性证明 |
| 多个独立、较强IV | IVW主估计＋少数假设不同的稳健敏感性 | median有效权重、mode有效簇、Egger InSIDE/NOME；数量下限不代表功效 |
| 较多弱IV | RAPS、dIVW/pIVW等 | 各自平衡/独立多效性模型；相关适应证路径未必可修复 |
| 相关cis工具 | 广义IVW/Egger/ML、PC-GMM或专用cis方法 | signedLD方向/祖源/版本、数值稳定与弱条件工具 |
| 相关多效性或全基因组混杂 | CAUSE、MR-APSS、MRAID、LHC-MR、MAPLE等有条件候选 | 完整背景GWAS、LD/误差结构、模型识别；不能套少数lead SNP |
| 样本重叠/赢家诅咒 | 优先独立发现/估计；条件满足考虑MRlap、CARE/RIVW、SimSS | 选择与误差分布、重叠结构；框架不替代有效暴露构念 |
| 多暴露/中介 | MVMR及相应弱工具/稳健模型 | 条件F、暴露效应协方差、直接效应假设；不能只看单变量F |
| 双性状区域近似单信号 | coloc.abf＋预设先验网格 | 完整区域、sdY/病例比例、H0–H4与信号强度；不强迫最大H4 |
| 区域多信号 | coloc-SuSiE；PWCoCo/其他模型为有条件备选 | 每个性状的signedLD/精细定位质量，信号对而非全区一个H4 |
| 多性状/组织或跨祖源 | HyPrColoc/moloc、CAFEH/ColocBoost、MACA等 | 每法单/多信号、样本重叠与祖源LD假设；非直接MR效应 |

结果同向或诊断P>.05不能证明IV有效；不通过方法多数投票选结论。遗传使用倾向、内源激素、分子表达与实际药物剂量要分别定义。具体临床项目的结局、数据人数、主次分析与未发表假设保留在其独立研究仓库。

## 更新与许可

1. 从原始方法论文和作者/官方软件确认身份、输入、假设及真实entrypoint。登记DOI、固定commit、版本与许可证据；综述或包装器仅作检索线索。
2. 合并既有算法和别名，再新增稳定ID；源码已获取不自动升级adapter/runtime/scientific状态。输入未精审的参考方法不加入自动候选。
3. 获取固定源码；检查receipt及每文件SHA。完整快照与代码子集分别统计；无代码、过大、未取得、失败及许可缺口均保留。运行 `python scripts/snapshot_mr_sources.py` 保存可移植的有日期receipt；CI无本地缓存时据此核对索引，不能据历史快照声称本机已有源码。实际字节校验见 [source integrity](source_integrity_20261021.json)，分类计数见 [coverage summary](../../../catalog/mr_coloc_summary.json)。
4. 自主适配器调用原生API；关键输入拒绝与输出保真用行为检查，真实研究另做分析级验证。不得把本地旧概念实现当作者算法等价。
5. `python scripts/mr_coloc.py build`，然后按AGENTS更新全库索引/README；核最终diff、索引、原生receipt和暂存范围再提交。

第三方源码不进入本Git仓库。SharePro_coloc核得CC BY-NC-ND 4.0；moloc等有缺失许可，部分包DESCRIPTION占位或与仓库LICENSE冲突。记录原文而不替作者选择许可。私有阅读缓存不表示获得再分发、改编或任何用途授权；对外复用按各来源许可处理。

## 持续缺口

“所有MR方法”保留为持续覆盖目标。尚需跟踪未公开代码的新方法、单样本/非线性/家系/交互设计、部分多变量弱工具方法、非标准结局与跨祖源推断的具体软件版本。MRZ早期页面等不能仅凭摘要编程。未提供每一种方法的统一运行适配器，也未完成所有方法的统计性能benchmark；参考条目的选择须回到原文和作者软件。
