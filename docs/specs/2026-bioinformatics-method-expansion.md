# 2026 生物信息学方法增补规范

日期冻结：2026-08-20

## 目标

在既有 11 个领域、197 个模块基础上，补齐当前目录中缺失的常规生物信息学领域，并把 2026 年新方法转成可审计、可复用、默认可运行的本地模块。此次增补不是声称在 Windows 工作站完整复现上游重型流程；默认路径只读取标准格式的代表性输出，完成输入校验、核心质量统计、图形与机器可读摘要。上游启动命令仅在官方 README 已核实参数时生成，不自动执行。

## 证据与时间范围

- 全量检索范围：2026-01-01 至 2026-08-20；增量重点为旧清单冻结日 2026-06-26 之后。
- 只把论文官网、DOI 页面、官方仓库与官方文档作为“已核实”来源。
- 代码许可不明确、上游 API 未核实、必须 GPU/大模型或尚无稳定输出契约的候选只进入冻结登记表，不落成可运行模块。
- 本次优先落地六个有公开论文、官方代码与稳定文件格式的工具：ClairS、EviAnn、OrthoFinder v3、CoverM、damidBind、DiaReport。

## 新领域与模块

| 领域 | 编号 | 模块 | 默认可运行边界 |
|---|---:|---|---|
| 基因组注释与变异 | 598 | ClairS long-read somatic audit | 读取 VCF，按过滤状态和变异类型汇总；生成已核实启动命令 |
| 基因组注释与变异 | 599 | EviAnn annotation audit | 读取 GFF3，检查层级、证据与完整性；生成已核实启动命令 |
| 比较基因组与系统发育 | 600 | OrthoFinder v3 orthogroup audit | 读取 Orthogroups.tsv，统计物种覆盖、直系同源组大小与单拷贝候选 |
| 微生物组与宏基因组 | 601 | CoverM abundance audit | 读取 genome abundance 长表，校验非负值、样本和相对丰度和 |
| 表观基因组与染色质 | 602 | damidBind result audit | 读取导出的差异结合结果，按 FDR/效应方向审计并绘图 |
| 蛋白质组与代谢组 | 603 | DiaReport result audit | 读取 DIA/MSqRob 差异蛋白表，按 FDR/效应方向审计并绘图 |

## 统一输入输出契约

每个模块包含主 Python 脚本、五段式 README、`example_data/` 合成数据、`assets/` 展示图和被忽略的 `results/`。默认命令为：

```powershell
python <NNN_module>.py
```

通用参数：`--input`、`--outdir`；598、599、600、601 额外提供 `--print-upstream-command`，只打印官方命令，不执行。成功运行必须产生：

- `<id>_summary.csv`
- `<id>_summary.json`
- 至少一张 PNG 和一张 PDF

错误输入（缺列、负值、无法解析、空表）必须返回非零退出码并把可操作错误写入 stderr。

## 科学与工程边界

- 示例数据只证明 I/O、审计计算和出图路径，不证明上游方法的科学性能。
- 不根据合成数据报告生物学结论。
- damidBind 与 DiaReport 模块不重写其统计模型；只审计上游导出结果，并在 README 指向官方分析入口。
- 每项统计在测试中使用手工推导的固定期望值；不以源代码文本搜索代替行为测试。
- 本次不安装上游工具，不下载模型或大型参考数据，不修改远端 GitHub 状态。

## 主要来源

- ClairS: https://doi.org/10.1038/s41592-026-03152-4 ; https://github.com/HKU-BAL/ClairS
- EviAnn: https://doi.org/10.1038/s41592-026-03156-0 ; https://github.com/alekseyzimin/EviAnn_release
- OrthoFinder: https://github.com/OrthoFinder/OrthoFinder ; https://orthofinder.github.io/OrthoFinder/
- CoverM: https://doi.org/10.1093/bioinformatics/btaf147 ; https://github.com/wwood/CoverM
- damidBind: https://doi.org/10.1093/bioinformatics/btag512 ; https://github.com/marshall-lab/damidBind
- DiaReport: https://doi.org/10.1093/bioinformatics/btag435 ; https://github.com/Gevaert-Lab/diareport
