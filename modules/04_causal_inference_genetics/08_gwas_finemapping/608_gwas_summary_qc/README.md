# 608 · GWAS summary statistics 与 LD 输入审计

用途：在区域 fine-mapping、共定位和 cis-MR 前，检查 **稠密区域**汇总统计、语义元数据与有符号 LD 是否相容。该入口不运行 GWAS、clumping、liftover 或等位基因翻转。

## 输入与运行

`summary.tsv` 必含 `variant_id,chr,pos,ea,nea,beta,se,p,eaf`；可选 `n`。列名用制表符分隔。`beta` 必须是线性效应或 log(OR)，不能直接放 OR。`eaf` 对应 `ea`。不能先只保留显著 SNP 或独立工具变量再做共定位。

```powershell
python gwas_summary_qc.py --summary summary.tsv --metadata metadata.json --ld-matrix ld.txt --ld-variants ld_variants.tsv --out-dir audit
```

`ld.txt` 是无表头的空白分隔方阵，内容为 signed `r`；`r²` 不足以替代。`ld_variants.tsv` 必含 `variant_id,ea,nea`，行序与矩阵及汇总统计完全相同。可不提供 LD，只运行 summary 审计；后续 SuSiE 仍需要 LD。

`metadata.json` 需要 `trait_id,phenotype,genome_build,ancestry,beta_unit,effect_type,source,region,dense_region_confirmed`。build 为 `GRCh37` 或 `GRCh38`；effect_type 为 `linear_beta` 或 `log_odds`；dense_region_confirmed 必须明确为 `true`。使用 LD 时另需 `ld_genome_build,ld_ancestry,ld_source,ld_coded_allele="ea",ld_metric="signed_r"`。region 应写清 chr:start-end，不能由文件名推测。

## 输出与边界

`audit.json` 保存输入 SHA256、所有问题及其行号、LD 检查状态。退出码 2 表示错误；0 表示无错误，仍需读取 warnings。numpy 可用时做特征值 PSD 检查，缺失时明确标记未完成；不会静默安装。P 与 beta/SE 不一致仅报警，因为原 GWAS 可能使用非 Wald 检验。回文 SNP 报警且不自动删除；需要独立参考支持方向。
空文件、重复表头和不齐整行会明确拒绝。输出目录须不存在或为空，避免覆盖已有证据。

build、祖源和 dense-region 声明来自用户元数据，脚本不能证明其真实性。通过输入审计不代表工具有效、因果成立或分析已经完成。大区域 LD 读入和 PSD 检查成本分别约 O(p²) 与 O(p³)。

`region` 使用 1-based 闭区间，如 `chr1:100-200`；`chr` 支持 1–22、X、Y、MT/M 及可选 `chr` 前缀。逐行拒绝无效染色体或超出声明区间的变异；区间包含核对不能证明所有区域变异已完整提供。

## 来源

- [GWASLab](https://github.com/Cloufield/gwaslab)：完整的汇总统计处理、参考校验和可视化工具；本入口是独立原创的轻量审计，没有复制其实现。
- [coloc 数据约定](https://chr1swallace.github.io/coloc/articles/a02_data.html)：区域内 SNP 覆盖及标准输入结构。
- [SuSiE RSS](https://stephenslab.github.io/susieR/reference/susie_rss.html)：z-score 与相关矩阵输入。

证据状态：2026-10-19 核验上游文档与真实仓库；仅做合成输入行为验证，未运行任何真实课题数据。
