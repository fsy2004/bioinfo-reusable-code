# 609 · SuSiE RSS fine-mapping + 多信号 coloc

用途：对两个已完成等位基因、build、区域和祖源核对的稠密区域汇总数据，调用真实 `susieR::susie_rss`、公开 `coloc::annotate_susie` 与 `coloc::coloc.susie`，输出信号对 H0–H4 及预先指定的 prior 敏感性。额外的 `coloc.abf` 是单因果变异假设下的对照。

```powershell
Rscript run_susie_coloc.R config.R
```

先复制并替换 `config.example.R` 的所有占位值。配置是普通受信任的本地 R 源文件，会执行其中代码。入口不会下载数据或安装包。
输出目录必须不存在或为空；非空目录会在写入前拒绝，复跑请指定新目录。

## 输入

两份 TSV 均含 `variant_id,ea,nea,beta,se,eaf`，可选常数列 `n`。它们必须保留同一稠密区域完全一致的 SNP 覆盖与顺序，且 ea/nea 方向相同；不得输入 clumped/显著 SNP 表。若需求交集，应在上游量化缺失信号覆盖，明确记录再生成输入，入口不静默删 SNP。

每个性状独立给出有符号 `r` LD 矩阵（无表头、空白分隔）和 `variant_id,ea,nea` 的矩阵行列顺序表。不同性状可使用各自匹配参考，单个祖源内分析；跨祖源任务使用 MACA/SuSiEx 等专用流程。文档里的 ancestry/build/LD 来源必须来自真实资料。当前入口要求常数 N；逐 SNP 样本量不一致会停下。

quant 类型须给出真实 `sdY`。病例对照类型给总 N 和样本病例比例 `s`，并明确允许 RSS 对 binary GWAS 的近似。不要把真实患病率代入 s。`estimate_residual_variance=FALSE` 适合外部 LD；入口不自动修复 LD 不匹配。

## 输出与解释

- `susie_fits.rds`、`variant_pip.tsv`：原始拟合与 SNP PIP；可信集保存在 fit 对象内。
- `susie_coloc_h0_h4_prior_sensitivity.tsv`、各 prior 的 RDS：信号对 H0–H4 和详细条件 SNP 后验。
- `single_causal_abf_baseline.tsv`：单因果变异对照，保留其假设标签。
- `status.txt`、`run_provenance.*`：状态、真实包版本、输入 MD5 与配置。无可信集时报告 `NO_ESTIMABLE_SIGNAL_PAIRS`，不能改写成“无生物效应”。

H0：无关联；H1/H2：仅一个性状有关联；H3：两个性状由不同变异驱动；H4：共享变异。H4 支持共享信号，不能单独证明表达介导、药物作用、受体激活方向或激素滥用。信号对结果与整区域 H4 不同。不同 prior 下结果不稳定应保留，不能挑选有利 prior。coloc 的 SNP posterior 是以 H4 为条件的后验，H4 支持不足时不用于宣称因果 SNP。

## 来源与当前验证

- [SuSiE RSS API](https://stephenslab.github.io/susieR/reference/susie_rss.html)、[coloc.susie API](https://chr1swallace.github.io/coloc/reference/coloc.susie.html)、[官方多信号教程](https://chr1swallace.github.io/coloc/articles/a06_SuSiE.html)。
- [2026 eQTL-GWAS coloc 设计评价](https://doi.org/10.1371/journal.pgen.1012141)：细胞状态、功效与共调控限制。
- [2026 multi-ancestry coloc](https://doi.org/10.1371/journal.pgen.1012221)：跨祖源需要相应 LD/分层数据，属于独立扩展路线。

2026-10-19 实际环境：susieR 0.14.2、coloc 5.2.3。合成信号用于 API 和输入守卫验证，不能作为课题分析或外部科学复制。上游开发版与本机版分别记录；未安装或升级依赖。
