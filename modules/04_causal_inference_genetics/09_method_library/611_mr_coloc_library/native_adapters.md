# Native MR and coloc.abf adapters

显式配置调用已经安装的官方R包。不自动安装、下载、执行缓存源码、harmonize、挑选工具、生成研究数据或绘图。输入检查和原生程序返回不证明科学有效性；用药liability、内源激素、受体表达与药物剂量分别解释。

```powershell
Rscript run_native_mr.R --config C:/analysis/mr_config.json
Rscript run_coloc_abf.R --config C:/analysis/coloc_config.json
```

所有相对路径相对于JSON所在目录；输出目录须为空，既有结果不会覆盖。成功退出0，配置拒绝、原生失败或有方法失败退出2。包缺失时报告错误，不安装。可从R中`source()`后调用`run_native_mr(config_path)`或`run_coloc_abf(config_path)`取得同一结果对象。

## MR schema and configuration

CSV每行一个已harmonized、预先选择的双样本工具：

| 列 | 约定 |
|---|---|
| `snp` | 非空唯一稳定ID；不自动去重 |
| `beta_exposure`, `se_exposure` | 有限beta、严格正SE；beta不能为0 |
| `beta_outcome`, `se_outcome` | 有限beta、严格正SE；二分类结局须logOR beta，不能填OR或把OR的SE直接冒充logOR SE |
| `effect_allele_exposure`, `other_allele_exposure` | 大写单碱基A/C/G/T；两列不同 |
| `effect_allele_outcome`, `other_allele_outcome` | 必须与暴露的EA/OA完全一致；不自动翻转或换链 |
| `eaf_exposure`, `eaf_outcome` | 非回文时可省；若提供须有限且在(0,1)。回文SNP两列必需，且需`palindromes_resolved=true`；0.42–0.58附近频率或两者差>0.1拒绝 |

此最小适配器不收indel、多等位形式、未liftOver或未对齐汇总；应在合适的准备模块中明确处理后再输入。列名和SE/variance语义由用户核实；原始来源、构念、单位、祖源、构建必须填写。

```json
{
  "engine": "TwoSampleMR",
  "input": "harmonized.csv",
  "methods": ["mr_ivw", "mr_egger_regression", "mr_weighted_median"],
  "seed": 20261021,
  "outdir": "results/mr_run01",
  "harmonized": true,
  "independent_instruments": true,
  "metadata": {
    "exposure": {
      "construct": "REPLACE_WITH_ACTUAL_EXPOSURE_CONSTRUCT",
      "unit": "REPLACE_WITH_BETA_SCALE",
      "ancestry": "EUR",
      "build": "GRCh37",
      "source": "REPLACE_WITH_RELEASE_URL_AND_INPUT_SHA256"
    },
    "outcome": {
      "construct": "REPLACE_WITH_ACTUAL_ENDPOINT",
      "unit": "log odds",
      "ancestry": "EUR",
      "build": "GRCh37",
      "source": "REPLACE_WITH_RELEASE_URL_AND_INPUT_SHA256"
    },
    "sample_overlap": "REPLACE_WITH_AUDIT_OR_EXPLICIT_UNKNOWN"
  }
}
```

构建仅接受`GRCh37`或`GRCh38`，两侧声明的build/ancestry须一致。相同祖源/build不证明样本无重叠；未知可明确写unknown，不自动当0。`independent_instruments=true`只记录已经完成独立性审计的声明，适配器不会替代LD剪枝或MHC/多效性审计。

### TwoSampleMR

本机0.7.5的19个原生method ID均可分派，调用`mr_method_list()`核实接口，逐方法直接调用官方函数，参数从`default_parameters()`取默认。

| 输入数 | 原生ID |
|---|---|
| 恰好1 | `mr_wald_ratio` |
| 至少2 | `mr_two_sample_ml`, `mr_ivw`, `mr_ivw_radial`, `mr_ivw_mre`, `mr_ivw_fe`, `mr_uwr` |
| 至少3 | `mr_egger_regression`, `mr_egger_regression_bootstrap`, `mr_simple_median`, `mr_weighted_median`, `mr_penalised_weighted_median`, `mr_simple_mode`, `mr_weighted_mode`, `mr_weighted_mode_nome`, `mr_simple_mode_nome`, `mr_raps`, `mr_grip` |
| 至少6 | `mr_sign`，输出为符号一致性评分，不是因果beta |

单SNP必须单独运行Wald，不能与Egger等一起请求；多SNP请求Wald也拒绝。数量只是软件输入下限，不代表足够功效或可靠估计。TwoSampleMR不接受本适配器的相关工具输入。

可选`parameters`覆盖官方默认中已有的键，例如`{"nboot":1000}`。正式分析须根据方法要求预设参数；测试用100次bootstrap只是行为验证，不能当科学标准。当前`mr.raps`未安装：请求`mr_raps`产生明确`failed`结果、错误信息和`partial_failure` provenance，其余合格请求继续保存；绝不把缺依赖记为成功。

### MendelianRandomization

将`engine`改为`MendelianRandomization`，支持11个真实单变量函数：

`mr_ivw`, `mr_egger`, `mr_median`, `mr_maxlik`, `mr_mbe`, `mr_conmix`, `mr_lasso`, `mr_divw`, `mr_pivw`, `mr_hetpen`, `mr_cML`。

IVW/maxlik至少2 SNP，其余至少3；无单SNP替代。`mr_hetpen`子集组合昂贵，适配器限制≤15 SNP，仅是计算保护。`mr_cML`需`metadata.exposure.sample_size`，适配器传给官方`n`。方法参数放`method_args`下对应ID，不接受未知参数或覆盖`object/correl/seed/random_seed/n`。

```json
{
  "method_args": {
    "mr_median": {"weighting": "weighted", "iterations": 10000},
    "mr_mbe": {"weighting": "weighted", "iterations": 10000},
    "mr_cML": {"DP": false, "MA": true}
  }
}
```

seed传给包含原生seed/random_seed参数的函数；其他函数用R统一种子。原生函数的其他选项按其安装版本解释，未在本测试覆盖的参数组合不称已验证。

仅`mr_ivw/mr_egger/mr_maxlik`支持显式signed LD。此时设置`independent_instruments=false`，并提供：

```json
{
  "ld": {
    "file": "signed_ld.csv",
    "signed_r": true,
    "effect_allele_aligned": true,
    "reference": "REPLACE_WITH_REFERENCE_RELEASE_AND_SOURCE",
    "ancestry": "EUR",
    "build": "GRCh37"
  }
}
```

LD CSV首列必须叫`snp`，其值及后续列名均按摘要输入SNP顺序，严格一致；不静默重排。值为按输入effect allele定向的signed correlation **r**，不能用r²；对称、单位对角、有限且正定。LD参考祖源/build须匹配；等位基因对应是用户经核验的声明，纯数值检查不能证明它。相关输入禁止robust/penalized选项，避免原生内部改选项而产生误读。

### MR outputs

- `native_results.rds`：原生`MRInput`（适用时）、每方法原生list/S4结果及provenance，保留native failure。
- `results.csv`：每个请求的returned/failed、原生estimate/SE/P/CI、输入SNP数、原生报告的nsnp、警告和错误。原生未提供CI或P则留空，不手算补齐；`mr_sign`明确为评分。
- `native_fields.csv`：原生字段的长表，完整复杂对象仍以RDS为准。
- `provenance.json`、`sessionInfo.txt`：配置、输入及适配代码MD5、实际版本、seed、LD/独立性声明、状态。

## coloc.abf schema and configuration

各CSV是同一预设完整基因组区域的未按P值筛选、未clump的摘要，列为：

`snp,chromosome,position,beta,se,effect_allele,other_allele,eaf`。

要求非空唯一SNP、正整数坐标、有限beta/正SE、EAF在(0,1)、大写单碱基且EA≠OA。两个文件共享SNP必须坐标及EA/OA一致；不自动翻转。`chromosome`用与config完全相同的字符串，例如`1`，不自动去chr前缀。只按SNP交集且记录每侧缺失，拒绝静默丢行。适配器以SE平方传官方`varbeta`，不将SE直接当variance。

```json
{
  "trait1": {
    "name": "quantitative_trait",
    "input": "trait1_full_region.csv",
    "type": "quant",
    "N": 100000,
    "sdY": 1,
    "sdY_source": "REPLACE_WITH_MEASURED_VARIANCE_OR_STANDARDIZATION_EVIDENCE",
    "full_region": true,
    "association_filtered": false,
    "unclumped": true,
    "metadata": {
      "construct": "REPLACE_WITH_TRAIT1",
      "unit": "SD",
      "ancestry": "EUR",
      "build": "GRCh37",
      "source": "REPLACE_WITH_FULL_SUMMARY_RELEASE_AND_SHA256"
    }
  },
  "trait2": {
    "name": "case_control_endpoint",
    "input": "trait2_full_region.csv",
    "type": "cc",
    "N": 120000,
    "case_fraction": 0.25,
    "n_cases": 30000,
    "full_region": true,
    "association_filtered": false,
    "unclumped": true,
    "metadata": {
      "construct": "REPLACE_WITH_PRECISE_ENDPOINT",
      "unit": "log odds",
      "ancestry": "EUR",
      "build": "GRCh37",
      "source": "REPLACE_WITH_FULL_SUMMARY_RELEASE_AND_SHA256"
    }
  },
  "region": {"chromosome": "1", "start": 1000000, "end": 2000000},
  "metadata": {"sample_overlap": "REPLACE_WITH_AUDIT_OR_EXPLICIT_UNKNOWN"},
  "alleles_aligned": true,
  "allow_partial_overlap": false,
  "min_variants": 100,
  "min_overlap_fraction": 0.8,
  "prior_grid": [
    {"p1": 0.0001, "p2": 0.0001, "p12": 0.00001},
    {"p1": 0.0001, "p2": 0.0001, "p12": 0.000001}
  ],
  "seed": 20261021,
  "outdir": "results/coloc_run01"
}
```

上面N、区域、case_fraction等都是格式示意，不是项目实际数据。quant须正sdY及其来源；若sdY确实未知，可去掉这两键并明确`estimate_sdY=true`，由官方coloc以N、MAF与varbeta估计，警告和假设保留。cc须显式病例比例，若还给n_cases则检查与N一致。群体患病率不能替代GWAS样本病例比例。

默认区域与交集至少100变异，允许显式设为≥50；这是拒绝lead-only等明显错误的输入保护，**不是区域完整性、精度或科学有效性证明**。full_region/unfiltered/unclumped及source记录需要用户依据真实全文件审计填写，不能将筛选表声明成完整GWAS。若两侧变异集不同，须主动允许partial overlap，且交集占每侧比例均达到阈值；被排除SNP完整列入provenance。不会删除坏SE、重复或不一致等位基因以“修好”输入。

`prior_grid`必须事先指定，分别运行官方`coloc.abf`，不自动挑最大H4。先验要求正、p12≤min(p1,p2)，且区域变异数×各先验<1，防止原生静默调整过大先验。ABF假定区域内每个性状至多一个因果变异；不需要LD输入，不能代表多信号SuSiE。需要signed LD与多信号时用已有609模块。H4不证明因果方向、介导或药效，两个性状LD结构/信号/功效差异仍需科学审查。

输出`native_coloc.rds`、`posterior_h0_h4.csv`、各prior的原生variant结果CSV、`provenance.json`及`sessionInfo.txt`；保留H0–H4，不只输出H4。官方函数仍可能因模型/数值条件失败，失败不生成伪科学结果。

## Actual behavior verification

```powershell
Rscript tests/check_mr_coloc_native.R
```

从仓库根目录运行；脚本使用可再生合成摘要和临时目录，调用真实R4.4.3、本机TwoSampleMR0.7.5、MendelianRandomization0.10.0、coloc5.2.3，无GC或人类数据。逐函数实际返回/失败、输出RDS重载、原生结果对应、CLI退出码、输入拒绝见`docs/methods/mr_coloc/native_validation_20261021.json/.md`。没有把exports数量当作运行验证，也没有验证所有参数组合、假设、统计性能或本课题的科学结论。

API依据：本机已安装包的`formals/getExportedValue`及官方实现；[TwoSampleMR reference](https://mrcieu.github.io/TwoSampleMR/reference/index.html)、[MendelianRandomization reference](https://amymariemason.github.io/MR/reference/index.html)、[coloc.abf reference](https://chr1swallace.github.io/coloc/reference/coloc.abf.html)。不复制估计器源码，不在适配器手写估计式。
