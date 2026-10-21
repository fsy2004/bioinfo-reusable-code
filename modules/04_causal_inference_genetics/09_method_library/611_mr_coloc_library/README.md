# 611 · MR与GWAS共定位方法库

入口是 [methods.json](methods.json)，按稳定方法ID检索算法、参数变体、设计、诊断及来源。源文件、假设、输入、许可、实际运行状态分别登记。覆盖目标是持续补全；当前条目不构成“所有MR算法均已实现”的声明。

```powershell
python scripts/mr_coloc.py query pleiotropy --category mr_estimator
python scripts/mr_coloc.py query SuSiE --category colocalization
python scripts/mr_coloc.py select modules/04_causal_inference_genetics/09_method_library/611_mr_coloc_library/examples/single_iv_context.json
python scripts/mr_coloc.py fetch-code --source SRC-mrcieu-twosamplemr
python scripts/mr_coloc.py verify-sources
python scripts/mr_coloc.py build --check
```

从仓库根目录运行。`query`不需要R；`select`只检查声明的操作性输入条件，输出候选和逐项阻断原因，不能证明IV有效或自动选“最佳结果”。未知完整条件的高级方法保持参考状态，需读其原文与官方输入契约。

原生运行入口见 [native_adapters.md](native_adapters.md)。支持TwoSampleMR、MendelianRandomization及coloc.abf，统一配置入口为 `python scripts/mr_coloc.py run --config <explicit.json>`；用此Python入口运行ABF时配置另加 `"engine":"coloc.abf"`。没有输入或缺包不会生成替代分析。相关IV只在支持的原生方法中使用signed r；多信号SuSiE/coloc沿用模块609。

代码保存在忽略的 `modules/_upstream_sources/mr_coloc/<source_id>/<commit>` 或 `<commit>-code`。`fetch`获取有大小上限的完整仓库常规文件快照，跳过软硬链接；`fetch-code`获取固定Git树中的代码/许可/README子集，排除数据和过大文件，明确记录遗漏。二者均不安装、不执行源码。每个文件有SHA256，子集另核Git blob SHA；`verify-sources`检查文件改动、缺失及额外文件。源码缓存不进入Git，公开库保存可重建清单、来源和自主适配器。

`record_kind=family_overview`仅作旧43方法族的导航；同一算法的多个包实现合并，参数变体单独标记。CLPP是eCAVIAR统计量别名；纯SuSiE/FINEMAP/CAVIAR等归fine_mapping；LDSC/LCV/SMR解释边界单列。获取某仓库代码不等于逐算法实现或验证完成。源码HEAD和本机已安装包版本分别记录，禁止混用两者API证据。

维护、覆盖计数与选择建议见 [方法库总览](../../../../docs/methods/mr_coloc/README.md)；机器入口为 `catalog/mr_coloc_methods.jsonl` 和 `mr_coloc_sources.jsonl`。
