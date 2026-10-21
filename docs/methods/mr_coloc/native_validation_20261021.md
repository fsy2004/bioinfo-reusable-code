# Native adapter behavior validation — 2026-10-21

最终源版本在真实本机R4.4.3中运行`Rscript tests/check_mr_coloc_native.R`，**34/34行为检查通过**，进程退出0。精确代码MD5、逐检查状态及逐方法返回状态见同名JSON。本记录使用可再生合成摘要与临时输出，不含GC研究数据；没有科学有效性、因果识别或方法统计性能验证。

| 原生接口 | 本次真实执行 | 结果 |
|---|---:|---|
| TwoSampleMR0.7.5 | 19个已核原生method ID；18个多SNP请求＋单独Wald | 18个可用方法返回；`mr_raps`因未装`mr.raps`明确failed，非成功 |
| MendelianRandomization0.10.0 | 11个不同原生函数，另对IVW/Egger/maxlik传signed LD | 11个函数均返回可重载S4；三个相关工具调用通过 |
| coloc5.2.3 | 官方`coloc.abf`，120变异quant＋cc、两组先验 | H0–H4全保留、和为1；RDS重载与直接原生结果一致 |
| coloc量化性状sdY | 显式允许原生估计 | 返回并保留相关警告；没有手写sdY估计器 |
| CLI及错误路径 | JSON→Rscript→缺RAPS | 退出2，持久化failed结果与provenance，未装依赖不伪装成功 |

TwoSampleMR实际返回的方法为`mr_wald_ratio`, `mr_two_sample_ml`, `mr_egger_regression`, `mr_egger_regression_bootstrap`, `mr_simple_median`, `mr_weighted_median`, `mr_penalised_weighted_median`, `mr_ivw`, `mr_ivw_radial`, `mr_ivw_mre`, `mr_ivw_fe`, `mr_simple_mode`, `mr_weighted_mode`, `mr_weighted_mode_nome`, `mr_simple_mode_nome`, `mr_sign`, `mr_uwr`, `mr_grip`。另一个已分派的`mr_raps`明确依赖失败。

MendelianRandomization实际返回的11个函数为`mr_ivw`, `mr_egger`, `mr_median`, `mr_maxlik`, `mr_mbe`, `mr_conmix`, `mr_lasso`, `mr_divw`, `mr_pivw`, `mr_hetpen`, `mr_cML`。验证用了100次median/MBE迭代、`mr_pivw(Boot.Fieller=FALSE)`、`mr_cML(DP=FALSE, MA=TRUE)`及较粗conmix/hetpen网格；没有宣称所有参数组合已验证。这些原生路径之间有方法重叠，不能将19＋11称30种独立算法。

拒绝路径实际覆盖：非正SE、重复SNP、未对齐EA/OA、零暴露beta、单SNP Egger、多SNP Wald、未知方法、未声明harmonization/独立性、build不一致、单位缺失或误填OR、模糊回文频率、非空输出目录。signed LD检查另实际拒绝r²/未声明signed矩阵、重排ID、奇异矩阵和不支持相关输入的TwoSampleMR。

coloc实际拒绝lead-only、声明为按关联筛选的表、坏SE、重复SNP、等位基因或坐标不一致、错误染色体、缺sdY/未允许估计、缺病例比例、无效先验、交集不足。不同变异集默认拒绝；显式允许后将119/120交集与被排除SNP记录在provenance。

首轮34检查扩展前曾有一项错误消息断言失败，由此定位R的`$`部分匹配会将缺失`sdY`误匹配`sdY_source`；已改为精确`[["sdY"]]`访问，并在最终源版本重跑全部检查通过。R启动有本机`LC_*`地区设置警告；实际包加载、调用和结果验收均完成，未为此安装或修改系统。

输入完整区域/独立工具/祖源/构建/重叠等声明仍要由真实输入审计支持。至少100变异（显式下限50）只拦截明显lead-only，不证明区域完整；F、非显著Egger或高H4也不能证明药效或介导。相关SuSiE共定位仍使用既有609，不重复实现。未实现MRPRESSO、未下载/执行第三方源码、未安装包、未修改研究仓库数据或运行GC科学分析。
