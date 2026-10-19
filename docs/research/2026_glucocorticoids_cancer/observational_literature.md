# 观察性数据库论文与方法准备记录

检索日：2026-10-19（Asia/Shanghai）；时间窗：2026-01-01 至 2026-10-19。当前只做方法准备，不下载个体数据、不完成数据库选型。此清单是有目的的方法采样，不是系统综述，也不声称覆盖今年全部论文。

## 论文叫什么

总称为“基于公共数据库的观察性流行病学研究”或“二次数据分析”。NHANES 常见是复杂抽样横断面研究；链接死亡文件后可开展队列生存分析。CHARLS 多波次是纵向队列研究，单波次是横断面分析。研究药物暴露属于药物流行病学。与 GWAS/MR 组合可表述为观察性与遗传流行病学证据三角互证，需要分别说明暴露、结局和估计目标。

## 最需要保留的边界

- 糖皮质激素“使用”“长期使用”“高剂量”“不适当使用/滥用”是不同暴露，不能互相替代。
- 2024 年已经有包含 glucocorticoids 的 19 类药物 + NHANES + MR 直接先例，笼统技术组合本身难以构成创新；应先精读该先例及其补充材料。
- 2026 CKM 论文同篇结合 UKBB、NHANES、CHARLS、MR、组学和中介，但它是致编辑短文；CHARLS 仅 133 例癌症事件，MR 没有一致显著因果证据。可以学习多来源证据组织，不能把方法堆叠当成证明。
- 全部论文仅作方法来源，所列局限用于决定能否复用；本轮没有用任何论文的显著性结果推断糖皮质激素的肿瘤效应。

## 2026 年论文清单

| ID | 日期 | 论文与一手来源 | 设计/人群 | 本轮读取 |
|---|---|---|---|---|
| OBS2026_01 | 2026-08-24 | [Cardiovascular-kidney-metabolic syndrome and cancer risk: Evidence from multi-omics analysis results](https://pmc.ncbi.nlm.nih.gov/articles/PMC13508829/)；DOI 10.1515/jtim-2026-0072；PMID 42657337 | 多队列观察性分析 + MR + 组学 + 中介；致编辑短文；UKBB n=182140；NHANES n=18079；CHARLS n=7944 | full_text_main_short_letter; supplement_not_read |
| OBS2026_02 | 2026-04-24 | [Hypocholesterolemia is a consequence, not a cause, of kidney cancer: A bidirectional Mendelian randomization study with NHANES 2001–2018](https://www.sciencedirect.com/science/article/pii/S1877782126000998)；DOI 10.1016/j.canep.2026.103084；PMID 42033882 | NHANES横断面 + 双向两样本MR的证据三角互证；NHANES 102例自报肾癌、45366无癌对照；ccRCC GWAS总样本752817 | publisher_abstract_and_section_snippets; pubmed_metadata |
| OBS2026_03 | 2026-01-30 | [Association of neutrophil count with lung cancer risk and prognosis: A Mendelian randomization and NHANES analysis](https://pmc.ncbi.nlm.nih.gov/articles/PMC12863823/)；DOI 10.1097/MD.0000000000047510；PMID 41630295 | 两样本MR评估肺癌风险 + NHANES肺癌患者生存分析；MR暴露GWAS n=170702；肺癌GWAS n=27209；NHANES肺癌生存分析正文称104人 | full_text_main_methods_and_results; supplement_not_read |
| OBS2026_04 | 2026-02-27 | [Sleep Duration and Depression as Predictors of Cancer Incidence in Middle-aged and Older Adults: A Longitudinal Analysis from the China Health and Retirement Longitudinal Study (CHARLS)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12949910/)；DOI 10.21873/invivo.14258；PMID 41760316 | CHARLS纵向队列二次分析；CHARLS 2011–2020，9400人，271例新发自报癌症，随访中位9年 | PMC_search_rendered_full_text_methods_results_and_metadata; supplement_not_read |
| OBS2026_05 | 2026-06-16 | [Interacting and joint effects of frailty and depressive symptoms in relation to risk of cardio-oncology comorbidity in older adults from three prospective cohorts](https://link.springer.com/article/10.1186/s12877-026-07834-8)；DOI 10.1186/s12877-026-07834-8；PMID 42304293 | CHARLS+ELSA+HRS多队列纵向研究；三队列合计10937名≥65岁者；CHARLS 2993人 | full_text_main_methods_and_limitations; supplement_not_read |
| OBS2026_06 | 2026-05-01 | [Association between self-purchased medication use and gastric disease among Chinese adults: a cross-sectional analysis based on CHARLS 2018 data](https://link.springer.com/article/10.1186/s12889-026-27616-6)；DOI 10.1186/s12889-026-27616-6；PMID 42067841 | CHARLS 2018横断面二次分析；药物流行病学邻近方法案例；合并后19752人，最终6087人，659例胃病 | full_text_main_methods_and_limitations; supplement_not_read |
| OBS2026_07 | 2026-03-19 | [The Relationship between Sarcopenia and All-Cause and Cardiovascular Mortality Risk among Middle-Aged and Older Adults across Stages 0–3 of Cardiovascular-Kidney-Metabolic Syndrome: Evidence from NHANES and CHARLS](https://karger.com/crm/article-pdf/16/1/138/4519649/000550891.pdf)；DOI 10.1159/000550891 | NHANES死亡随访 + CHARLS外部队列复核；NHANES最终3724人，85例全因死亡/12例CVD死亡；CHARLS 2011–2020复核 | publisher_pdf_full_main_methods_and_limitations_read; supplement_not_read |
| OBS2026_08 | 2026-03-18 | [Machine learning-based association analysis of triglyceride-glucose index with melanoma prevalence and all-cause mortality: insights from cross-sectional NHANES 1999–2018 data and an external hospital-based dataset](https://www.frontiersin.org/journals/nutrition/articles/10.3389/fnut.2026.1726865/full)；DOI 10.3389/fnut.2026.1726865 | NHANES横断面患病率 + 链接死亡随访 + 内部ML预测；医院队列描述；NHANES 21360人、127例黑色素瘤；医院475例肢端黑色素瘤用于描述 | full_text_main_methods_and_abstract; supplement_not_read |

每篇逐项 methods、详细边界、DOI/PMID 和原始 URL 均见同目录 observational_papers.json 的 papers 数组。日期优先采用首个电子发表日；Sleep 论文 PMC 显示 2026-02-27，PubMed 期号为 2026 Mar-Apr，已区分。未核实全文参数的条目明确标记，不将摘要读取写成全文精读。

## 逐篇方法与复用限制

### OBS2026_01

方法：UKBB/CHARLS HR、NHANES OR；分层；蛋白与代谢风险分数；通路富集；炎性指标与MetRS中介。短文未披露完整MR估计器与所有建模参数，需补充材料。

限制：CHARLS 133例癌症事件，检验效能有限；NHANES/CHARLS结局为自报，UKBB为登记；MR未发现显著一致的遗传因果；本文属于Letter，不宜按完整方法学文章直接复刻。

### OBS2026_02

方法：按权重、分层、PSU建模；加权线性回归调整年龄、性别、BMI、他汀、吸烟、高血压和糖尿病；双向IVW及稳健MR方法。

限制：NHANES肾癌未分组织学，MR结局是ccRCC，估计对象并不完全一致；病例仅102例；摘要/片段可读，完整参数与补充材料未核读。

### OBS2026_03

方法：P<5e-8，10Mb窗口、LD阈值按文献原文需校核方向，排除中等频率回文SNP，F>100；KM、Cox、年龄分层；TwoSampleMR。

限制：肺癌发生MR与确诊者生存的估计目标不同；原文134例减29例却称104例，存在算术不一致，需核对原始流程；文中将无显著多效性当作假设被证明，不能沿用该解释；未在方法段明确完整复杂抽样处理。

### OBS2026_04

方法：Shapiro-Wilk、Mann-Whitney U/χ²；三级Cox；KM/log-rank；RCS；两段线性与似然比比较；VIF；抑郁及其他协变量分层。R 4.5.1与Zstats。

限制：PMC检索渲染可读方法与结果；最初直接打开遭验证码，后已取得正文方法；原文KM累计发生率叙述与保护方向存在文字不一致，复刻前需核图表；自报癌症/睡眠；大量完全病例排除；数据驱动阈值应探索性解释。

### OBS2026_05

方法：MICE；Cox+Schoenfeld；乘积项；RERI/AP/SI；四组联合暴露；排除首轮事件、完全病例及追加药物调整。

限制：结局是CVD与癌症共病，不等于单独肿瘤发生；队列间CES-D版本/阈值与测量存在异质性；自报疾病、基线单次暴露、残余混杂。

### OBS2026_06

方法：EF001_W4自购药二值；DA007_10自报胃病；tableone、survey、pROC、nricens、PredictABEL；三级调整模型。

限制：自购药未区分药类、剂量、疗程，不能表示糖皮质激素滥用；结局为胃病，不能当胃癌研究；大量完全病例排除；原文二值暴露却出现first quartile用语，存在描述不一致；ROC/NRI/IDI为样本内区分度，不是外部预测验证。

### OBS2026_07

方法：70/30训练测试；ABESS与Boruta选变量、VIF；10种ML、10折CV+网格调参；ROC/AUC；三级Cox、RCS、KM/log-rank，CHARLS复核全因死亡。

限制：DXA与体测方程跨队列肌肉量测量异质；CVD仅12事件，复杂模型稳定性需审查；基线单次暴露、活动/心肺适能残余混杂；结局为死亡而非癌症发生。

### OBS2026_08

方法：MICE m=5并m=20敏感性；70/30划分、10折CV；7种算法；SHAP；AUROC/PR-AUC/Brier/校准/DCA；ML只纳入样本权重，未直接纳入分层PSU。

限制：医院队列仅描述，不能称模型外部验证；横断面癌症史为患病率；空腹TyG分析的权重选择需依NHANES官方原则重新审查，不能照抄WTMEC2YR；罕见结局和缺失处理需防信息泄漏。

## 历史直接先例（不计入 2026 年数量）

[Association between 19 medication use and risk of common cancers: A cross-sectional and Mendelian randomisation study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10940964/)，2024-03-15，DOI 10.7189/jogh.14.04057，PMID 38487860。

方法：ATC分类、PSM、复杂抽样Logistic、用药时长RCS；用药开始≥癌症前6月且调查前月仍用；FinnGen复核/meta、遗传相关及反向MR。

本题边界：早已有糖皮质激素这一药类的直接先例，不能将笼统NHANES+MR当创新；用药倾向GWAS可能反映适应证/健康行为，不能等同药理效应或滥用；现用药和癌症回忆的时间重构仍受选择与回忆偏倚；公开RX无每日剂量，不能解释剂量滥用。

## 方法 → 论文 → GitHub 线索

GitHub 线索是经现场读取确认的包/框架入口；除明确说“论文使用”的包以外，只表示可选实现。没有安装或执行。具体 commit、许可证和适配检查由主任务统一收口。

| 方法 | 论文 ID | GitHub | 状态/边界 |
|---|---|---|---|
| complex_survey_design / weighted descriptive / logistic / Cox | OBS2026_02, OBS2026_06, OBS2026_08, DIRECT2024_01 | [cran/survey](https://github.com/cran/survey) | verified_GitHub_CRAN_read_only_mirror；官方主页在R-Forge；GitHub镜像不是开发上游。 |
| multiple_imputation_mice / Rubin pooling | OBS2026_05, OBS2026_08 | [amices/mice](https://github.com/amices/mice) | verified_upstream |
| Cox / KM / Schoenfeld / survival workflows | OBS2026_03, OBS2026_04, OBS2026_05, OBS2026_07, OBS2026_08 | [therneau/survival](https://github.com/therneau/survival) | verified_upstream |
| restricted cubic spline / regression modeling | OBS2026_04, OBS2026_07, OBS2026_08, DIRECT2024_01 | [harrelfe/rms](https://github.com/harrelfe/rms) | verified_upstream；复杂抽样下使用样条基函数与survey模型配合，不能直接丢弃抽样设计。 |
| propensity_score_matching | DIRECT2024_01 | [kosukeimai/MatchIt](https://github.com/kosukeimai/MatchIt) | verified_repository_redirect_or_current_maintainer_to_check；README安装开发版指向ngreifer/MatchIt；代码代理应以当前Git/API确认canonical。 |
| additive interaction RERI AP SI | OBS2026_05 | [tunsmart/interactionR](https://github.com/tunsmart/interactionR) | verified_upstream；适配coxph/glm；复杂抽样对象适配须另核。 |
| ROC/AUC | OBS2026_06, OBS2026_08 | [xrobin/pROC](https://github.com/xrobin/pROC) | verified_upstream |
| ML resampling framework | OBS2026_08 | [tidymodels/tidymodels](https://github.com/tidymodels/tidymodels) | verified_framework_option_not_claimed_paper_package |
| SHAP/model explanation | OBS2026_08 | [ModelOriented/DALEX](https://github.com/ModelOriented/DALEX) | verified_framework_option_not_claimed_paper_package |
| Firth rare-event logistic | 邻近扩展；未在本8篇中确认 | [cran/logistf](https://github.com/cran/logistf) | verified_GitHub_CRAN_mirror_adjacent_extension；上游georgheinze/logistf由镜像声明；本8篇未确认使用，不能标记为已提取方法。 |

ABESS、Boruta、SVM/神经网络/QDA/LDA/Naive Bayes/KNN、NRI/IDI、DCA、蛋白/代谢评分、中介、富集和 MR 稳健估计器的代码映射需合入主任务的 GitHub 方法目录。本文件不猜测其版本或假称已纳入代码库。分层、交互、阈值、早期事件排除和完全病例敏感性属于分析设计，通常应由本库包装脚本实现，而非为每个名词复制一个上游仓库。

## 仅限方法边界的官方变量核查

| 对象 | 官方文件 | 已核实与影响 |
|---|---|---|
| NHANES | [官方问卷/codebook](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/RXQ_RX_J.htm) | RXDUSE过去30天处方药；RXDDRUG/RXDDRGID药名/码；RXDDAYS自述使用时长；公开字段无每日剂量、频次、累计泼尼松当量；药物层需先汇总到人再按SEQN合并。 因此现有字段可研究使用/时长关联，不能直接构造临床意义的滥用。 仅以2017–2018官方codebook作方法边界示例，未完成跨周期字段审计。 |
| NHANES | [官方问卷/codebook](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/MCQ_J.htm) | MCQ220曾被医生/卫生专业人员告知癌症；MCQ230a-c癌种；MCD240a-c诊断年龄。 主要支持癌症史/患病分析；常规死亡链接不是癌症新发登记。 |
| CHARLS | [官方问卷/codebook](https://charls.pku.edu.cn/wenjian/jiahuwenjuan.pdf) | 官方问卷可读癌症部位、治疗大类；没有在本轮确证糖皮质激素具体药名/剂量字段。 不将自购药或西药治疗大类等同糖皮质激素；完整选型暂缓。 该官方PDF历史波次；不能据未检出断言所有波次均无具体用药。 |

NHANES 旧式 /2017-2018/RXQ_RX_J.htm 链接本轮返回404，已改用官方 /Nchs/Data/Nhanes/Public/2017/DataFiles/ 路径。未因未检出就断言 CHARLS 所有波次没有某字段。后续真正选库时仍须逐波次正式变量字典审计。

## 检索记录与排除

- 2026-10-19：`site:pmc.ncbi.nlm.nih.gov 2026 NHANES cancer risk`
- 2026-10-19：`site:pmc.ncbi.nlm.nih.gov 2026 CHARLS cancer risk`
- 2026-10-19：`2026 NHANES cancer "Mendelian" study`
- 2026-10-19：`2026 CHARLS cancer "cohort" study`
- 2026-10-19：`"Cardiovascular-kidney-metabolic syndrome and cancer risk"`
- 2026-10-19：`"Sleep Duration and Depression as Predictors of Cancer Incidence" 2026 DOI`
- 2026-10-19：`"Association between self-purchased medication use and gastric disease" DOI 2026`
- 2026-10-19：`site:wwwn.cdc.gov "RXQ_RX_J" "RXDDAYS"`
- 2026-10-19：`site:charls.pku.edu.cn questionnaire cancer medication`
- 2026-10-19：`GitHub primary package repositories: survey, survival, mice, rms, MatchIt, interactionR, pROC, tidymodels, DALEX`
- 2026-10-19：`"Sleep Duration and Depression" "14258" "Statistical" → PMC full-text search-rendered`

- 搜索引擎的Published/Crawled相对时间不用于正式日期；以正文/期刊/PubMed日期为准。
- 2025 cholecystectomy+CRC+NHANES+MR (10.21037/tcr-2025-1493) 排除2026清单。
- 2025 metformin+pan-cancer+NHANES+MR (10.1007/s12672-025-02021-4) 排除2026清单。
- Pre-diagnostic low muscle mass...10.1038/s41598-026-74449-9 搜到候选但本轮主站读取失败、确切日期未核，不计正式8篇。
- 未对论文做撤稿/勘误数据库专项核验；方法摘录不代表质量背书。
- GitHub仓库可读取不等于已安装、已执行或适配本题；许可、版本、计算验证由主任务统一处理。

## 下一次方法完善的优先顺序

1. 精读 2024 药物先例的糖皮质激素具体结果、工具变量 GWAS 的表型来源及补充材料，确认拟研究问题与既有工作差别。
2. 在独立 protocol 中定义“使用/长期/高剂量/滥用”、潜伏期、混杂结构和适应证；此步骤不要求现在收集数据。
3. 以复杂抽样、缺失数据、基本关联、非线性、稳健性、MR、定位/共定位建立可组合模块。预测 ML、中介及组学仅在回答额外问题时使用。
4. 不把基于用药倾向 GWAS 的 MR 解读为药物剂量或滥用效应；靶点 cis 工具、药物使用工具和内源皮质醇工具是不同研究问题。
