# 工具选型指南 · 何时用哪个（避免乱套乱用）

> 按分析目的分组;每组给候选工具的**优劣 + 选型建议 + 对应本库模块**。
> ✓ = 本库已有 turnkey 模块;⏳ = 可补(中优先,非缺失);⚠️ = 避坑。
> DL/foundation-model 的优劣判断有 2025–26 基准实证(见末尾"避坑" + 桌面 `research-frameworks`)。

## 0. 三条总原则
1. **经典稳工具做主力,DL/foundation model 当假设生成器/加分项**(基准:FM 常打不过 PCA/Harmony/线性基线)。
2. **关键结论多工具取交集**(通讯/扰动/轨迹),比单工具抗审稿。
3. **一切虚拟/in-silico 用 hypothesis-gen 措辞 + 配简单基线**(底线见记忆 no-unfounded-claims)。

---

## 1. 单细胞:整合 / 注释 / 聚类
| 任务 | 首选 | 备选 / 何时用 | 模块 |
|---|---|---|---|
| 全流程 | **Seurat** | 标准主力,生态最全 | ✓026/046/049 |
| 批次整合 | **Harmony** | 快稳默认;scVI 仅超大/复杂时;CCA/RPCA 跨条件 | ✓ |
| 去双细胞 | scDblFinder | 标配 | ✓ |
| 注释 | **marker 手工(canonical)** | 最可靠必做;SingleR 自动作辅助 | ✓046 |
| ⚠️ | — | Geneformer/scGPT zero-shot 注释/整合**常打不过 Harmony**,搁置 | research-frameworks |

## 2. 轨迹 / 拟时序 / 命运
| 任务 | 工具 | 优劣 / 何时用 | 模块 |
| 拟时序 | Monocle2/3 (DDRTree) | 经典,分支轨迹 | ✓044/082 |
| 干性/起点 | **CytoTRACE2** | 无需预设根,定起始细胞 | ✓082 |
| 分支+沿轨迹DE | Slingshot/tradeSeq | 稳 | ✓082 |
| 命运概率 | Palantir / CellRank2 | 终态概率 + 驱动基因 | ✓082/072 |
| 矢量场 | scTour(DL) / VECTOR | 方向流场 | ✓062 / ⏳VECTOR |
| ★选型 | **三算法收敛定起始细胞**(CytoTRACE+Monocle+VECTOR/scTour 取一致) — PDAC 范文打法,抗单算法质疑 | |

## 3. 细胞通讯
| 工具 | 优劣 / 何时用 | 模块 |
| **CellChat** | 主力,DB 全,可视化好(圈/弦/气泡) | ✓051 |
| CellPhoneDB | 置换检验严格;与 CellChat 双跑取交集更稳 | |
| NicheNet | 配体→下游靶基因(调控潜能),补因果 | |
| ★选型 | 常规 CellChat;抗质疑 = **CellChat+CellPhoneDB 交集**;机制闭环 + NicheNet(→UCell→GSEA→Venn 重叠) | ⏳闭环 |

## 4. 虚拟扰动（你的强项方向）
| 工具 | 优劣 / 何时用 | 模块 |
| **CellOracle** | GRN + in-silico KO,转录层主力,成熟 | ✓069 |
| GenKI / scTenifoldKnk | 共表达网络 KO | ✓494/495/026 |
| scTenifoldNet | 过表达 OE | |
| GEARS / Geneformer-scGPT in-silico | DL 扰动;⚠️新基准下**不稳超线性基线**,搁置/必配 baseline | ⏳068 / research-frameworks |
| ★选型 | 转录层 KO = CellOracle;**双引擎取共识**(CellOracle + FM)= AI 路线;任何 DL 扰动**必配简单基线** | |

## 5. 空间转录组
| 任务 | 首选 | 备选 | 模块 |
| 基础 | Seurat/Semla | 聚类/注释/可视化 | ✓027/050 |
| 反卷积 | **RCTD** | cell2location(贝叶斯)/Tangram(映射) | ✓505 |
| 空间域 | SpaGCN/GraphST(DL) | BayesSpace | |
| 多视图空间依赖 | MISTy | 多尺度关系 | ✓505 |
| 生态位 | NMF(分解 RCTD 比例) | 数据驱动 niche | ✓505 |
| 邻接共定位 | CellDegree(KNN) | 物理界面 | ✓505 |
| 空间多组学 | SpatialGlue | RNA+蛋白 | ⏳ |

## 6. MR / GWAS（因果,你的强项）
| 任务 | 工具 | 模块 |
| 主体 MR | **TwoSampleMR**(IVW/Egger/加权中位/众数) | ✓032/033 |
| 异常值 | MR-PRESSO | ✓ |
| 共定位 | coloc | ✓075 |
| 多变量/中介 | MVMR / 两步中介 MR(Sobel/Delta) | ✓079 / ⏳两步中介 |
| 细胞型特异 | csMR(单核 eQTL) | ⏳ |
| 单细胞 TWAS | sc-TWMR | (SSc 用过) |
| ★选型 | 基础因果 TwoSampleMR;中介 = 两步 MR;细胞分辨率 = csMR/sc-TWMR |

## 7. ML 特征 / 诊断 / 预后
| 任务 | 工具 | 模块 |
| 特征选择 | LASSO/SVM-RFE/RF/Boruta/12-ML | ✓012-015/034 |
| ★稳健筛选 | **三法投票**(拓扑×相关×Boruta,取≥2) | ✓502 |
| 诊断模型 | ROC/校准/DCA/nomogram | ✓016/063 |
| 多算法预后 | Mime(101 combo)/IRLS | ✓496 |
| 解释 | SHAP(按亚群分层) | ✓052 |
| ★泛化诚实 | **REML meta + LODO + meta-score** | ✓503 |

## 8. 富集 / WGCNA / TF / 多组学 / 免疫 / 代谢 / 药敏 / 对接
| 任务 | 工具 | 模块 |
| 富集 | clusterProfiler(GO/KEGG)/GSEA;通路活性 GSVA/UCell/decoupler | ✓007/076 |
| 共表达 | bulk WGCNA / 单细胞 hdWGCNA | ✓054 / ✓504 |
| TF 调控 | SCENIC/pySCENIC;★+ DepMap CRISPR + JASPAR 位点收口 | ✓047/081 / ⏳收口 |
| 多组学整合 | MOFA/DIABLO;NMF/ConsensusCluster 分型 | ✓083/084 |
| 免疫浸润 | CIBERSORT/IOBR(多方法) | ✓017-021/492 |
| 代谢 | scMetabolism / scFEA | ⏳ |
| 药敏 | beyondcell / oncoPredict | ⏳ |
| 对接 / MD | AutoDock Vina + GROMACS + MM-PBSA(PCA-FEL) | ✓022/086 |

## 9. ⚠️ 避坑清单（2025–26 基准实证）
- **DL foundation model = 假设生成器,非主力**:zero-shot 嵌入常打不过 PCA/Harmony/scVI(Genome Biol 2025);扰动预测打不过线性基线(Nat Methods 2025)。→ 用 FM **必配简单基线 + 别把结论压单一 FM**。
- **单细胞 FM 无 scaling law、早饱和** → 不必自研预训练。
- 通讯/扰动/轨迹:**多工具取交集 > 单工具**。
- 注释:**marker 验证必做**,别裸信自动/FM。
- 一切 in-silico:**hypothesis-gen 措辞**。

## 10. 还缺什么?（诚实评估,2026-06-25）
- **核心栈已完整**:单细胞 / 空间 / 轨迹 / 通讯 / 虚拟扰动 / MR / ML / 富集 / WGCNA / TF / 多组学 / 免疫 / 对接-MD,加新补的 502–505,常见纯干二区打法全覆盖。
- **⏳ 中优先(非缺失,锦上添花,需要时再建)**:两步中介 MR/csMR、细胞通讯功能复现闭环、TF DepMap/JASPAR 收口、scMetabolism/scFEA、beyondcell/BayesPrism、VECTOR、SpatialGlue。
- **AI/DL 路线**:方案已备(`research-frameworks`),SSc 投稿处理后再启动。
