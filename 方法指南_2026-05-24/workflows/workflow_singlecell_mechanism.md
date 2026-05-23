# workflow_singlecell_mechanism

用途：单细胞注释、疾病相关细胞识别、细胞通讯、轨迹/命运推断、虚拟扰动和空间验证路线。

## 适用场景

- 已有 10x、Seurat RDS 或 h5ad 单细胞数据。
- 需要解释候选基因在哪类细胞、在哪个状态、通过哪些细胞互作发挥作用。
- 想把文章从 UMAP 描述推进到通信、轨迹、调控和扰动。

## 输入规范

| 输入 | 推荐格式 | 关键字段 |
|---|---|---|
| 单细胞对象 | Seurat RDS 或 h5ad | counts、metadata、sample、group |
| marker 表 | CSV/TSV | celltype、gene |
| bulk 表型 | 表达矩阵+表型向量 | sample、phenotype |
| 空间对象 | Visium 文件或 AnnData | spatial coordinates |
| 候选基因 | TXT/CSV | gene |

## 推荐顺序

| 步骤 | 脚本 | 核心输出 | 说明 |
|---:|---|---|---|
| 1 | [023_RDS对象_结构检查.R](<C:/Users/fsy/Desktop/原始代码/按分析用途分类/08_单细胞_空间转录组_细胞轨迹/023_RDS对象_结构检查.R>) | 对象结构 | 先确认 assay、metadata、降维是否存在 |
| 2 | [024_单细胞RDS_数据整理.R](<C:/Users/fsy/Desktop/原始代码/按分析用途分类/08_单细胞_空间转录组_细胞轨迹/024_单细胞RDS_数据整理.R>) 或 [025_单细胞_数据整理.R](<C:/Users/fsy/Desktop/原始代码/按分析用途分类/08_单细胞_空间转录组_细胞轨迹/025_单细胞_数据整理.R>) | 标准 Seurat 对象 | 数据入口 |
| 3 | [049_单细胞_手工注释_CellChat_轨迹.R](<C:/Users/fsy/Desktop/原始代码/按分析用途分类/08_单细胞_空间转录组_细胞轨迹/049_单细胞_手工注释_CellChat_轨迹.R>) | 注释 UMAP、marker、轨迹 | 主模板 |
| 4 | [046_单细胞_单基因发表级图.R](<C:/Users/fsy/Desktop/原始代码/按分析用途分类/08_单细胞_空间转录组_细胞轨迹/046_单细胞_单基因发表级图.R>) | 目标基因图 | 候选基因定位 |
| 5 | [058_单细胞_Scissor疾病相关细胞.R](<C:/Users/fsy/Desktop/原始代码/按分析用途分类/08_单细胞_空间转录组_细胞轨迹/058_单细胞_Scissor疾病相关细胞.R>) | Scissor+/- 细胞 | bulk 到 scRNA 桥接 |
| 6 | [051_单细胞_CellChat细胞通讯.R](<C:/Users/fsy/Desktop/原始代码/按分析用途分类/08_单细胞_空间转录组_细胞轨迹/051_单细胞_CellChat细胞通讯.R>) | LR 网络和 bubble | 细胞通讯 |
| 7 | [077_NicheNet_配体靶基因通信推断.R](<C:/Users/fsy/Desktop/原始代码/按分析用途分类/16_空间通讯_细胞命运/077_NicheNet_配体靶基因通信推断.R>) | ligand activity、ligand-target | 通讯机制深化 |
| 8 | [062_scTour_拟时序向量场教程.py](<C:/Users/fsy/Desktop/原始代码/按分析用途分类/08_单细胞_空间转录组_细胞轨迹/062_scTour_拟时序向量场教程.py>) 或 [072_CellRank_命运概率与驱动基因.py](<C:/Users/fsy/Desktop/原始代码/按分析用途分类/16_空间通讯_细胞命运/072_CellRank_命运概率与驱动基因.py>) | pseudotime、fate、driver genes | 轨迹/命运 |
| 9 | [081_pySCENIC_Regulon_TF活性.py](<C:/Users/fsy/Desktop/原始代码/按分析用途分类/13_转录因子调控_基因组圈图/081_pySCENIC_Regulon_TF活性.py>) 或 [076_decoupler_TF通路活性评分.py](<C:/Users/fsy/Desktop/原始代码/按分析用途分类/16_空间通讯_细胞命运/076_decoupler_TF通路活性评分.py>) | TF/regulon 活性 | 调控解释 |
| 10 | [068_GEARS_单细胞组合扰动预测.py](<C:/Users/fsy/Desktop/原始代码/按分析用途分类/14_单细胞虚拟扰动_扰动数据库/068_GEARS_单细胞组合扰动预测.py>)、[069_CellOracle_GRN虚拟扰动.py](<C:/Users/fsy/Desktop/原始代码/按分析用途分类/14_单细胞虚拟扰动_扰动数据库/069_CellOracle_GRN虚拟扰动.py>) 或 [085_Squidiff_扩散模型单细胞扰动预测.py](<C:/Users/fsy/Desktop/原始代码/按分析用途分类/14_单细胞虚拟扰动_扰动数据库/085_Squidiff_扩散模型单细胞扰动预测.py>) | 扰动预测 | 高级创新 |
| 11 | [074_Tangram_单细胞到空间映射.py](<C:/Users/fsy/Desktop/原始代码/按分析用途分类/16_空间通讯_细胞命运/074_Tangram_单细胞到空间映射.py>)、[080_cell2location_Squidpy_空间生态位.py](<C:/Users/fsy/Desktop/原始代码/按分析用途分类/16_空间通讯_细胞命运/080_cell2location_Squidpy_空间生态位.py>) 和 [073_COMMOT_空间细胞通讯.py](<C:/Users/fsy/Desktop/原始代码/按分析用途分类/16_空间通讯_细胞命运/073_COMMOT_空间细胞通讯.py>) | 空间映射、邻域、空间通讯 | 空间验证 |

## 结果链条写法

细胞分群 -> marker 手工注释 -> 疾病相关细胞 -> 细胞通讯 -> 配体-受体-靶基因 -> 轨迹/命运 -> TF/regulon -> 虚拟扰动 -> 空间生态位验证。

## 热点组合

- 巨噬细胞/成纤维细胞/T 细胞互作：CellChat + NicheNet + COMMOT。
- 细胞命运和虚拟扰动：CellRank/scTour + CellOracle/GEARS/Squidiff。
- 空间生态位：Tangram/cell2location + Squidpy + COMMOT。

## 质量控制检查

- 每个样本和分组在 metadata 中保留。
- 自动注释必须用 marker 手工校验。
- 通讯分析需要足够细胞数，低丰度细胞类型应合并或谨慎解释。
- 轨迹 root/terminal state 要有生物学依据。
- 虚拟扰动结果必须用外部扰动数据、空间定位或实验设计支撑。
