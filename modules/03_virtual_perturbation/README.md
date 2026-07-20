# 03 · 虚拟扰动技术 — Virtual perturbation

本域共 26 个条目。完整字段见 [`../CATALOG.md`](../CATALOG.md)。


## 虚拟敲除与扰动模拟 — In-silico knockout

- [026_scrna_seurat_sctenifoldknk_ko.R](01_insilico_knockout/026_scrna_seurat_sctenifoldknk_ko.R) — Full QC/cluster/annotate + virtual KO pipeline
- [067_scperturb_etest.R](01_insilico_knockout/067_scperturb_etest.R) — scPerturb perturbation distance + E-test (light)
- [068_gears_combo_perturbation.py](01_insilico_knockout/068_gears_combo_perturbation.py) — GEARS single/combo perturbation prediction (GPU)
- [069_celloracle_grn_perturbation.py](01_insilico_knockout/069_celloracle_grn_perturbation.py) — CellOracle GRN virtual knockout (heavy)
- [085_squidiff_diffusion_perturbation.py](01_insilico_knockout/085_squidiff_diffusion_perturbation.py) — Squidiff/PerturbDiff diffusion perturbation (GPU)
- [494_genki_vgae_ko.py](01_insilico_knockout/494_genki_vgae_ko.py) — GenKI graph-VGAE virtual KO (KL ranking)
- [495_bulkvgk_sctenifoldknk_di.R](01_insilico_knockout/495_bulkvgk_sctenifoldknk_di.R) — Bulk co-expression virtual KO + differential influence
- [507_geneformer_insilico](01_insilico_knockout/507_geneformer_insilico) — Geneformer zero-shot embedding + in-silico deletion (baseline local)
- [561_regvelo_grn_velocity](01_insilico_knockout/561_regvelo_grn_velocity) — GRN 约束的 RNA 速率 + 逐 TF in-silico 调控子敲除,经 CellRank 命运概率重分配打分筛查转录因子
- [581_veloagent_velocity](01_insilico_knockout/581_veloagent_velocity) — 空间信息驱动的 RNA velocity 与 in-silico 敲除:scVelo + 空间 kNN 平滑基线 + veloAgent 守卫式封装,出速度场 quiver / raincloud / lollipop

## 基因调控网络推断 — GRN inference

- [047_rcistarget_tf_motif_network.R](02_grn_inference/047_rcistarget_tf_motif_network.R) — RcisTarget motif/TF enrichment + regulatory network
- [081_pyscenic_regulon_tf_activity.py](02_grn_inference/081_pyscenic_regulon_tf_activity.py) — pySCENIC GRN + ctx + AUCell wrapper
- [511_tf_convergence_depmap_jaspar](02_grn_inference/511_tf_convergence_depmap_jaspar) — Three-evidence convergence to core TFs
- [582_dspin_regulatory_network](02_grn_inference/582_dspin_regulatory_network) — 从多重扰动 scRNA-seq 反推程序级 Ising 自旋网络(共享耦合矩阵 J + 每扰动场向量 h),自带混池相关 / 朴素平均场两条带真值评分的基线,D-SPIN 正式伪似然求解器走守卫式调用。
- [583_kegni_knowledge_grn](02_grn_inference/583_kegni_knowledge_grn) — 知识图增强 GRN 推断:5 种本机可跑基线(Pearson/Spearman/PCA 嵌入点积/纯知识先验/知识-表达秩融合)按 BEELINE 口径出 EPR/AUPRC/AUROC 榜单,并对上游 KEGNI 深度模型做守卫式 CLI 封装(不臆造 Python API)。
- [584_cellpolaris_grn_transfer](02_grn_inference/584_cellpolaris_grn_transfer) — 在已有 GRN 上建高斯概率图模型做 TF 虚拟敲除(ΔX),并用 ΔX 与真实相邻状态表达差的余弦相似度沿分化轨迹排主控 TF;上游迁移学习生成 GRN 段为守卫式封装。
- [585_ignite_grn_inference](02_grn_inference/585_ignite_grn_inference) — 用非对称动力学 Ising 模型的反问题（IGNITE，PLoS Comput Biol 2026）从未扰动的伪时序单细胞数据反推有向有符号 GRN 并模拟基因敲除，自带 3 个本机可跑的朴素 GRN 基线（Pearson / GraphicalLassoCV 偏相关 / 滞后岭回归）做 AUROC-AUPRC 边恢复对照。
- [586_psgrn_grn_inference](02_grn_inference/586_psgrn_grn_inference) — 586 · PSGRN — 从带干预标签的单细胞扰动矩阵(CRISPRi/Perturb-seq 风格,含 non-targeting 对照)推断有向基因调控网络。忠实复现上游算法:相关性造合成金标准 → 4 个扰动特征 → LightGBM 自训练重排全部有序基因对;内置两条朴素基线(共表达 |Pearson|、单变量干预效应)做强制对照,出 PR 曲线 / precision@K dumbbell / 打分热图。上游官方 CausalBench 评测入口做守卫式封装,缺包时打印真实命令、不伪造返回值。
- [587_regformer_grn_mamba_fm](02_grn_inference/587_regformer_grn_mamba_fm) — RegFormer GRN 重建评测台：把上游「基因嵌入→余弦相似度 TF 有向图→谱聚类模块」下游链路本地复刻，用共表达/PCA 朴素嵌入作必跑基线，可插入任意外部 gene_embedding.npy 在同一口径下对比。

## 因果表示与反事实 — Causal & counterfactual

- [588_sccausalvi_causal_perturbation](03_causal_perturbation/588_sccausalvi_causal_perturbation) — 案例-对照单细胞扰动响应的因果解耦模块：默认跑可复现线性基线（条件中心化 PCA 背景表示 + 全局/细胞类型特异 Δ 反事实 + kNN 响应细胞打分），scCausalVI 深度模型为守卫式可选路径。

## 药物扰动与响应 — Drug perturbation

- [070_chemcpa_drug_perturbation.py](04_drug_perturbation/070_chemcpa_drug_perturbation.py) — chemCPA drug-perturbation expression prediction (GPU)
- [071_scdrug_response_prediction.py](04_drug_perturbation/071_scdrug_response_prediction.py) — scDrug single-cell drug response (heavy)
- [518_beyondcell_drug_response](04_drug_perturbation/518_beyondcell_drug_response) — beyondcell core re-impl: BCS + therapeutic clusters
- [589_scdruglink_drug_response](04_drug_perturbation/589_scdruglink_drug_response) — 按 scDrugLink 上游源码复刻的单细胞药物重定位打分:Drug2Cell 靶点臂(促进/抑制)× 扰动签名臂(敏感/耐药)在细胞类型层面 exp(weight) 相乘串联,输出全图谱与细胞类型两级治疗评分排序 + AUROC/AUPR 对照。

## 扰动预测基准 — Perturbation benchmarks

- [590_scperturbench_generalization](05_benchmark/590_scperturbench_generalization) — 把自己的单细胞扰动预测按 scPerturBench (Nat Methods 2026) 的同一套指标打分，并强制与上游口径的朴素基线 (controlBaseline / trainMean) 对照，判断深度模型是否真的赢过"什么都不预测"。
- [591_scarchon_perturbation_benchmark](05_benchmark/591_scarchon_perturbation_benchmark) — 留一批次(leave-one-batch-out)的单细胞扰动响应预测基准骨架:按 scArchon 口径评估预测的扰动后表达,强制与 control/mean 朴素基线同台对比
