# =============================================================================
# run_pipeline.R · 标准分析主流程骨架 (以单细胞为例;其它分析改 Step 内容即可)
# -----------------------------------------------------------------------------
# 设计:分步 + 断点续跑(cache_step) + 每步内嵌【铁律自查点】+ 矢量出图 + 落盘统计。
# 运行:Rscript run_pipeline.R   (或在 VSCode 工作区里逐块运行)
# 这是骨架:标 TODO 处填你的真实分析;但"自查点"不要删,逐条做到再往下。
# =============================================================================
source("00_setup.R")                      # 种子/路径/框架/cache_step/log_stat

## ── Step 1. 读入数据 ────────────────────────────────────────────────────────
cat("\n== Step 1 读入数据 ==\n")
obj <- cache_step("01_raw", {
  # TODO: Seurat::Read10X(...) / readRDS(...) / read_table_smart(file.path(PROJ_ROOT,DIR_DATA,"x.csv"))
  stop("TODO: 在此读入你的数据")
})

## ── Step 2. 质控 QC ─────────────────────────────────────────────────────────
# 【自查·铁律2】先画 nFeature / nCount / mt% 分布,按【本数据】分位数/拐点定阈值;
#              勿套通用 200/20%。把最终阈值写回 config.R 的 PARAMS 并注明依据。
cat("\n== Step 2 QC ==\n")
obj <- cache_step("02_qc", {
  # TODO: 画分布图 -> save_fig(p_qc, file.path(PROJ_ROOT,DIR_FIGURES,"qc_dist"))
  stopifnot(!is.na(PARAMS$min_genes), !is.na(PARAMS$max_mito_pct))  # 逼自己先定阈值
  # TODO: subset(obj, nFeature_RNA > PARAMS$min_genes & percent.mt < PARAMS$max_mito_pct)
  obj
})
log_stat("cells_after_QC", NA)            # TODO: ncol(obj) —— 铁律6:数字由代码记录

## ── Step 3. 标准化 + 降维 ───────────────────────────────────────────────────
# 【自查·铁律1】所有随机步骤传 seed.use=SEED;n_pcs 用 ElbowPlot 核实后改 config。
# 【自查·铁律5】大矩阵:用稀疏矩阵,及时 rm()+gc(),勿 as.matrix() 全量展开。
cat("\n== Step 3 标准化/降维 ==\n")
obj <- cache_step("03_dimred", {
  # TODO: SCTransform/NormalizeData -> RunPCA -> RunUMAP(seed.use = SEED)
  obj
})

## ── Step 4. 聚类 ────────────────────────────────────────────────────────────
# 【自查·铁律3】分辨率用 clustree 看稳定性后定;勿默认 0.8 一把过。
cat("\n== Step 4 聚类 ==\n")
obj <- cache_step("04_cluster", {
  # TODO: FindNeighbors -> FindClusters(resolution = PARAMS$cluster_res, random.seed = SEED)
  obj
})

## ── Step 5. 细胞注释 ────────────────────────────────────────────────────────
# 【自查·铁律2】注释必须 marker 验证:每个簇查 top marker + 已知 canonical marker 比对;
#              自动注释(SingleR 等)只作参考,最终人工+marker 双确认。可疑簇标 "unknown"。
cat("\n== Step 5 注释 ==\n")
obj <- cache_step("05_annot", {
  # TODO: FindAllMarkers(only.pos=TRUE) -> 比对 canonical marker -> 命名
  obj
})

## ── Step 6. 下游分析 ────────────────────────────────────────────────────────
# 【自查·铁律2/3】结果要有生物学支撑;过滤假阳性(校正 p、effect size、最小细胞数);
#   · 轨迹:起点/细胞选型要有依据;  · 虚拟扰动:必须验证调控网络,不只看表型;
#   · SHAP/重要性:按亚群【分层】算,勿全样本混算。
cat("\n== Step 6 下游分析 ==\n")
res6 <- cache_step("06_downstream", {
  # TODO: 你的核心分析(DEG/通讯/轨迹/虚拟扰动/ML...)
  list()
})

## ── Step 7. 出图 (铁律4:矢量 + 期刊配色 + 要素齐全) ─────────────────────────
# 用框架:theme_pub() + pal_pub("npg"/...) + save_fig()(一次出 PDF+PNG);
# 同类细胞配色跨图统一;比例尺/统计标注/图例/分组齐全;图中文字英文。
cat("\n== Step 7 出图 ==\n")
# TODO: p <- ...+ theme_pub(); save_fig(p, file.path(PROJ_ROOT, DIR_FIGURES, "umap_annotated"))

## ── Step 8. 可复现快照 (铁律6) ──────────────────────────────────────────────
save_session()                            # 写 logs/sessionInfo.txt
cat("\n[pipeline] 完成。关键统计见 logs/key_stats.tsv,图见 figures/。\n")
cat("[pipeline] 收尾自查:对照 _framework/QUALITY_CHECKLIST.md 逐条核对后再写文稿。\n")
