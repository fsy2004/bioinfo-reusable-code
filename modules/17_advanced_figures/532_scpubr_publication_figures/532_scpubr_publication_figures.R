# =============================================================================
# 编号   : 532
# 脚本名 : SCpubr 单细胞出版级图一键出(UMAP/dotplot/violin/alluvial)+ 色盲安检
# 分类   : 17_advanced_figures
# 用途   : 用 SCpubr(v3.x)把单细胞标准图一键升级为顶刊审美、色盲友好:
#          do_DimPlot(UMAP)/ do_DotPlot(marker 点图)/ do_FeaturePlot(精修表达)/
#          do_ViolinPlot / do_AlluvialPlot(细胞流向桑基),并跑 do_ColorBlindCheck
#          对所用配色做 deutan/protan/tritan 三型色盲模拟安检。
# ★诚实基线 : 这是【出图/审美标准化】工具,不是【分析方法】。底层聚类/表达值完全不变,
#          SCpubr 只替换各包简陋默认图(Seurat::DimPlot 灰底彩点 / DotPlot 默认主题)
#          为统一刊级排版 + 色盲安全配色。本脚本同时画 Seurat 默认 UMAP 作为对照
#          (assets/00_baseline_seurat_default_umap),让"标准化前后"可并排目检 ——
#          诚实地说:它不改变结论,只改变图的可读性与投稿通过率。审美与库 theme_pub 一致。
# 高级图  : 全部非平凡条形图 —— 点图(dot)/小提琴(violin)/桑基(alluvial)/散点降维(UMAP)。
# 依赖   : SCpubr(>=3.0) · Seurat · ggplot2 ;均已安装。
# 运行   : Rscript 532_scpubr_publication_figures.R                       # 合成示例,零改动即跑
#          Rscript 532_scpubr_publication_figures.R --input my_seurat.rds --celltype celltype --condition condition
# 输入   : --input = Seurat 对象 .rds(含 UMAP 降维 + 元数据列);留空=脚本内生成合成 demo。
#          --celltype/--condition = 元数据中的细胞类型 / 分组列名。
# =============================================================================

## ---- 定位框架并载入顶刊主题 -------------------------------------------------
.find_fw <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  d <- if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
  for (i in 0:6) { cand <- file.path(d, paste(rep("..", i), collapse = "/"), "_framework", "theme_pub.R")
    if (file.exists(cand)) return(normalizePath(cand)) }
  stop("未找到 _framework/theme_pub.R")
}
source(.find_fw())
set.seed(42)
suppressWarnings(suppressMessages({ library(Seurat); library(SCpubr); library(ggplot2) }))

SCRIPT_DIR <- bio_script_dir()
DDAT   <- file.path(SCRIPT_DIR, "example_data")
ASSETS <- file.path(SCRIPT_DIR, "assets")
args <- bio_args(list(
  input     = file.path(DDAT, "synthetic_seurat.rds"),
  celltype  = "celltype",     # 元数据中的细胞类型列
  condition = "condition",    # 元数据中的分组列(用于 split / 桑基)
  outdir    = file.path(SCRIPT_DIR, "results")))
for (d in c(DDAT, ASSETS, args$outdir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

## 期刊离散配色(NPG),色盲友好 —— 同时供 SCpubr colors.use 与基线对照复用
pal_for <- function(levels) { v <- pal_pub(length(levels), "npg"); stats::setNames(v, levels) }

## ---- 0. 合成示例 Seurat 对象(synthetic, for demo only)----------------------
# 4 类细胞(Tcell/Bcell/Mono/NK),各有 5 个特异 marker 基因;2 个分组条件;
# 标准 Seurat 流程(Normalize→PCA→UMAP)出真实 UMAP 坐标。规模小,CPU 数秒。
if (!file.exists(args$input)) {
  cat("Step 0: 生成合成 Seurat 对象 (synthetic demo only)...\n")
  ng <- 80; nc <- 480
  cts   <- c("Tcell", "Bcell", "Mono", "NK")
  ct    <- rep(cts, length.out = nc)
  genes <- paste0("Gene", sprintf("%02d", seq_len(ng)))
  # 每类前 20 基因里各占 5 个为特异 marker
  marker_idx <- list(Tcell = 1:5, Bcell = 6:10, Mono = 11:15, NK = 16:20)
  m <- matrix(rpois(ng * nc, 2), ng, nc, dimnames = list(genes, paste0("cell", seq_len(nc))))
  for (k in names(marker_idx)) {
    cols <- which(ct == k)
    m[marker_idx[[k]], cols] <- m[marker_idx[[k]], cols] + rpois(length(marker_idx[[k]]) * length(cols), 18)
  }
  so <- CreateSeuratObject(counts = m)
  so$celltype  <- factor(ct, levels = cts)
  so$condition <- factor(sample(c("Healthy", "Disease"), nc, TRUE), levels = c("Healthy", "Disease"))
  so <- NormalizeData(so, verbose = FALSE)
  so <- FindVariableFeatures(so, verbose = FALSE)
  so <- ScaleData(so, features = genes, verbose = FALSE)
  so <- RunPCA(so, npcs = 15, verbose = FALSE)
  so <- RunUMAP(so, dims = 1:10, verbose = FALSE)
  saveRDS(so, args$input)
  cat(sprintf("  合成对象: %d 细胞 × %d 基因, 4 类细胞, 2 条件; UMAP 已算\n", nc, ng))
}

## ---- 1. 读对象 + 取 marker 列表 --------------------------------------------
cat("Step 1: 读 Seurat 对象...\n")
so <- readRDS(args$input)
stopifnot(args$celltype %in% colnames(so@meta.data))
Idents(so) <- args$celltype
ct_levels <- levels(factor(so@meta.data[[args$celltype]]))
ct_cols   <- pal_for(ct_levels)
# marker:每类取该类平均表达最高的 3 个基因(真实数据这里换成 FindAllMarkers 结果)
cat("  按类挑 marker 基因(每类 top-3 by 平均表达)...\n")
expr <- GetAssayData(so, layer = "data")
pick_markers <- function() {
  out <- c()
  for (k in ct_levels) {
    cells <- colnames(so)[so@meta.data[[args$celltype]] == k]
    mu <- Matrix::rowMeans(expr[, cells, drop = FALSE])
    out <- c(out, names(sort(mu, decreasing = TRUE))[1:3])
  }
  unique(out)
}
markers <- pick_markers()
cat(sprintf("  marker 基因 %d 个: %s\n", length(markers), paste(head(markers, 12), collapse = ", ")))

## ---- 2. 诚实基线对照:Seurat 默认 UMAP(简陋默认图)-------------------------
cat("Step 2: 诚实基线 — Seurat 默认 DimPlot(对照,非分析差异,仅审美)...\n")
p_base <- Seurat::DimPlot(so, group.by = args$celltype, cols = ct_cols) +
  ggtitle("Seurat default DimPlot (baseline)")
save_fig(p_base, file.path(ASSETS, "00_baseline_seurat_default_umap"), width = 6, height = 5.2)

## ---- 3. SCpubr UMAP(do_DimPlot)——刊级排版 + 色盲友好 ----------------------
cat("Step 3: SCpubr do_DimPlot(UMAP, 色盲友好)...\n")
p_umap <- SCpubr::do_DimPlot(sample = so, group.by = args$celltype,
                             colors.use = ct_cols, label = TRUE, repel = TRUE,
                             legend.position = "right",
                             plot.title = "Cell types (SCpubr do_DimPlot)")
save_fig(p_umap, file.path(ASSETS, "01_scpubr_umap_celltype"), width = 6.2, height = 5.2)

# split by condition —— 一行展示各条件 UMAP
if (args$condition %in% colnames(so@meta.data)) {
  p_split <- SCpubr::do_DimPlot(sample = so, group.by = args$celltype, split.by = args$condition,
                                colors.use = ct_cols, legend.position = "bottom",
                                plot.title = "Cell types split by condition")
  save_fig(p_split, file.path(ASSETS, "02_scpubr_umap_split_condition"), width = 9, height = 5)
}

## ---- 4. SCpubr DotPlot(marker 点图,z-score)——非条形 ---------------------
cat("Step 4: SCpubr do_DotPlot(marker 点图, z-score)...\n")
p_dot <- SCpubr::do_DotPlot(sample = so, features = markers, group.by = args$celltype,
                            zscore.data = TRUE, legend.position = "right",
                            plot.title = "Marker expression (z-scored)")
save_fig(p_dot, file.path(ASSETS, "03_scpubr_dotplot_markers"), width = 8.5, height = 4.6)

## ---- 5. SCpubr FeaturePlot(精修表达 UMAP)——viridis,非灰红 ---------------
cat("Step 5: SCpubr do_FeaturePlot(精修表达)...\n")
p_feat <- SCpubr::do_FeaturePlot(sample = so, features = markers[1:4], ncol = 2,
                                 order = TRUE,
                                 plot.title = "Marker feature plots")
save_fig(p_feat, file.path(ASSETS, "04_scpubr_featureplot"), width = 8, height = 7)

## ---- 6. SCpubr ViolinPlot(小提琴 + 箱线)——非条形 -------------------------
# 注:多 feature 时 SCpubr 自动按基因分面并以 y 轴基因名区分,不传 plot.title 以免
#     各分面标题相互重叠;legend.position="none"(分面已用颜色+x 轴标明细胞类型)。
cat("Step 6: SCpubr do_ViolinPlot(小提琴, 箱线叠加)...\n")
p_vln <- SCpubr::do_ViolinPlot(sample = so, features = markers[1:3], group.by = args$celltype,
                               colors.use = ct_cols, plot_boxplot = TRUE,
                               legend.position = "none")
save_fig(p_vln, file.path(ASSETS, "05_scpubr_violin"), width = 8.5, height = 4.2)

## ---- 7. SCpubr AlluvialPlot(细胞流向桑基)——非条形 ------------------------
if (args$condition %in% colnames(so@meta.data)) {
  cat("Step 7: SCpubr do_AlluvialPlot(condition → celltype 桑基)...\n")
  p_all <- SCpubr::do_AlluvialPlot(sample = so, first_group = args$condition,
                                   last_group = args$celltype, colors.use = ct_cols,
                                   fill.by = args$celltype, use_geom_flow = TRUE,
                                   plot.title = "Cell composition flow")
  save_fig(p_all, file.path(ASSETS, "06_scpubr_alluvial_flow"), width = 6.5, height = 5.5)
}

## ---- 8. 色盲安检:do_ColorBlindCheck(三型色盲模拟)------------------------
cat("Step 8: SCpubr do_ColorBlindCheck(deutan/protan/tritan 色盲安检)...\n")
p_cb <- SCpubr::do_ColorBlindCheck(colors.use = ct_cols)
save_fig(p_cb, file.path(ASSETS, "07_scpubr_colorblind_check"), width = 7.5, height = 4)

## ---- 收尾:统计落盘 + 依赖快照 ---------------------------------------------
summ <- data.frame(
  item  = c("cells", "genes", "celltypes", "markers_used", "palette", "colorblind_safe_check"),
  value = c(ncol(so), nrow(so), length(ct_levels), length(markers), "npg (color-blind friendly)",
            "see assets/07_scpubr_colorblind_check"))
write.csv(summ, file.path(args$outdir, "figure_summary.csv"), row.names = FALSE)
cat("\n完成。图见 assets/(00=基线对照, 01-07=SCpubr 刊级图);摘要见 results/figure_summary.csv\n")
sink(file.path(args$outdir, "sessionInfo.txt")); print(sessionInfo()); sink()   # 依赖版本快照(铁律6)
