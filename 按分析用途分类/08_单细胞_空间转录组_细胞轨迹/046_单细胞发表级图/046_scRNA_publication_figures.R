# =============================================================================
# 编号       : 046
# 脚本名     : 单细胞发表级图 (turnkey + 顶刊图)
# 分类       : 08_单细胞_空间转录组_细胞轨迹
# 用途       : Seurat 标准流程(QC→归一→PCA→聚类→UMAP→marker),输出发表级
#              UMAP / marker 点图 / marker 热图 / 目标基因 FeaturePlot 与小提琴图。
# 方法/包    : Seurat;主题 theme_pub.R(viridis 连续色 + 期刊离散配色)
# 结果图     : UMAP_clusters;Marker_dotplot;Marker_heatmap;<gene>_FeaturePlot;<gene>_violin
# -----------------------------------------------------------------------------
# 运行(示例): Rscript 046_scRNA_publication_figures.R
# 运行(自己): Rscript 046_scRNA_publication_figures.R --input data/counts.csv --genes "CD3D,MS4A1"
# 可选参数 : --resolution 0.5 --npcs 20 --minfeature 100
# 输入规格 : 计数矩阵 CSV(首列基因名,其余列=细胞;原始 counts)。
# 整理日期 : 2026-06-23(turnkey 重构;Set3→期刊配色,grey-red→viridis)
# =============================================================================

.find_fw <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  d <- if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
  for (i in 0:6) { cand <- file.path(d, paste(rep("..", i), collapse = "/"), "_framework", "theme_pub.R")
    if (file.exists(cand)) return(normalizePath(cand)) }
  stop("未找到 _framework/theme_pub.R")
}
source(.find_fw())
suppressWarnings(suppressMessages({ library(Seurat); library(ggplot2); library(dplyr) }))

SCRIPT_DIR <- bio_script_dir()
args <- bio_args(list(input = file.path(SCRIPT_DIR, "example_data", "counts.csv"), outdir = file.path(SCRIPT_DIR, "results"),
                      genes = "", resolution = "0.5", npcs = "20", minfeature = "100"))
ASSETS <- file.path(SCRIPT_DIR, "assets")
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE); dir.create(ASSETS, recursive = TRUE, showWarnings = FALSE)
PAL <- pal_pub(name = "npg")

cat("Step 1/5: 读取 counts + 建 Seurat 对象...\n")
cnt <- read_table_smart(args$input, row_names = TRUE)
so <- CreateSeuratObject(counts = as.matrix(cnt), min.cells = 3, min.features = as.integer(args$minfeature))
cat("  ", ncol(so), "细胞 x", nrow(so), "基因\n")

cat("Step 2/5: 归一 → HVG → PCA → 聚类 → UMAP...\n")
so <- NormalizeData(so, verbose = FALSE) |> FindVariableFeatures(nfeatures = 2000, verbose = FALSE) |> ScaleData(verbose = FALSE)
so <- RunPCA(so, npcs = as.integer(args$npcs), verbose = FALSE)
so <- FindNeighbors(so, dims = 1:as.integer(args$npcs), verbose = FALSE) |> FindClusters(resolution = as.numeric(args$resolution), verbose = FALSE)
so <- tryCatch(RunUMAP(so, dims = 1:as.integer(args$npcs), verbose = FALSE), error = function(e) { cat("  UMAP 失败,改用 tSNE\n"); RunTSNE(so, dims = 1:as.integer(args$npcs)) })
red <- if ("umap" %in% names(so@reductions)) "umap" else "tsne"
ncl <- length(levels(so)); cols <- pal_pub(ncl, "npg")
cat("  聚类数:", ncl, "\n")

cat("Step 3/5: UMAP 聚类图...\n")
p_umap <- DimPlot(so, reduction = red, label = TRUE, label.size = 5, cols = cols, pt.size = 0.6) +
  labs(title = "Cell clusters (UMAP)") + theme_pub(base_size = 12, border = TRUE) + theme(legend.position = "right")
save_fig(p_umap, file.path(ASSETS, "UMAP_clusters"), 6.5, 5.5); save_fig(p_umap, file.path(args$outdir, "UMAP_clusters"), 6.5, 5.5)

cat("Step 4/5: marker 识别 + 点图/热图...\n")
mk <- FindAllMarkers(so, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, verbose = FALSE)
write.csv(mk, file.path(args$outdir, "cluster_markers.csv"), row.names = FALSE)
top <- mk |> group_by(cluster) |> slice_max(avg_log2FC, n = 3) |> ungroup()
p_dot <- DotPlot(so, features = unique(top$gene), cols = c("lightgrey", "#E64B35")) +
  labs(title = "Top markers", x = NULL, y = "Cluster") + theme_pub(base_size = 10, border = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
save_fig(p_dot, file.path(ASSETS, "Marker_dotplot"), 9, 5); save_fig(p_dot, file.path(args$outdir, "Marker_dotplot"), 9, 5)
top5 <- mk |> group_by(cluster) |> slice_max(avg_log2FC, n = 5) |> ungroup()
p_hm <- DoHeatmap(so, features = unique(top5$gene), group.colors = cols, size = 3) +
  scale_fill_viridis_c(option = "D") + theme(axis.text.y = element_text(size = 6))
save_fig(p_hm, file.path(ASSETS, "Marker_heatmap"), 9, 7); save_fig(p_hm, file.path(args$outdir, "Marker_heatmap"), 9, 7)

cat("Step 5/5: 目标基因 FeaturePlot + 小提琴...\n")
genes <- if (nzchar(args$genes)) trimws(strsplit(args$genes, ",")[[1]]) else unique(top$gene)[1:2]
genes <- intersect(genes, rownames(so))
for (g in genes) {
  pf <- FeaturePlot(so, features = g, reduction = red, pt.size = 0.6) +
    scale_colour_viridis_c(option = "C") + labs(title = paste0(g, " (UMAP)")) + theme_pub(base_size = 12, border = TRUE)
  save_fig(pf, file.path(ASSETS, paste0(g, "_FeaturePlot")), 6, 5.5); save_fig(pf, file.path(args$outdir, paste0(g, "_FeaturePlot")), 6, 5.5)
  pv <- VlnPlot(so, features = g, cols = cols, pt.size = 0) +
    labs(title = paste0(g, " by cluster"), x = "Cluster") + theme_pub(base_size = 12, border = TRUE) + theme(legend.position = "none")
  save_fig(pv, file.path(ASSETS, paste0(g, "_violin")), 6, 4.5); save_fig(pv, file.path(args$outdir, paste0(g, "_violin")), 6, 4.5)
}
saveRDS(so, file.path(args$outdir, "seurat_object.rds"))
cat("完成。单细胞图见", normalizePath(ASSETS), "\n")
