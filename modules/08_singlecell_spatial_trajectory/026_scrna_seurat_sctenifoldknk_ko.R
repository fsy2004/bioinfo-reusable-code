# =============================================================================
# 编号       : R026
# 脚本名     : 单细胞数据分析.R
# 分类       : 08_singlecell_spatial_trajectory
# 项目来源   : 网络毒理学_scTenifoldKnk_单细胞_空间转录组_免疫_诊断
# 用途       : 执行单细胞质控、降维聚类、细胞注释、标志基因和可视化分析。
# 结果图     : 热图；PCA图；条形图/柱状图；气泡图/点图；箱线图；小提琴图；散点图；UMAP/t-SNE降维图；单细胞基因表达/Marker图
# 主要 R 包  : celda; CellChat; celldex; cluster; ClusterGVis; clusterProfiler; clustree; ComplexHeatmap; DoubletFinder; dplyr; ggalluvial; ggplot2; ggrepel; ggsci; jsonlite; limma; magrittr; Matrix; monocle3; NMF; org.Hs.eg.db; patchwork; RColorBrewer; scRNAtoolVis; scTenifoldKnk; Seurat; SingleR; svglite; tidyr
# 整理时间   : 2026-05-10
# =============================================================================
# ===============================================================================
# 单细胞RNA测序综合分析流程 - 整合版
# 包含：质量控制、聚类分析、聚类树、差异表达、细胞注释、富集分析、
#       轨迹分析、目标基因分析等


# ------------------- 1. 预加载R包及依赖 --------------------
# 基础包
library(Seurat)           # 主要的单细胞分析包
library(dplyr)            # 数据整理、管道操作
library(ggplot2)          # 画图
library(magrittr)         # 管道符 %>%
library(RColorBrewer)     # 颜色方案
library(tidyr)            # 数据整理（长宽表互转等）
library(Matrix)           # 稀疏矩阵
library(patchwork)        # 图像拼接
library(svglite)          # 导出svg格式图
library(limma)            # 差异分析
library(NMF)              # 非负矩阵分解相关（部分包依赖）
library(CellChat)         # 细胞通讯分析
library(ggalluvial)       # Alluvial 图
library(celldex)        # SingleR注释参考数据库（手动注释不需要）
library(SingleR)        # 自动细胞类型注释（手动注释不需要）
library(monocle3)         # 轨迹分析 
library(clustree)         # 聚类树图
library(cluster)          # 轮廓系数计算
library(scRNAtoolVis)     # UMAP美化可视化
library(ggrepel)          # 避免标签重叠
library(ClusterGVis)      # 基因聚类可视化
library(org.Hs.eg.db)     # 人类基因注释数据库
library(clusterProfiler)  # 功能富集分析
library(ComplexHeatmap)   # 复杂热图绘制
library(ggsci)            # 科学期刊配色方案

# 双细胞检测和去污染包
library(DoubletFinder)    # 双细胞检测
library(celda)            # decontX去除环境RNA污染

# ------------------- 2. 设置工作目录和参数 --------------------
workDir <- "C:/Users/fsy/Desktop/网络毒理学+单细胞+scTenifoldKnk模拟基因敲除+空间转录组+细胞轨迹+机器学习+分子对接+免疫分析+诊断模型/33.单细胞数据分析"
setwd(workDir)

# 分析参数
logFC_filter       <- 1        # 差异分析logFC阈值
p_adj_filter       <- 0.05     # 差异分析校正p值阈值
min_cells_gene     <- 5        # 至少在5个细胞出现的基因参与分析
min_genes_per_cell <- 200      # 每个细胞至少检测到的基因数
post_filter_cells  <- 300      # 二次过滤时，细胞内基因数大于此值
post_filter_mito   <- 20       # 线粒体比例不得高于此值（%）
n_top_var_features <- 2500     # 选取前2500个高度变异基因
n_pcs              <- 22       # PCA降维时使用的主成分数
n_topmarker_heat   <- 10       # 每簇热图展示top10 marker
cluster_resolution <- 0.3      # 聚类分辨率
neighbor_dims      <- 15       # 聚类与UMAP时用多少PC
qc_vln_width       <- 15       # 质控小提琴图宽度
qc_vln_height      <- 7        # 质控小提琴图高度

# 输出目录设置
output_dir         <- "analysis_results"
input_expr_file    <- "GSE296117_RA_geo.rds"

# 目标基因设置
target_genes <- c("AIM2")

# GSE296117_RA_geo.rds 已经是作者整理好的 Seurat 对象，优先保留原始注释和降维结果
use_author_celltype_annotation <- TRUE
author_celltype_col <- "celltype"
skip_decontX_for_preprocessed_seurat <- TRUE
skip_doubletfinder_for_preprocessed_seurat <- TRUE

# 基因聚类富集分析中的重要标记基因
important_marker_genes <- c("AIM2", "LTF", "ELANE", "MPO", "PGLYRP1", "BPI", "CEACAM8", "RETN", "CD177", "MMP9", "SLPI",
                           "CD3D", "CD3E", "CD8A", "CD4", "CD19", "CD14", "FCGR3A",
                           "LYZ", "MS4A1", "GNLY", "PRF1", "CST3", "FCER1A")

# 基因敲除分析参数
ko_n_hvg           <- 2000      # 用于构建网络的高变基因数量
ko_qc_mtThreshold  <- 0.1       # 线粒体基因比例阈值
ko_qc_minLSize     <- 1000      # 文库大小阈值（细胞测到的基因总数）
ko_nc_nNet         <- 3         # 子网络数量（降低以加速，影响最大）
ko_nc_nCells       <- 300       # 每个网络中随机抽取的细胞数（降低以加速）
ko_pval_threshold  <- 0.05      # 显著性p值阈值

# ------------------- 3. 创建输出文件夹结构 --------------------
# 创建主输出目录
if(!dir.exists(output_dir)) dir.create(output_dir, recursive=TRUE)

# 创建子目录结构
subdirs <- c(
  "01_Quality_Control",         # 质量控制
  "02_Clustering_Analysis",     # 聚类分析
  "03_Cell_Type_Annotation",    # 细胞类型注释
  "04_Differential_Expression", # 差异表达分析
  "05_Gene_Cluster_Enrichment", # 基因聚类富集分析
  "06_Trajectory_Analysis",     # 轨迹分析
  "07_Target_Gene_Analysis",    # 目标基因分析
  "08_Data_Export",             # 数据导出
  "09_Statistics_Plots",        # 统计可视化
  "10_Gene_Knockout_Analysis"   # 基因敲除分析
)

for(subdir in subdirs) {
  dir_path <- file.path(output_dir, subdir)
  if(!dir.exists(dir_path)) dir.create(dir_path, recursive=TRUE)
}

cat("输出目录结构已创建：", output_dir, "\n")


# 第一部分：数据读取与预处理


cat("Step 1: 读取表达矩阵...\n")
counts <- readRDS(input_expr_file)
input_is_seurat <- inherits(counts, "Seurat")
raw_seurat_object <- NULL

# 检查读取的数据结构并提取表达矩阵
cat("读取的数据类型:", class(counts), "\n")

if (inherits(counts, "Seurat")) {
  cat("检测到Seurat对象，提取counts矩阵\n")
  raw_seurat_object <- counts
  if (packageVersion("Seurat") >= "5.0.0") {
    cat("检测到Seurat v5，使用新的API\n")
    expr_matrix <- tryCatch({
      joined_counts <- JoinLayers(counts)
      raw_seurat_object <<- joined_counts
      GetAssayData(joined_counts, assay = "RNA", layer = "counts")
    }, error = function(e) {
      cat("无法合并layers，尝试直接提取counts layer\n")
      tryCatch({
        LayerData(counts, assay = "RNA", layer = "counts")
      }, error = function(e2) {
        cat("LayerData提取失败，尝试GetAssayData直接提取counts layer\n")
        GetAssayData(counts, assay = "RNA", layer = "counts")
      })
    })
  } else {
    cat("检测到Seurat v4或更早版本\n")
    expr_matrix <- GetAssayData(counts, assay = "RNA", slot = "counts")
  }
} else if (inherits(counts, "SingleCellExperiment")) {
  cat("检测到SingleCellExperiment对象，提取counts矩阵\n")
  expr_matrix <- assay(counts, "counts")
} else if (inherits(counts, "dgCMatrix") || inherits(counts, "dgTMatrix") || 
           inherits(counts, "Matrix")) {
  cat("检测到稀疏矩阵格式\n")
  expr_matrix <- counts
} else if (is.list(counts)) {
  cat("数据是列表格式，包含以下元素:", names(counts), "\n")
  if ("RNA" %in% names(counts)) {
    expr_matrix <- counts$RNA
  } else if ("counts" %in% names(counts)) {
    expr_matrix <- counts$counts
  } else if (length(counts) == 1) {
    expr_matrix <- counts[[1]]
  } else {
    stop("无法确定表达矩阵的位置，请检查数据结构")
  }
} else if (is.matrix(counts) || is.data.frame(counts)) {
  cat("检测到矩阵或数据框格式\n")
  expr_matrix <- counts
} else {
  expr_matrix <- counts
}

cat("表达矩阵维度：", nrow(expr_matrix), "行（基因）", ncol(expr_matrix), "列（细胞）\n")


# 第二部分：质量控制


cat("Step 2: 创建Seurat对象与质量控制...\n")

# 处理矩阵格式并检查NA值
if (any(is.na(expr_matrix))) {
  cat("警告：数据中包含NA值，将替换为0\n")
  expr_matrix[is.na(expr_matrix)] <- 0
}

# 创建Seurat对象；如果输入本来就是Seurat对象，则保留作者metadata、celltype和UMAP
if (!is.null(raw_seurat_object)) {
  cat("使用输入Seurat对象作为分析对象，保留作者metadata和已有降维结果\n")
  scObject <- raw_seurat_object
  DefaultAssay(scObject) <- "RNA"
} else {
  scObject <- CreateSeuratObject(
    counts = expr_matrix,
    project = "scRNAseqProject",
    min.cells = min_cells_gene,
    min.features = min_genes_per_cell
  )
}

# 统一线粒体比例列名，兼容percent.mt和percent.mito
if ("percent.mt" %in% colnames(scObject@meta.data)) {
  scObject$percent.mito <- scObject$percent.mt
} else if (!"percent.mito" %in% colnames(scObject@meta.data)) {
  scObject[["percent.mito"]] <- PercentageFeatureSet(scObject, pattern = "^MT-")
}

# 保存过滤前的对象用于对比
scObject_before_filter <- scObject
cell_num_before <- ncol(scObject)

# QC小提琴图（过滤前）
pdf(file.path(output_dir, "01_Quality_Control", "QC_violin_before_filtering.pdf"), width=12, height=6)
p_before <- VlnPlot(
  scObject_before_filter,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mito"),
  pt.size = 0,
  group.by = "orig.ident"
) +
  ggtitle(paste0("Before Filtering (n = ", cell_num_before, " cells)")) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.position = "right",
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.text.x  = element_text(size = 12, face = "bold", angle = 45, vjust=1, hjust=1),
    axis.text.y  = element_text(size = 12, face = "bold"),
    strip.text   = element_text(size = 15, face = "bold", color = "black")
  )
print(p_before)
dev.off()

# 二次过滤与散点图
cat("Step 3: 二次细胞质控...\n")
scObject <- subset(scObject, subset = nFeature_RNA > post_filter_cells & percent.mito < post_filter_mito)
cell_num_after <- ncol(scObject)
cat(sprintf("细胞过滤: 原%d，剩%d。\n", cell_num_before, cell_num_after))

if(cell_num_after < 10) stop("过滤后过少细胞，检查阈值！")

# QC小提琴图（过滤后）
pdf(file.path(output_dir, "01_Quality_Control", "QC_violin_after_filtering.pdf"), width=12, height=6)
p_after <- VlnPlot(
  scObject,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mito"),
  pt.size = 0,
  group.by = "orig.ident"
) +
  ggtitle(paste0("After Filtering (n = ", cell_num_after, " cells)")) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.position = "right",
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.text.x  = element_text(size = 12, face = "bold", angle = 45, vjust=1, hjust=1),
    axis.text.y  = element_text(size = 12, face = "bold"),
    strip.text   = element_text(size = 15, face = "bold", color = "black")
  )
print(p_after)
dev.off()

# QC小提琴图对比（过滤前后合并到一张图）
pdf(file.path(output_dir, "01_Quality_Control", "QC_violin_basicMetrics.pdf"), width=14, height=12)

# 重新创建过滤前的小提琴图用于对比
p_before_comp <- VlnPlot(
  scObject_before_filter,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mito"),
  pt.size = 0,
  group.by = "orig.ident"
) +
  ggtitle(paste0("Before Filtering (n = ", cell_num_before, " cells)")) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", color = "#E74C3C"),
    legend.position = "none",
    axis.title.x = element_text(size = 13, face = "bold"),
    axis.title.y = element_text(size = 13, face = "bold"),
    axis.text.x  = element_text(size = 11, face = "bold", angle = 45, vjust=1, hjust=1),
    axis.text.y  = element_text(size = 11, face = "bold"),
    strip.text   = element_text(size = 13, face = "bold", color = "black")
  )

p_after_comp <- VlnPlot(
  scObject,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mito"),
  pt.size = 0,
  group.by = "orig.ident"
) +
  ggtitle(paste0("After Filtering (n = ", cell_num_after, " cells)")) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", color = "#27AE60"),
    legend.position = "none",
    axis.title.x = element_text(size = 13, face = "bold"),
    axis.title.y = element_text(size = 13, face = "bold"),
    axis.text.x  = element_text(size = 11, face = "bold", angle = 45, vjust=1, hjust=1),
    axis.text.y  = element_text(size = 11, face = "bold"),
    strip.text   = element_text(size = 13, face = "bold", color = "black")
  )

# 合并两个图（上下排列）使用patchwork
comparison_plot <- p_before_comp / p_after_comp +
  patchwork::plot_annotation(
    title = "Quality Control: Before vs After Filtering",
    subtitle = sprintf("Filtered %d cells (%.1f%% removed)",
                      cell_num_before - cell_num_after,
                      (cell_num_before - cell_num_after) / cell_num_before * 100),
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 14, color = "gray40")
    )
  )

print(comparison_plot)
dev.off()

cat(sprintf("QC小提琴图已保存（过滤前: %d细胞 → 过滤后: %d细胞，移除: %.1f%%）\n",
            cell_num_before, cell_num_after,
            (cell_num_before - cell_num_after) / cell_num_before * 100))

pdf(file.path(output_dir, "01_Quality_Control", "QC_scatter_metrics.pdf"), width=13, height=7)
plot1 <- FeatureScatter(scObject, feature1 = "nCount_RNA", feature2 = "percent.mito", pt.size = 1.5)
plot2 <- FeatureScatter(scObject, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 1.5)
print(plot1 + plot2)
dev.off()

# ===============================================================================
# 第2.5部分：decontX去除环境RNA污染
# ===============================================================================

if (!(input_is_seurat && skip_decontX_for_preprocessed_seurat)) {

cat("Step 3.5: 使用decontX去除环境RNA污染...\n")

# 获取原始count矩阵
raw_counts <- GetAssayData(scObject, layer = "counts")

# 运行decontX
decontX_results <- decontX(raw_counts)

# 获取校正后的矩阵和污染比例
decontXcounts <- decontX_results$decontXcounts
contamination <- decontX_results$contamination

# 将污染比例添加到metadata
scObject$decontX_contamination <- contamination

# 统计污染情况
cat(sprintf("  平均污染比例: %.2f%%\n", mean(contamination) * 100))
cat(sprintf("  中位污染比例: %.2f%%\n", median(contamination) * 100))
cat(sprintf("  最大污染比例: %.2f%%\n", max(contamination) * 100))

# 可视化污染分布
pdf(file.path(output_dir, "01_Quality_Control", "decontX_Contamination_Distribution.pdf"), width = 8, height = 6)
hist(contamination, breaks = 50, main = "decontX Contamination Distribution",
     xlab = "Contamination Proportion", col = "steelblue", border = "white")
abline(v = mean(contamination), col = "red", lwd = 2, lty = 2)
legend("topright", legend = paste("Mean:", round(mean(contamination), 3)), col = "red", lty = 2, lwd = 2)
dev.off()

# 用校正后的矩阵替换原始数据
scObject[["RNA"]]$counts <- decontXcounts

# 保存decontX结果摘要
decontX_summary <- data.frame(
  Metric = c("Mean_Contamination", "Median_Contamination", "Max_Contamination", "Min_Contamination"),
  Value = c(round(mean(contamination), 4), round(median(contamination), 4),
            round(max(contamination), 4), round(min(contamination), 4))
)
write.csv(decontX_summary, file.path(output_dir, "01_Quality_Control", "decontX_Summary.csv"), row.names = FALSE)

# 过滤高污染细胞 (污染比例 > 50%)
high_contamination_threshold <- 0.5
cells_before_contam_filter <- ncol(scObject)

scObject <- subset(
  scObject,
  subset = decontX_contamination < high_contamination_threshold
)

cells_after_contam_filter <- ncol(scObject)
cat(sprintf("  高污染细胞过滤: 去除 %d 个细胞 (污染>50%%)\n",
            cells_before_contam_filter - cells_after_contam_filter))
cat(sprintf("  过滤后剩余: %d 个细胞\n", cells_after_contam_filter))

cat("decontX环境RNA去除完成！\n")

} else {
  cat("Step 3.5: 检测到预处理Seurat对象，跳过decontX环境RNA去除。\n")
}

# ===============================================================================
# 第2.6部分：DoubletFinder双细胞检测和去除
# ===============================================================================

if (!(input_is_seurat && skip_doubletfinder_for_preprocessed_seurat)) {

cat("Step 3.6: 使用DoubletFinder进行双细胞检测和去除...\n")

# 记录去除前细胞数
cells_before_doublet <- ncol(scObject)

# DoubletFinder需要预处理数据
scObject <- NormalizeData(scObject)
scObject <- FindVariableFeatures(scObject, selection.method = "vst", nfeatures = 2000)
scObject <- ScaleData(scObject)
scObject <- RunPCA(scObject)
scObject <- RunUMAP(scObject, dims = 1:10)
scObject <- FindNeighbors(scObject, dims = 1:10)
scObject <- FindClusters(scObject, resolution = 0.5)

# 参数扫描确定最优pK
sweep.res.list <- paramSweep(scObject, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK_value <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
cat(sprintf("  最优pK值: %s\n", pK_value))

# 根据细胞数自动估算双细胞率 (每1000细胞约0.8%)
nCells <- ncol(scObject)
doublet_rate <- nCells * 0.8 / 100000
cat(sprintf("  细胞数: %d, 预估双细胞率: %.2f%%\n", nCells, doublet_rate * 100))

homotypic.prop <- modelHomotypic(scObject$seurat_clusters)
nExp_poi <- round(doublet_rate * nCells)
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

# 运行DoubletFinder
scObject <- doubletFinder(scObject, PCs = 1:10, pN = 0.25, pK = pK_value, nExp = nExp_poi.adj, sct = FALSE)

# 获取分类列名并去除双细胞
df_col <- grep("^DF.classifications", colnames(scObject@meta.data), value = TRUE)[1]
scObject$Doublet_Classification <- scObject@meta.data[[df_col]]

# 可视化双细胞分布
pdf(file.path(output_dir, "01_Quality_Control", "DoubletFinder_UMAP.pdf"), width = 10, height = 8)
print(DimPlot(scObject, group.by = "Doublet_Classification", cols = c("Singlet" = "grey", "Doublet" = "red")))
dev.off()

# 去除双细胞
scObject <- subset(scObject, subset = Doublet_Classification == "Singlet")

# 统计结果
cells_after_doublet <- ncol(scObject)
cat(sprintf("  双细胞去除前: %d 细胞\n", cells_before_doublet))
cat(sprintf("  双细胞去除后: %d 细胞\n", cells_after_doublet))
cat(sprintf("  去除双细胞: %d (%.2f%%)\n", cells_before_doublet - cells_after_doublet,
            (cells_before_doublet - cells_after_doublet)/cells_before_doublet*100))

# 保存双细胞检测结果摘要
doublet_summary <- data.frame(
  Metric = c("Cells_Before_Doublet_Removal", "Cells_After_Doublet_Removal",
             "Doublets_Removed", "Doublet_Rate_Percent", "Optimal_pK"),
  Value = c(cells_before_doublet, cells_after_doublet,
            cells_before_doublet - cells_after_doublet,
            round((cells_before_doublet - cells_after_doublet)/cells_before_doublet*100, 2),
            pK_value)
)
write.csv(doublet_summary, file.path(output_dir, "01_Quality_Control", "DoubletFinder_Summary.csv"), row.names = FALSE)

# 重置对象用于后续标准分析
scObject@reductions <- list()
scObject@graphs <- list()
scObject$seurat_clusters <- NULL

cat("DoubletFinder双细胞检测和去除完成！\n")

} else {
  cat("Step 3.6: 检测到预处理Seurat对象，跳过DoubletFinder双细胞检测。\n")
}

# ===============================================================================
# 第三部分：数据标准化与降维
# ===============================================================================

cat("Step 4: 归一化+高变基因...\n")
scObject <- NormalizeData(scObject, normalization.method="LogNormalize", scale.factor=10000)
scObject <- FindVariableFeatures(scObject, selection.method="vst", nfeatures=n_top_var_features)

pdf(file.path(output_dir, "01_Quality_Control", "VarGenes_overview.pdf"), width=10, height=6)
VariableFeaturePlot(scObject)
dev.off()

cat("Step 5: PCA降维...\n")
scObject <- ScaleData(scObject)
scObject <- RunPCA(scObject, npcs=n_pcs, features=VariableFeatures(scObject))

pdf(file.path(output_dir, "02_Clustering_Analysis", "PCA_geneLoadings.pdf"))
VizDimLoadings(scObject, dims=1:3, reduction="pca", nfeatures=25)
dev.off()

pdf(file.path(output_dir, "02_Clustering_Analysis", "PCA_heatmap_topGenes.pdf"), width=8, height=7)
DimHeatmap(scObject, dims=1:3, cells=400, balanced=TRUE)
dev.off()

# ===============================================================================
# 第四部分：聚类分析（包含聚类树）
# ===============================================================================

cat("Step 6: 聚类与UMAP降维...\n")
scObject <- FindNeighbors(scObject, dims=1:neighbor_dims)

# 多分辨率聚类用于clustree分析
resolutions <- seq(0, 1, 0.1)  # 从0到1，步长0.1的分辨率
for (res in resolutions) {
  scObject <- FindClusters(scObject, resolution = res)
}

# ----------- Clustree图 - 展示不同分辨率下的聚类关系 ----------
cat("Step 7: 生成Clustree图...\n")
pdf(file.path(output_dir, "02_Clustering_Analysis", "Clustree_resolution_comparison.pdf"), width=15, height=10)
clustree_plot <- clustree(scObject@meta.data, prefix = "RNA_snn_res.") +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.key.size = unit(0.3, "cm"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.border = element_blank(),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.margin = unit(c(1, 4, 1, 1), "cm")
  ) +
  labs(title = "Clustering Tree Across Different Resolutions") +
  guides(
    colour = guide_legend(ncol = 2, title = "Cluster", override.aes = list(size = 2)),
    size = guide_legend(ncol = 2, title = "Sample count", override.aes = list(alpha = 0.7))
  )
print(clustree_plot)
dev.off()

# ----------- 自动分析不同分辨率特征 ----------
cat("Step 8: 分析不同分辨率特征...\n")

# 创建分辨率分析结果表
resolution_analysis <- data.frame(
  Resolution = numeric(),
  N_Clusters = integer(),
  Min_Cluster_Size = integer(),
  Max_Cluster_Size = integer(),
  Mean_Cluster_Size = numeric(),
  Small_Clusters_Count = integer(),
  Silhouette_Score = numeric(),
  stringsAsFactors = FALSE
)

# 分析每个分辨率
for (res in resolutions) {
  col_name <- paste0("RNA_snn_res.", res)
  if (col_name %in% colnames(scObject@meta.data)) {
    clusters <- scObject@meta.data[[col_name]]
    cluster_sizes <- table(clusters)
    
    # 计算基本统计
    n_clusters <- length(unique(clusters))
    min_size <- min(cluster_sizes)
    max_size <- max(cluster_sizes)
    mean_size <- mean(cluster_sizes)
    small_clusters <- sum(cluster_sizes < 50)  # 小于50个细胞的聚类数量
    
    # 计算轮廓系数（可选，计算较慢）
    if (n_clusters > 1 && n_clusters < 20) {  # 避免过多聚类导致计算过慢
      tryCatch({
        # 使用PCA降维后的数据计算轮廓系数
        pca_data <- Embeddings(scObject, reduction = "pca")[, 1:10]
        sil_score <- mean(cluster::silhouette(as.numeric(clusters), dist(pca_data))[, 3])
      }, error = function(e) {
        sil_score <- NA
      })
    } else {
      sil_score <- NA
    }
    
    # 添加到结果表
    resolution_analysis <- rbind(resolution_analysis, data.frame(
      Resolution = res,
      N_Clusters = n_clusters,
      Min_Cluster_Size = min_size,
      Max_Cluster_Size = max_size,
      Mean_Cluster_Size = round(mean_size, 1),
      Small_Clusters_Count = small_clusters,
      Silhouette_Score = round(sil_score, 3),
      stringsAsFactors = FALSE
    ))
    
    cat(sprintf("分辨率 %.1f: %d个聚类 (最小:%d, 最大:%d, 平均:%.1f, 小聚类:%d个)\n", 
                res, n_clusters, min_size, max_size, mean_size, small_clusters))
  }
}

# 保存分析结果
write.csv(resolution_analysis, file.path(output_dir, "02_Clustering_Analysis", "Resolution_Analysis_Summary.csv"), row.names = FALSE)

cat("详细分析结果已保存至: Resolution_Analysis_Summary.csv\n")

# 使用默认分辨率进行最终聚类
scObject <- FindClusters(scObject, resolution=cluster_resolution) 
scObject <- RunUMAP(scObject, dims=1:neighbor_dims)

nColors <- length(unique(scObject$seurat_clusters))
myColors <- colorRampPalette(brewer.pal(12, "Set3"))(nColors)

pdf(file.path(output_dir, "02_Clustering_Analysis", "UMAP_clustered_samples.pdf"), width=7, height=5)
DimPlot(scObject, reduction="umap", label=TRUE, cols=myColors)
dev.off()

write.csv(data.frame(Cell=colnames(scObject), Cluster=scObject$seurat_clusters),
          file=file.path(output_dir, "08_Data_Export", "CellCluster_UMAP_assignments.csv"), row.names=FALSE)

# ===============================================================================
# 第四部分补充：高质量UMAP可视化（Nature风格）
# Publication-quality UMAP visualization
# ===============================================================================

cat("Step 6.5: 生成高质量UMAP可视化图（Nature风格）...\n")

# 创建UMAP可视化专用目录
umap_viz_dir <- file.path(output_dir, "02_Clustering_Analysis", "UMAP_Publication_Style")
if (!dir.exists(umap_viz_dir)) dir.create(umap_viz_dir, recursive = TRUE)

# 提取UMAP坐标
umap_coords <- as.data.frame(Embeddings(scObject, "umap"))
colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
umap_coords$Cluster <- as.character(scObject$seurat_clusters)
umap_coords$Cell <- rownames(umap_coords)

# 获取细胞总数
total_cells <- ncol(scObject)
total_cells_label <- format(total_cells, big.mark = ",")

# 计算聚类中心用于标签
cluster_centers <- umap_coords %>%
  group_by(Cluster) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2),
    n_cells = n(),
    .groups = "drop"
  )

# 设置聚类颜色（柔和的颜色方案）
n_clusters <- length(unique(umap_coords$Cluster))

# 使用柔和的颜色方案（类似图片风格）
soft_colors <- c(
  "#8FBCDB", "#F4A582", "#92C5DE", "#D6604D", "#B2DF8A",
  "#CAB2D6", "#FDBF6F", "#A6CEE3", "#FB9A99", "#E31A1C",
  "#33A02C", "#FF7F00", "#6A3D9A", "#B15928", "#1F78B4",
  "#FFFF99", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD"
)
cluster_colors <- setNames(soft_colors[1:n_clusters], sort(unique(umap_coords$Cluster)))

# ----------- 图1: 聚类UMAP（Nature风格）-----------
cat("  生成聚类UMAP图...\n")

p_cluster <- ggplot(umap_coords, aes(x = UMAP_1, y = UMAP_2, color = Cluster)) +
  geom_point(size = 0.1, alpha = 0.6) +
  scale_color_manual(values = cluster_colors) +
  geom_text(data = cluster_centers,
            aes(x = UMAP_1, y = UMAP_2, label = Cluster),
            color = "black", size = 4, fontface = "bold") +
  labs(
    title = NULL,
    x = "UMAP_1",
    y = "UMAP_2",
    caption = paste0(total_cells_label, " cells")
  ) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.caption = element_text(hjust = 0.5, size = 12, color = "gray30",
                                margin = margin(t = 10)),
    plot.margin = margin(10, 10, 10, 10),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  coord_fixed(ratio = 1)

# 保存聚类UMAP
pdf(file.path(umap_viz_dir, "UMAP_Cluster_NatureStyle.pdf"), width = 8, height = 8)
print(p_cluster)
dev.off()

png(file.path(umap_viz_dir, "UMAP_Cluster_NatureStyle.png"), width = 2400, height = 2400, res = 300)
print(p_cluster)
dev.off()

cat("  聚类UMAP图已保存\n")

# ----------- 图2: 细胞类型UMAP（注释后生成）-----------
# 注意：此图需要在细胞类型注释完成后生成
# 在注释完成后会自动调用下面的函数

# 定义生成细胞类型UMAP的函数（供后续调用）
generate_celltype_umap <- function(seurat_obj, output_dir, celltype_col = "cell_type") {

  cat("  生成细胞类型UMAP图（Nature风格）...\n")

  # 提取UMAP坐标
  umap_coords <- as.data.frame(Embeddings(seurat_obj, "umap"))
  colnames(umap_coords) <- c("UMAP_1", "UMAP_2")

  # 获取细胞类型
  if (celltype_col %in% colnames(seurat_obj@meta.data)) {
    umap_coords$CellType <- as.character(seurat_obj@meta.data[[celltype_col]])
  } else {
    umap_coords$CellType <- as.character(Idents(seurat_obj))
  }

  # 获取细胞总数
  total_cells <- ncol(seurat_obj)
  total_cells_label <- format(total_cells, big.mark = ",")

  # 计算细胞类型中心用于标签
  celltype_centers <- umap_coords %>%
    group_by(CellType) %>%
    summarise(
      UMAP_1 = median(UMAP_1),
      UMAP_2 = median(UMAP_2),
      n_cells = n(),
      .groups = "drop"
    )

  # 设置细胞类型颜色（柔和的颜色方案）
  n_celltypes <- length(unique(umap_coords$CellType))

  # 使用柔和的颜色方案
  soft_colors <- c(
    "#8FBCDB", "#F4A582", "#92C5DE", "#D6604D", "#B2DF8A",
    "#CAB2D6", "#FDBF6F", "#A6CEE3", "#FB9A99", "#E31A1C",
    "#33A02C", "#FF7F00", "#6A3D9A", "#B15928", "#1F78B4",
    "#FFFF99", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD",
    "#80B1D3", "#FDB462", "#BEBADA", "#FB8072", "#8DD3C7"
  )
  celltype_colors <- setNames(soft_colors[1:n_celltypes], sort(unique(umap_coords$CellType)))

  # 生成细胞类型UMAP图
  p_celltype <- ggplot(umap_coords, aes(x = UMAP_1, y = UMAP_2, color = CellType)) +
    geom_point(size = 0.1, alpha = 0.6) +
    scale_color_manual(values = celltype_colors) +
    ggrepel::geom_text_repel(
      data = celltype_centers,
      aes(x = UMAP_1, y = UMAP_2, label = CellType),
      color = "black",
      size = 3.5,
      fontface = "bold",
      box.padding = 0.5,
      point.padding = 0.3,
      segment.color = "gray50",
      segment.size = 0.3,
      max.overlaps = 30,
      force = 2
    ) +
    labs(
      title = NULL,
      x = "UMAP_1",
      y = "UMAP_2",
      caption = paste0(total_cells_label, " cells")
    ) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.caption = element_text(hjust = 0.5, size = 12, color = "gray30",
                                  margin = margin(t = 10)),
      plot.margin = margin(10, 10, 10, 10),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    ) +
    coord_fixed(ratio = 1)

  # 保存细胞类型UMAP
  pdf(file.path(output_dir, "UMAP_CellType_NatureStyle.pdf"), width = 8, height = 8)
  print(p_celltype)
  dev.off()

  png(file.path(output_dir, "UMAP_CellType_NatureStyle.png"), width = 2400, height = 2400, res = 300)
  print(p_celltype)
  dev.off()

  cat("  细胞类型UMAP图已保存\n")

  # ----------- 图3: 聚类和细胞类型并排对比图 -----------
  cat("  生成聚类与细胞类型对比图...\n")

  # 获取聚类信息
  umap_coords$Cluster <- as.character(seurat_obj$seurat_clusters)

  cluster_centers <- umap_coords %>%
    group_by(Cluster) %>%
    summarise(
      UMAP_1 = median(UMAP_1),
      UMAP_2 = median(UMAP_2),
      .groups = "drop"
    )

  n_clusters <- length(unique(umap_coords$Cluster))
  cluster_colors <- setNames(soft_colors[1:n_clusters], sort(unique(umap_coords$Cluster)))

  # 聚类图
  p_cluster_compare <- ggplot(umap_coords, aes(x = UMAP_1, y = UMAP_2, color = Cluster)) +
    geom_point(size = 0.1, alpha = 0.6) +
    scale_color_manual(values = cluster_colors) +
    geom_text(data = cluster_centers,
              aes(x = UMAP_1, y = UMAP_2, label = Cluster),
              color = "black", size = 3.5, fontface = "bold") +
    labs(title = "Clusters", caption = paste0(total_cells_label, " cells")) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.caption = element_text(hjust = 0.5, size = 10, color = "gray30"),
      plot.margin = margin(5, 5, 5, 5)
    ) +
    coord_fixed(ratio = 1)

  # 细胞类型图
  p_celltype_compare <- ggplot(umap_coords, aes(x = UMAP_1, y = UMAP_2, color = CellType)) +
    geom_point(size = 0.1, alpha = 0.6) +
    scale_color_manual(values = celltype_colors) +
    ggrepel::geom_text_repel(
      data = celltype_centers,
      aes(x = UMAP_1, y = UMAP_2, label = CellType),
      color = "black", size = 3, fontface = "bold",
      box.padding = 0.3, max.overlaps = 30
    ) +
    labs(title = "Cell Types", caption = paste0(total_cells_label, " cells")) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.caption = element_text(hjust = 0.5, size = 10, color = "gray30"),
      plot.margin = margin(5, 5, 5, 5)
    ) +
    coord_fixed(ratio = 1)

  # 合并两图
  combined_plot <- p_cluster_compare + p_celltype_compare +
    patchwork::plot_layout(ncol = 2)

  pdf(file.path(output_dir, "UMAP_Cluster_vs_CellType_NatureStyle.pdf"), width = 16, height = 8)
  print(combined_plot)
  dev.off()

  png(file.path(output_dir, "UMAP_Cluster_vs_CellType_NatureStyle.png"), width = 4800, height = 2400, res = 300)
  print(combined_plot)
  dev.off()

  cat("  聚类与细胞类型对比图已保存\n")

  return(list(p_cluster = p_cluster_compare, p_celltype = p_celltype_compare))
}

cat("高质量UMAP可视化函数已定义，将在细胞类型注释完成后自动生成细胞类型图\n")

# ===============================================================================
# 第五部分：差异表达分析
# ===============================================================================

cat("Step 9: 差异表达分析...\n")
marker_df <- tryCatch({
  FindAllMarkers(scObject, only.pos=TRUE, min.pct=0.2, logfc.threshold=logFC_filter)
}, error=function(e){cat("Marker分析异常\n"); NULL})

if(!is.null(marker_df) && nrow(marker_df)>0){
  sigMarkers <- marker_df[abs(marker_df$avg_log2FC)>logFC_filter & marker_df$p_val_adj<p_adj_filter,]
  write.csv(sigMarkers, file=file.path(output_dir, "04_Differential_Expression", "Cluster_Markers_DEGs.csv"), row.names=FALSE)
  
  topmarker <- marker_df %>% group_by(cluster) %>% top_n(n_topmarker_heat, avg_log2FC)
  valid_genes <- topmarker$gene[topmarker$gene %in% rownames(scObject)]
  if(length(valid_genes) > 0){
    pdf(file.path(output_dir, "04_Differential_Expression", "Markers_DoHeatmap.pdf"))
    print(DoHeatmap(scObject, features=valid_genes, size=4) + NoLegend())
    dev.off()
  }
  
  # 保存每个聚类的marker
  split_markers <- split(marker_df, marker_df$cluster)
  for (cid in names(split_markers)) {
    write.csv(split_markers[[cid]], 
              file=file.path(output_dir, "04_Differential_Expression", paste0("Markers_Cluster_", cid, ".csv")), 
              row.names=FALSE)
  }
} else {
  cat('未检出显著marker。\n')
}

# ===============================================================================
# 第六部分：细胞类型自动注释 (SingleR)
# Automatic Cell Type Annotation using SingleR
# ===============================================================================

cat("Step 10: 细胞类型自动注释 (SingleR)...\n")

# 创建注释分析目录
anno_dir <- file.path(output_dir, "03_Cell_Type_Annotation")
if (!dir.exists(anno_dir)) dir.create(anno_dir, recursive = TRUE)

# =============================================
# 10.1 使用SingleR进行自动注释
# =============================================

cat("进行细胞类型注释：优先使用作者celltype，缺失时再使用SingleR...\n")

# 初始化变量
singler_success <- FALSE
cluster_annotation_table <- NULL

# GSE296117对象已包含作者注释celltype，优先使用它，避免SingleR粗略注释覆盖原结果
if (use_author_celltype_annotation && author_celltype_col %in% colnames(scObject@meta.data)) {
  cat(sprintf("检测到作者细胞类型注释列：%s，跳过SingleR并直接使用作者注释。\n", author_celltype_col))

  author_labels <- as.character(scObject@meta.data[[author_celltype_col]])
  author_labels[is.na(author_labels) | author_labels == ""] <- "Unknown"
  names(author_labels) <- colnames(scObject)

  if ("seurat_clusters" %in% colnames(scObject@meta.data)) {
    cluster_annotation_table <- data.frame(
      Cluster = sort(unique(as.character(scObject$seurat_clusters))),
      CellType = NA_character_,
      stringsAsFactors = FALSE
    )

    for (i in seq_len(nrow(cluster_annotation_table))) {
      cl <- cluster_annotation_table$Cluster[i]
      labels_in_cluster <- author_labels[as.character(scObject$seurat_clusters) == cl]
      cluster_annotation_table$CellType[i] <- names(sort(table(labels_in_cluster), decreasing = TRUE))[1]
    }
  } else {
    cluster_annotation_table <- data.frame(
      Cluster = sort(unique(author_labels)),
      CellType = sort(unique(author_labels)),
      stringsAsFactors = FALSE
    )
  }

  singler_result <- list(
    success = TRUE,
    cell_types = author_labels,
    cluster_table = cluster_annotation_table
  )
  singler_success <- TRUE
}

# 检查SingleR和celldex是否可用
if (!singler_success && (!requireNamespace("SingleR", quietly = TRUE) || !requireNamespace("celldex", quietly = TRUE))) {
  cat("警告：SingleR 或 celldex 包未安装，使用默认注释\n")
  singler_success <- FALSE
} else if (!singler_success) {
  # 尝试运行SingleR注释
  singler_result <- tryCatch({
    # 步骤1：提取表达数据和聚类信息
    cat("提取表达数据和聚类信息...\n")

    counts_data <- GetAssayData(object = scObject, layer = "data")
    clusters <- scObject@meta.data$seurat_clusters

    cat(sprintf("表达矩阵维度: %d 基因 x %d 细胞\n", nrow(counts_data), ncol(counts_data)))
    cat(sprintf("聚类数量: %d\n", length(unique(clusters))))

    # 步骤2：加载参考数据集
    cat("Loading SingleR reference dataset...\n")
    celldex_exports <- getNamespaceExports("celldex")
    reference_candidates <- c("MonacoImmuneData", "BlueprintEncodeData", "HumanPrimaryCellAtlasData")
    ref <- NULL
    ref_name <- NA_character_

    for (ref_fun in reference_candidates) {
      if (!(ref_fun %in% celldex_exports)) next

      ref_try <- tryCatch({
        do.call(getExportedValue("celldex", ref_fun), list())
      }, error = function(e) {
        cat(sprintf("SingleR reference %s failed: %s\n", ref_fun, e$message))
        NULL
      })

      if (!is.null(ref_try)) {
        ref <- ref_try
        ref_name <- ref_fun
        break
      }
    }

    if (is.null(ref)) {
      stop("No celldex reference dataset could be loaded")
    }

    ref_label_cols <- colnames(SummarizedExperiment::colData(ref))
    if ("label.fine" %in% ref_label_cols) {
      ref_labels <- ref$label.fine
      ref_label_name <- "label.fine"
    } else if ("label.main" %in% ref_label_cols) {
      ref_labels <- ref$label.main
      ref_label_name <- "label.main"
    } else {
      stop("SingleR参考数据集中未找到label.fine或label.main列")
    }

    cat(sprintf("SingleR reference: %s, labels: %s\n", ref_name, ref_label_name))

    cat(sprintf("参考数据集已加载: %d 基因 x %d 样本\n", nrow(ref), ncol(ref)))

    # 步骤3：进行聚类水平注释（更稳健）
    cat("运行SingleR聚类水平注释...\n")

    singler_clusters <- SingleR(
      test = counts_data,
      ref = ref,
      labels = ref_labels,
      clusters = clusters
    )

    cluster_labels <- singler_clusters$labels
    if ("pruned.labels" %in% colnames(singler_clusters)) {
      pruned_labels <- singler_clusters$pruned.labels
      cluster_labels <- ifelse(is.na(pruned_labels), cluster_labels, pruned_labels)
    }

    # 保存聚类注释结果
    clusterAnn <- data.frame(
      Cluster = rownames(singler_clusters),
      CellType = cluster_labels,
      stringsAsFactors = FALSE
    )

    write.csv(clusterAnn,
              file = file.path(anno_dir, "Cluster_to_CellType_Mapping.csv"),
              row.names = FALSE)

    cat("聚类注释完成。结果:\n")
    for(i in 1:nrow(clusterAnn)) {
      cat(sprintf("  Cluster %s -> %s\n", clusterAnn$Cluster[i], clusterAnn$CellType[i]))
    }

    # 步骤4：将聚类注释映射到单个细胞
    cat("将聚类注释映射到单个细胞...\n")

    # 创建聚类到细胞类型的映射
    cluster_to_celltype <- setNames(clusterAnn$CellType, clusterAnn$Cluster)

    # 获取每个细胞的聚类信息
    cell_clusters <- scObject@meta.data$seurat_clusters

    # 将聚类注释映射到每个细胞
    cell_types_from_clusters <- cluster_to_celltype[as.character(cell_clusters)]
    names(cell_types_from_clusters) <- colnames(scObject)

    # 检查是否有缺失注释
    if(any(is.na(cell_types_from_clusters))) {
      cat("警告：部分细胞缺失注释，使用聚类名称作为备选\n")
      cell_types_from_clusters[is.na(cell_types_from_clusters)] <- paste0("Cluster_",
        cell_clusters[is.na(cell_types_from_clusters)])
    }

    # 返回结果列表
    list(
      success = TRUE,
      cell_types = cell_types_from_clusters,
      cluster_table = clusterAnn
    )

  }, error = function(e) {
    cat(sprintf("错误：SingleR注释失败。错误信息: %s\n", e$message))
    list(success = FALSE, error = e$message)
  })

  singler_success <- singler_result$success
}

# 根据SingleR结果或使用默认注释来更新scObject
if (singler_success) {
  cat("将细胞类型注释添加到Seurat对象...\n")

  # 将基于聚类的细胞类型注释添加到元数据中
  scObject <- AddMetaData(
    scObject,
    metadata = singler_result$cell_types,
    col.name = "cell_type"
  )

  # 同时添加为SingleR_labels和SingleR_celltype以保持兼容性
  scObject$SingleR_labels <- singler_result$cell_types
  scObject$SingleR_celltype <- singler_result$cell_types

  # 设置Idents为细胞类型用于下游分析
  Idents(scObject) <- singler_result$cell_types

  # 保存聚类-细胞类型对应表
  cluster_annotation_table <- singler_result$cluster_table
  colnames(cluster_annotation_table) <- c("Cluster", "CellType")

  # 打印细胞类型统计
  cell_type_table <- table(singler_result$cell_types)
  cat("\n细胞类型分布:\n")
  for(ct_name in names(cell_type_table)) {
    cat(sprintf("  %s: %d 细胞\n", ct_name, cell_type_table[ct_name]))
  }

  cat("细胞类型注释成功完成！\n")

} else {
  # 创建默认注释
  cat("创建默认的基于聚类的注释...\n")

  default_labels <- paste0("Cluster_", scObject$seurat_clusters)

  # 添加到元数据
  scObject <- AddMetaData(scObject, metadata = default_labels, col.name = "cell_type")
  scObject$SingleR_labels <- default_labels
  scObject$SingleR_celltype <- default_labels
  Idents(scObject) <- default_labels

  # 创建默认的聚类注释表
  cluster_annotation_table <- data.frame(
    Cluster = sort(unique(as.character(scObject$seurat_clusters))),
    CellType = paste0("Cluster_", sort(unique(as.character(scObject$seurat_clusters)))),
    stringsAsFactors = FALSE
  )

  cat("默认注释已创建\n")
}

# 验证注释是否成功添加
cat("\n验证注释结果:\n")
cat(sprintf("  SingleR_labels列存在: %s\n", "SingleR_labels" %in% colnames(scObject@meta.data)))
cat(sprintf("  cell_type列存在: %s\n", "cell_type" %in% colnames(scObject@meta.data)))
cat(sprintf("  SingleR_celltype列存在: %s\n", "SingleR_celltype" %in% colnames(scObject@meta.data)))
cat(sprintf("  唯一细胞类型数量: %d\n", length(unique(scObject$cell_type))))

cat(sprintf("自动注释完成，识别到 %d 种细胞类型\n",
            length(unique(Idents(scObject)))))

# 打印细胞类型统计
cell_type_counts <- table(Idents(scObject))
cat("\n细胞类型统计:\n")
for(ct in names(sort(cell_type_counts, decreasing = TRUE))) {
  cat(sprintf("  %s: %d 细胞\n", ct, cell_type_counts[ct]))
}

# =============================================
# 10.2 可视化注释结果
# =============================================

cat("\n生成注释结果可视化...\n")

# UMAP可视化
pdf(file = file.path(anno_dir, "01_Auto_Annotation_UMAP.pdf"), width = 12, height = 8)
print(
  DimPlot(
    scObject,
    reduction = "umap",
    pt.size = 1.5,
    label = TRUE,
    repel = TRUE
  ) +
  ggtitle("Cell Type Annotation") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black", face = "bold")
  )
)
dev.off()

# 保存注释结果
cat("保存注释结果...\n")

# 确保cluster_annotation_table存在
if(!exists("cluster_annotation_table")) {
  cluster_annotation_table <- data.frame(
    Cluster = sort(unique(as.character(scObject$seurat_clusters))),
    CellType = as.character(Idents(scObject))[match(sort(unique(as.character(scObject$seurat_clusters))),
                                                     as.character(scObject$seurat_clusters))],
    stringsAsFactors = FALSE
  )
  # 去重
  cluster_annotation_table <- cluster_annotation_table[!duplicated(cluster_annotation_table$Cluster), ]
}

write.csv(cluster_annotation_table,
          file = file.path(anno_dir, "CellCluster_AutoAnno_SingleR.csv"),
          row.names = FALSE)

# 生成细胞ID到细胞类型映射表
cat("生成细胞ID到细胞类型映射表...\n")

cell_id_to_celltype_mapping <- data.frame(
  Cell_ID = colnames(scObject),
  Cluster_ID = as.character(scObject$seurat_clusters),
  Cell_Type = as.character(Idents(scObject)),
  nFeature_RNA = scObject$nFeature_RNA,
  nCount_RNA = scObject$nCount_RNA,
  stringsAsFactors = FALSE
)

# 添加线粒体百分比
if ("percent.mito" %in% colnames(scObject@meta.data)) {
  cell_id_to_celltype_mapping$percent_mito <- scObject$percent.mito
} else if ("percent.mt" %in% colnames(scObject@meta.data)) {
  cell_id_to_celltype_mapping$percent_mito <- scObject$percent.mt
}

# 按聚类ID排序
cell_id_to_celltype_mapping <- cell_id_to_celltype_mapping[order(as.numeric(cell_id_to_celltype_mapping$Cluster_ID)), ]

write.csv(cell_id_to_celltype_mapping,
          file = file.path(output_dir, "08_Data_Export", "Cell_ID_to_CellType_Mapping.csv"),
          row.names = FALSE)

# 输出映射表信息
cat("细胞ID到细胞类型映射表已生成\n")
cat("- 总细胞数:", nrow(cell_id_to_celltype_mapping), "\n")
cat("- 细胞类型数:", length(unique(cell_id_to_celltype_mapping$Cell_Type)), "\n")
cat("- 细胞类型:", paste(unique(cell_id_to_celltype_mapping$Cell_Type), collapse = ", "), "\n")

# 显示每个细胞类型的细胞数量
cat("\n=== 各细胞类型细胞数量统计 ===\n")
celltype_counts <- table(cell_id_to_celltype_mapping$Cell_Type)
for(celltype in names(celltype_counts)) {
  cat(celltype, ":", celltype_counts[celltype], "个细胞\n")
}

# =============================================
# 生成Nature风格的细胞类型UMAP图
# =============================================
cat("\n生成Nature风格的细胞类型UMAP图...\n")
tryCatch({
  generate_celltype_umap(scObject, umap_viz_dir, celltype_col = "cell_type")
  cat("Nature风格UMAP图生成完成！\n")
}, error = function(e) {
  cat("生成Nature风格UMAP图时出错:", e$message, "\n")
})

# =============================================
# 生成简洁风格UMAP图（类似参考图片风格）
# 特点：柔和配色、标签直接标注、底部显示细胞数、显示坐标轴
# =============================================
cat("\n生成简洁风格UMAP图（Publication Style）...\n")

tryCatch({
  # 提取UMAP坐标
  umap_data <- as.data.frame(Embeddings(scObject, "umap"))
  colnames(umap_data) <- c("UMAP_1", "UMAP_2")
  umap_data$cell_type <- as.character(scObject$cell_type)
  umap_data$cluster <- as.character(scObject$seurat_clusters)

  # 获取细胞总数
  total_cells <- ncol(scObject)
  total_cells_label <- format(total_cells, big.mark = ",")

  # 定义柔和的颜色方案（类似参考图片）
  soft_palette <- c(
    "#7EB6D9",  # 浅蓝色 (T cell)
    "#F4A460",  # 橙色 (B cell)
    "#98D98E",  # 浅绿色 (Cancer cell)
    "#DDA0DD",  # 浅紫色 (Myeloid cell)
    "#F0E68C",  # 浅黄色 (Plasma cell)
    "#E9967A",  # 深橙色 (Mast cell)
    "#87CEEB",  # 天蓝色 (NK)
    "#D3D3D3",  # 浅灰色 (Neutrophil)
    "#FFB6C1",  # 浅粉色 (Endothelial cell)
    "#B8860B",  # 深金色 (Epithelial cell)
    "#9ACD32",  # 黄绿色 (Stromal cell)
    "#BC8F8F",  # 玫瑰棕 (Schwann cell)
    "#ADD8E6",  # 淡蓝色
    "#FFDAB9",  # 桃色
    "#E6E6FA",  # 淡紫色
    "#F5DEB3",  # 小麦色
    "#D8BFD8",  # 蓟色
    "#AFEEEE",  # 淡青色
    "#FFC0CB",  # 粉红色
    "#F0FFF0"   # 蜜瓜色
  )

  # ----------- 图1: 细胞类型UMAP图 -----------
  cat("  生成细胞类型UMAP图...\n")

  # 获取唯一细胞类型并分配颜色
  unique_celltypes <- unique(umap_data$cell_type)
  n_celltypes <- length(unique_celltypes)
  celltype_colors <- setNames(soft_palette[1:n_celltypes], unique_celltypes)

  # 计算每个细胞类型的中心位置用于标签
  celltype_centers <- umap_data %>%
    group_by(cell_type) %>%
    summarise(
      UMAP_1 = median(UMAP_1),
      UMAP_2 = median(UMAP_2),
      n_cells = n(),
      .groups = "drop"
    )

  # 绘制细胞类型UMAP图（带坐标轴）
  p_celltype_simple <- ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = cell_type)) +
    geom_point(size = 0.3, alpha = 0.6) +
    scale_color_manual(values = celltype_colors) +
    # 添加细胞类型标签（黑色文字）
    geom_text(data = celltype_centers,
              aes(x = UMAP_1, y = UMAP_2, label = cell_type),
              color = "black",
              size = 3.5,
              fontface = "bold",
              inherit.aes = FALSE) +
    # 坐标轴标签和底部细胞数
    labs(x = "UMAP_1", y = "UMAP_2", caption = paste0(total_cells_label, " cells")) +
    # 简洁主题：显示坐标轴
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      axis.line = element_line(color = "black", size = 0.5),
      plot.caption = element_text(hjust = 0.5, size = 12, color = "black",
                                  margin = margin(t = 15)),
      plot.margin = margin(20, 20, 20, 20)
    )

  # 保存细胞类型UMAP图（仅PDF）
  ggsave(file.path(anno_dir, "UMAP_CellType_SimpleStyle.pdf"),
         plot = p_celltype_simple, width = 8, height = 8)

  cat("  细胞类型UMAP图已保存: UMAP_CellType_SimpleStyle.pdf\n")

  # ----------- 图2: 聚类UMAP图 -----------
  cat("  生成聚类UMAP图...\n")

  # 获取唯一聚类并分配颜色
  unique_clusters <- sort(unique(as.numeric(umap_data$cluster)))
  n_clusters <- length(unique_clusters)
  cluster_colors <- setNames(soft_palette[1:n_clusters], as.character(unique_clusters))

  # 计算每个聚类的中心位置用于标签
  cluster_centers <- umap_data %>%
    group_by(cluster) %>%
    summarise(
      UMAP_1 = median(UMAP_1),
      UMAP_2 = median(UMAP_2),
      n_cells = n(),
      .groups = "drop"
    )

  # 绘制聚类UMAP图（带坐标轴）
  p_cluster_simple <- ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
    geom_point(size = 0.3, alpha = 0.6) +
    scale_color_manual(values = cluster_colors) +
    # 添加聚类标签（黑色文字）
    geom_text(data = cluster_centers,
              aes(x = UMAP_1, y = UMAP_2, label = cluster),
              color = "black",
              size = 4,
              fontface = "bold",
              inherit.aes = FALSE) +
    # 坐标轴标签和底部细胞数
    labs(x = "UMAP_1", y = "UMAP_2", caption = paste0(total_cells_label, " cells")) +
    # 简洁主题：显示坐标轴
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      axis.line = element_line(color = "black", size = 0.5),
      plot.caption = element_text(hjust = 0.5, size = 12, color = "black",
                                  margin = margin(t = 15)),
      plot.margin = margin(20, 20, 20, 20)
    )

  # 保存聚类UMAP图（仅PDF）
  ggsave(file.path(anno_dir, "UMAP_Cluster_SimpleStyle.pdf"),
         plot = p_cluster_simple, width = 8, height = 8)

  cat("  聚类UMAP图已保存: UMAP_Cluster_SimpleStyle.pdf\n")

  # ----------- 图3: 聚类与细胞类型并排对比图 -----------
  cat("  生成聚类与细胞类型对比图...\n")

  # 合并两图
  combined_simple <- p_cluster_simple + p_celltype_simple +
    patchwork::plot_layout(ncol = 2)

  ggsave(file.path(anno_dir, "UMAP_Cluster_vs_CellType_SimpleStyle.pdf"),
         plot = combined_simple, width = 16, height = 8)

  cat("  对比图已保存: UMAP_Cluster_vs_CellType_SimpleStyle.pdf\n")

  cat("简洁风格UMAP图生成完成！\n")

}, error = function(e) {
  cat("生成简洁风格UMAP图时出错:", e$message, "\n")
})

# 合并marker和注释数据
if(!is.null(marker_df)) {
  marker_df$cluster <- as.character(marker_df$cluster)
  marker_df_anno <- marker_df %>%
    left_join(cluster_annotation_table, by = c("cluster" = "Cluster"))
  write.csv(marker_df_anno, file=file.path(output_dir, "04_Differential_Expression","Cluster_Markers_DEGs_with_CellType.csv"), row.names=FALSE)
}

# ===============================================================================
# PCA散点图（PCA Scatter Plot）- 展示细胞在PCA降维空间中的分布
# ===============================================================================

cat("Step 10.5: 绘制PCA散点图...\n")

# 检查cell_type列是否存在，如果不存在则报错中断
if (!"cell_type" %in% colnames(scObject@meta.data)) {
  stop("错误: cell_type列不存在！请检查细胞类型注释步骤是否成功完成。")
}

# 获取细胞类型颜色
unique_celltypes <- unique(as.character(scObject$cell_type))
celltype_colors <- colorRampPalette(brewer.pal(12, "Set3"))(length(unique_celltypes))
names(celltype_colors) <- unique_celltypes

# PCA散点图 - 按细胞类型着色
pdf(file.path(output_dir, "03_Cell_Type_Annotation", "PCA_scatter_by_celltype.pdf"), width = 10, height = 8)
p_pca_celltype <- DimPlot(
  object = scObject,
  reduction = "pca",
  group.by = "cell_type",
  cols = celltype_colors,
  pt.size = 1.5,
  label = FALSE
) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = ggplot2::margin(b = 15)),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    axis.line = element_line(color = "black", size = 0.8),
    axis.ticks = element_line(color = "black", size = 0.5),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  ) +
  labs(
    title = "PCA Scatter Plot - Colored by Cell Type",
    x = "PC1",
    y = "PC2",
    color = "Cell Type"
  )
print(p_pca_celltype)
dev.off()
cat("已保存: PCA_scatter_by_celltype.pdf\n")

# PCA散点图 - 按聚类着色
cluster_levels <- sort(unique(as.character(scObject$seurat_clusters)))
cluster_colors <- colorRampPalette(brewer.pal(12, "Paired"))(length(cluster_levels))
names(cluster_colors) <- cluster_levels

pdf(file.path(output_dir, "02_Clustering_Analysis", "PCA_scatter_by_cluster.pdf"), width = 10, height = 8)
p_pca_cluster <- DimPlot(
  object = scObject,
  reduction = "pca",
  group.by = "seurat_clusters",
  cols = cluster_colors,
  pt.size = 1.5,
  label = TRUE,
  label.size = 4
) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = ggplot2::margin(b = 15)),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    axis.line = element_line(color = "black", size = 0.8),
    axis.ticks = element_line(color = "black", size = 0.5),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  ) +
  labs(
    title = "PCA Scatter Plot - Colored by Cluster",
    x = "PC1",
    y = "PC2",
    color = "Cluster"
  )
print(p_pca_cluster)
dev.off()
cat("已保存: PCA_scatter_by_cluster.pdf\n")

# PCA散点图 - 合并版本（细胞类型和聚类左右排列）
pdf(file.path(output_dir, "03_Cell_Type_Annotation", "PCA_scatter_combined.pdf"), width = 18, height = 8)
combined_pca <- p_pca_celltype + p_pca_cluster +
  patchwork::plot_layout(ncol = 2) +
  patchwork::plot_annotation(
    title = "PCA Sample Distribution",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
  )
print(combined_pca)
dev.off()
cat("已保存: PCA_scatter_combined.pdf\n")

cat("PCA散点图绘制完成！\n\n")

# ===============================================================================
# 美化UMAP细胞类型图 - 带圆形边界（使用scRNAtoolVis包）
# ===============================================================================

cat("绘制美化版UMAP细胞类型图...\n")

# 检查细胞类型数据
cat("检查细胞类型数据...\n")
cat("细胞类型列表：", unique(scObject$SingleR_celltype), "\n")
cat("细胞类型数量：", length(unique(scObject$SingleR_celltype)), "\n")

# 处理NA值
if(any(is.na(scObject$SingleR_celltype))) {
  cat("发现NA值，将其替换为'Unknown'...\n")
  scObject$SingleR_celltype[is.na(scObject$SingleR_celltype)] <- "Unknown"
  cat("处理后的细胞类型列表：", unique(scObject$SingleR_celltype), "\n")
}

# 使用scRNAtoolVis绘制带圆形边界的UMAP图
cat("开始绘制UMAP图...\n")

# 第一个图：聚类边界
cat("绘制聚类边界版本...\n")
pdf(file = file.path(output_dir, "03_Cell_Type_Annotation", "cluster_with_circle_boundary.pdf"), width = 12, height = 9)

clusterCornerAxes(object = scObject,
                  reduction = 'umap',
                  clusterCol = "seurat_clusters",
                  noSplit = TRUE,
                  cellLabel = TRUE,
                  cellLabelSize = 4,
                  cornerTextSize = 3,
                  themebg = 'bwCorner',
                  show.legend = TRUE,
                  keySize = 6,
                  aspect.ratio = 1,
                  relLength = 0.3,
                  addCircle = FALSE) +
  ggtitle("Seurat Clusters")

dev.off()
cat("聚类边界图绘制完成！\n")

# 第二个图：细胞类型边界（按聚类画圈，显示细胞类型）
cat("绘制细胞类型边界版本...\n")
pdf(file = file.path(output_dir, "03_Cell_Type_Annotation", "celltype_with_circle_boundary.pdf"), width = 12, height = 9)

# 从CSV文件读取正确的注释
anno_path <- file.path(output_dir, "03_Cell_Type_Annotation", "CellCluster_AutoAnno_SingleR.csv")
if(file.exists(anno_path)) {
  # 读取正确的注释
  anno_df <- read.csv(anno_path, stringsAsFactors = FALSE)
  cat("从CSV读取的注释：\n")
  print(anno_df)
  
  # 创建聚类到细胞类型的正确映射
  cluster_to_celltype <- setNames(anno_df$CellType, as.character(anno_df$Cluster))
  
  # 为每个细胞分配正确的细胞类型
  correct_celltype <- cluster_to_celltype[as.character(scObject$seurat_clusters)]
  
  # 最简单直接的方法：先按聚类画圆圈，然后手动修改颜色和图例
  p_base <- clusterCornerAxes(object = scObject,
                             reduction = 'umap',
                             clusterCol = "seurat_clusters",      # 按聚类显示
                             noSplit = TRUE,
                             cellLabel = FALSE,                   # 暂不显示标签
                             cornerTextSize = 3,
                             themebg = 'bwCorner',
                             show.legend = FALSE,                 # 暂不显示图例
                             aspect.ratio = 1,
                             relLength = 0.3,
                             addCircle = FALSE)                   # 不添加圆圈边界
  
  # 获取数据并添加正确的细胞类型信息
  umap_data <- data.frame(
    UMAP_1 = Embeddings(scObject, "umap")[,1],
    UMAP_2 = Embeddings(scObject, "umap")[,2],
    CellType = correct_celltype  # 使用从CSV读取的正确注释
  )
  
  # 计算细胞类型中心位置
  celltype_centers <- aggregate(cbind(UMAP_1, UMAP_2) ~ CellType, data = umap_data, FUN = mean)
  
  # 设置颜色
  unique_celltypes <- unique(correct_celltype)
  celltype_colors <- colorRampPalette(brewer.pal(12, "Set3"))(length(unique_celltypes))
  names(celltype_colors) <- unique_celltypes
  
  # 在基础图上添加细胞类型着色点和标签
  p_final <- p_base + 
    geom_point(data = umap_data, aes(x = UMAP_1, y = UMAP_2, color = CellType), 
               size = 0.6, alpha = 0.8) +
    scale_color_manual(values = celltype_colors, name = "Cell Type") +
    # 使用ggrepel避免标签重叠
    ggrepel::geom_text_repel(data = celltype_centers, 
                            aes(x = UMAP_1, y = UMAP_2, label = CellType), 
                            color = "black", 
                            size = 3.5, 
                            fontface = "bold",
                            box.padding = 0.5,        # 标签与点的距离
                            point.padding = 0.3,      # 标签之间的距离
                            segment.color = "grey50",  # 连接线颜色
                            segment.size = 0.5,       # 连接线粗细
                            max.overlaps = 20) +      # 允许的最大重叠数
    theme(legend.position = "right", legend.key.size = unit(0.8, "cm")) +
    ggtitle("Cell Types with All Cluster Boundaries")
  
  print(p_final)
  
  cat("使用CSV文件中的正确注释绘制完成！\n")
  
} else {
  cat("警告：CSV文件不存在，使用Seurat对象中的注释\n")
  
  # 备用方法（如果CSV不存在）
  p_base <- clusterCornerAxes(object = scObject,
                             reduction = 'umap',
                             clusterCol = "seurat_clusters",
                             noSplit = TRUE,
                             cellLabel = TRUE,
                             cellLabelSize = 4,
                             cornerTextSize = 3,
                             themebg = 'bwCorner',
                             show.legend = TRUE,
                             keySize = 6,
                             aspect.ratio = 1,
                             relLength = 0.3,
                             addCircle = TRUE,
                             cicAlpha = 0.4,
                             cicDelta = 2.0,
                             nbin = 500) +
    ggtitle("Cell Types with Cluster Boundaries (Fallback)")
  
  print(p_base)
}

dev.off()
cat("细胞类型边界图绘制完成！\n")

# ===============================================================================
# 第七部分：基因聚类富集分析
# ===============================================================================

cat("Step 11: 基因聚类富集分析...\n")

# 检查是否有显著的marker基因进行聚类分析
if(!is.null(marker_df) && nrow(marker_df) > 0 && exists("sigMarkers") && nrow(sigMarkers) > 0) {
  
  cat("正在进行基因聚类分析...\n")
  
  # 获取top marker基因用于聚类分析
  top_markers <- marker_df %>%
    group_by(cluster) %>%
    arrange(desc(abs(avg_log2FC)), p_val_adj) %>%
    slice_head(n = 20)  # 每个cluster取top 20个基因
  
  # 获取标准化的表达数据
  normalized_data <- GetAssayData(scObject, layer = "data")
  
  # 筛选marker基因的表达矩阵
  marker_genes <- unique(top_markers$gene)
  marker_genes <- intersect(marker_genes, rownames(normalized_data))
  cat("用于聚类分析的标记基因数量:", length(marker_genes), "\n")
  
  if(length(marker_genes) > 10) {  # 至少需要10个基因进行聚类
    
    # 提取标记基因的表达矩阵
    gene_expression_matrix <- as.matrix(normalized_data[marker_genes, ])
    
    # 强制使用细胞类型注释而不是聚类信息
    # 从CSV注释文件获取正确的细胞类型信息
    anno_path <- file.path(output_dir, "03_Cell_Type_Annotation", "CellCluster_AutoAnno_SingleR.csv")
    use_celltype <- FALSE
    
    if(file.exists(anno_path)) {
      # 读取注释文件
      anno_df <- read.csv(anno_path, stringsAsFactors = FALSE)
      
      # 创建聚类到细胞类型的映射
      cluster_to_celltype <- setNames(anno_df$CellType, anno_df$Cluster)
      
      # 为每个细胞分配正确的细胞类型
      cell_clusters <- as.character(scObject@meta.data$seurat_clusters)
      cell_types_corrected <- cluster_to_celltype[cell_clusters]
      names(cell_types_corrected) <- colnames(scObject)
      
      # 检查是否成功获取细胞类型
      if(!all(is.na(cell_types_corrected))) {
        use_celltype <- TRUE
        cell_types <- cell_types_corrected
        unique_cell_types <- unique(cell_types[!is.na(cell_types)])
        
        cat("使用细胞类型注释进行基因聚类分析，包含", length(unique_cell_types), "种细胞类型\n")
        cat("细胞类型：", paste(unique_cell_types, collapse = ", "), "\n")
        
        # 构建平均表达矩阵（按细胞类型）
        avg_expression <- matrix(0, nrow = length(marker_genes), ncol = length(unique_cell_types))
        rownames(avg_expression) <- marker_genes
        colnames(avg_expression) <- unique_cell_types
        
        for (ct in unique_cell_types) {
          cells_of_type <- which(cell_types == ct & !is.na(cell_types))
          if (length(cells_of_type) > 1) {
            avg_expression[, ct] <- rowMeans(gene_expression_matrix[, cells_of_type])
          } else if (length(cells_of_type) == 1) {
            avg_expression[, ct] <- gene_expression_matrix[, cells_of_type]
          }
        }
      }
    }
    
    # 如果没有成功获取细胞类型注释，回退到聚类信息
    if(!use_celltype) {
      cat("警告：无法获取细胞类型注释，使用聚类信息进行基因聚类分析\n")
      clusters <- scObject@meta.data$seurat_clusters
      unique_clusters <- sort(unique(clusters))
      
      # 构建平均表达矩阵（按聚类）
      avg_expression <- matrix(0, nrow = length(marker_genes), ncol = length(unique_clusters))
      rownames(avg_expression) <- marker_genes
      colnames(avg_expression) <- paste0("Cluster_", unique_clusters)
      
      for (cl in unique_clusters) {
        cells_of_cluster <- which(clusters == cl)
        if (length(cells_of_cluster) > 1) {
          avg_expression[, paste0("Cluster_", cl)] <- rowMeans(gene_expression_matrix[, cells_of_cluster])
        } else if (length(cells_of_cluster) == 1) {
          avg_expression[, paste0("Cluster_", cl)] <- gene_expression_matrix[, cells_of_cluster]
        }
      }
    }
    
    cat("已构建平均表达矩阵，维度:", dim(avg_expression), "\n")
    
    # 使用ClusterGVis进行基因聚类分析
    cat("正在进行基因聚类分析...\n")
    
    # 确保数据格式正确
    avg_expression <- as.data.frame(avg_expression)
    
    # 聚类数量（自动适配基因数，避免报错）
    max_cluster <- min(12, nrow(avg_expression) %/% 5)
    cluster_number <- max(4, max_cluster)
    
    cat("基因数量:", nrow(avg_expression), "，细胞类型/簇数量:", ncol(avg_expression), "\n")
    cat("设置聚类数量为:", cluster_number, "\n")
    
    # 进行基因聚类（稳定版）
    tryCatch({
      gene_clusters <- clusterData(obj = avg_expression,
                                   cluster.method = "kmeans",
                                   cluster.num = cluster_number)
    }, error = function(e) {
      cat("kmeans 失败，改用 hclust ...\n")
      gene_clusters <<- clusterData(obj = avg_expression,
                                    cluster.method = "hclust",
                                    cluster.num = cluster_number)
    })
    
    cat("基因聚类完成，共分为", cluster_number, "个模块\n")
    
    # GO富集分析（稳定版，不会崩）
    cat("正在进行GO富集分析...\n")
    gene_enrich <- NULL
    
    tryCatch({
      library(clusterProfiler)
      library(org.Hs.eg.db)
      
      gene_enrich <- enrichCluster(
        object = gene_clusters,
        OrgDb = org.Hs.eg.db,
        type = "BP",
        organism = "hsa",
        pvalueCutoff = 0.05,
        topn = 5,
        seed = 5201314
      )
      
      if(!is.null(gene_enrich) && nrow(gene_enrich) > 0) {
        cat("GO富集完成，发现", nrow(gene_enrich), "个显著条目\n")
      } else {
        cat("GO富集完成，但无显著条目\n")
        gene_enrich <- NULL
      }
      
    }, error = function(e) {
      cat("GO富集失败，原因：", e$message, "\n")
      cat("仅绘制基础热图\n")
      gene_enrich <- NULL
    })
    
    # 标记基因（稳定版）
    mark_genes <- c()
    if(exists("important_marker_genes")) {
      mark_genes <- intersect(important_marker_genes, rownames(avg_expression))
    }
    
    if(length(mark_genes) < 10) {
      n_add <- min(15 - length(mark_genes), nrow(avg_expression))
      if(n_add > 0) {
        add_genes <- sample(setdiff(rownames(avg_expression), mark_genes), n_add)
        mark_genes <- c(mark_genes, add_genes)
      }
    }
    cat("将标记", length(mark_genes), "个重要基因\n")
    
    # 输出文件夹
    out_dir <- file.path(output_dir, "05_Gene_Cluster_Enrichment")
    if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    
    # 绘图（稳定不报错）
    cat("正在生成热图...\n")
    pdf(file.path(out_dir, "gene_cluster_enrichment_heatmap.pdf"),
        width = 16, height = 12, onefile = F)
    
    if(!is.null(gene_enrich)) {
      visCluster(
        object = gene_clusters,
        plot.type = "both",
        column_names_rot = 45,
        show_row_dend = FALSE,
        markGenes = mark_genes,
        markGenes.side = "left",
        genes.gp = c('italic', fontsize = 10, col = "black"),
        annoTerm.data = gene_enrich,
        line.side = "left"
      )
    } else {
      visCluster(
        object = gene_clusters,
        plot.type = "both",
        column_names_rot = 45,
        show_row_dend = FALSE,
        markGenes = mark_genes,
        markGenes.side = "left",
        genes.gp = c('italic', fontsize = 10, col = "black"),
        line.side = "left"
      )
    }
    
    dev.off()
    cat("热图已保存\n")
    
    # 保存结果
    write.table(gene_clusters$wide.res,
                file.path(out_dir, "gene_cluster_results.txt"),
                sep = "\t", quote = F, row.names = T)
    
    if(!is.null(gene_enrich)) {
      write.table(gene_enrich,
                  file.path(out_dir, "gene_enrichment_results.txt"),
                  sep = "\t", quote = F, row.names = F)
    }
    
    write.table(avg_expression,
                file.path(out_dir, "average_expression_matrix.txt"),
                sep = "\t", quote = F, row.names = T)
    
    cat("基因聚类富集分析完成！\n")
    
  } else {
    cat("警告：标记基因数量不足（<10），跳过\n")
  }
  
} else {
  cat("警告：无足够marker基因，跳过\n")
}

cat("Step 12: 各种统计可视化输出...\n")
meta_df <- scObject@meta.data

# 从CSV文件读取正确的细胞类型注释
anno_path <- file.path(output_dir, "03_Cell_Type_Annotation", "CellCluster_AutoAnno_SingleR.csv")
if(file.exists(anno_path)) {
  anno_df <- read.csv(anno_path, stringsAsFactors = FALSE)
  cluster_to_celltype <- setNames(anno_df$CellType, anno_df$Cluster)
  meta_df$CellType_corrected <- cluster_to_celltype[as.character(meta_df$seurat_clusters)]
  
  # 使用正确的细胞类型注释进行统计
  celltype_stats <- as.data.frame(table(
    Sample = meta_df$orig.ident, 
    CellType = meta_df$CellType_corrected
  ))
} else {
  # 如果CSV文件不存在，使用原始注释
  celltype_stats <- as.data.frame(table(
    Sample = meta_df$orig.ident, 
    CellType = meta_df$SingleR_celltype
  ))
}

celltype_stats <- celltype_stats %>%
  group_by(Sample) %>%
  mutate(Ratio = Freq / sum(Freq)) %>%
  ungroup()

celltype_stats$Sample <- factor(celltype_stats$Sample, levels = unique(celltype_stats$Sample))
celltype_stats$CellType <- factor(celltype_stats$CellType, levels = unique(celltype_stats$CellType))
cell_types <- levels(celltype_stats$CellType)
nColors <- length(cell_types)
myColors <- colorRampPalette(brewer.pal(12, "Set3"))(nColors)
names(myColors) <- cell_types
pdf(file.path(output_dir, "09_Statistics_Plots", "Sample_CellType_Composition_barplot.pdf"), width=7, height=6)
ggplot(celltype_stats, aes(x = Sample, y = Ratio, fill = CellType)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7, color = "grey70") +
  scale_y_continuous(expand = c(0,0), labels = scales::percent, name = "Proportion") +
  scale_fill_manual(values = myColors, name="Cell Type") +
  labs(x = "Sample") +
  theme_classic(base_size=16) +
  coord_flip() +  
  theme(
    legend.title=element_text(size=14),
    legend.text=element_text(size=13),
    axis.text = element_text(size=13),
    axis.title = element_text(size=14)
  )
dev.off()

if(exists("anno_df") && !is.null(anno_df)) {
  cell_counts <- as.data.frame(table(Sample = meta_df$orig.ident, CellType = meta_df$CellType_corrected))
}
write.csv(cell_counts, file.path(output_dir, "09_Statistics_Plots", "CellNumber_perSample_CellType.csv"), row.names = FALSE)

celltype_matrix <- celltype_stats %>% dplyr::select(Sample, CellType, Ratio) %>%
  tidyr::spread(CellType, Ratio, fill = 0)
write.csv(celltype_matrix, file.path(output_dir, "09_Statistics_Plots", "CellType_Proportion_Matrix.csv"), row.names=FALSE)
if(exists("anno_df") && !is.null(anno_df)) {
  cell_annot <- data.frame(CellBarcode = rownames(meta_df), Sample = meta_df$orig.ident,
                           Cluster = meta_df$seurat_clusters, CellType = meta_df$CellType_corrected)
}
write.csv(cell_annot, file.path(output_dir, "08_Data_Export", "Cell_Full_Annotation.csv"), row.names=FALSE)

# marker计数及绘图
if (exists("marker_df_anno") && nrow(marker_df_anno) > 0) {
  gene_count_by_celltype <- marker_df_anno %>%
    group_by(CellType) %>%
    summarise(gene_number = n_distinct(gene)) %>%
    arrange(desc(gene_number))
  write.csv(gene_count_by_celltype, file = file.path(output_dir, "09_Statistics_Plots", "GeneNumber_per_CellType.csv"), row.names = FALSE)
}

# PCA方差解释
pca_stdev <- scObject[["pca"]]@stdev
pca_var_explained <- (pca_stdev^2) / sum(pca_stdev^2)
cumulative_var <- cumsum(pca_var_explained)
pca_table <- data.frame(
  PC = seq_along(pca_var_explained),
  Variance_Explained = pca_var_explained,
  Cumulative_Variance = cumulative_var
)
write.csv(pca_table, file = file.path(output_dir, "08_Data_Export", "PCA_Variance_Explained.csv"), row.names=FALSE)

opt_cutoff <- 0.8
best_n_pcs <- which(cumulative_var >= opt_cutoff)[1]
write.csv(
  data.frame(Recommended_n_PCs = best_n_pcs, Cumulative_Variance = cumulative_var[best_n_pcs]),
  file = file.path(output_dir, "08_Data_Export", "PCA_Recommended_nPCs.csv"),
  row.names=FALSE
)

# ===============================================================================
# 第九部分：目标基因分析
# ===============================================================================

cat("Step 13: 目标基因分析...\n")

for (target_gene in target_genes) {
  if(! target_gene %in% rownames(scObject)) {
    cat("警告：指定基因", target_gene, "不在表达矩阵中！\n")
    next
  }
  
  # 创建基因专属目录
  gene_dir <- file.path(output_dir, "07_Target_Gene_Analysis", target_gene)
  if (!dir.exists(gene_dir)) dir.create(gene_dir, recursive = TRUE)
  
  # UMAP表达分布
  pdf(file.path(gene_dir, paste0(target_gene, "_UMAP_FeaturePlot.pdf")), width=7, height=6)
  print(
    FeaturePlot(
      scObject, 
      features = target_gene, 
      reduction = "umap", 
      cols = c("lightgrey", "red"), 
      pt.size = 1,
      label = TRUE
    ) + 
      ggtitle(paste(target_gene, "Expression (UMAP)"))
  )
  dev.off()
  
  # 按聚类的小提琴图
  pdf(file.path(gene_dir, paste0(target_gene, "_VlnPlot_byCluster.pdf")), width=8, height=6)
  p2 <- VlnPlot(
    scObject, 
    features = target_gene, 
    group.by = "seurat_clusters", 
    pt.size = 0,
    slot = "data"
  ) +
    geom_boxplot(width=0.15, outlier.size=0, fill=NA, color="black", lwd=0.7) +
    theme_classic(base_size=16) +
    theme(
      legend.position = "right",
      axis.title.x = element_text(size = 15, face = "bold"),
      axis.title.y = element_text(size = 15, face = "bold"),
      axis.text.x  = element_text(size = 12, face = "bold", angle = 45, vjust=1, hjust=1),
      axis.text.y  = element_text(size = 12, face = "bold"),
      strip.text   = element_text(size = 15, face = "bold", color = "black")
    ) +
    labs(x = NULL, y = "Normalized Expression") +
    ggtitle(paste(target_gene, "by Cluster"))
  print(p2)
  dev.off()
  
  # 导出分组均值
  avg_by_celltype <- AverageExpression(scObject, features = target_gene, group.by = "SingleR_celltype")$RNA
  avg_by_cluster  <- AverageExpression(scObject, features = target_gene, group.by = "seurat_clusters")$RNA
  
  write.csv(avg_by_celltype, file.path(gene_dir, paste0(target_gene, "_Average_By_CellType.csv")))
  write.csv(avg_by_cluster,  file.path(gene_dir, paste0(target_gene, "_Average_By_Cluster.csv")))
  
  # 导出每个细胞表达量
  gene_expr_vec <- FetchData(scObject, vars = target_gene)
  cell_anno <- data.frame(
    CellBarcode = colnames(scObject),
    CellType = scObject$SingleR_celltype,
    Cluster = scObject$seurat_clusters,
    Gene_Expression = as.vector(gene_expr_vec[,1])
  )
  write.csv(cell_anno, file.path(gene_dir, paste0(target_gene, "_Expression_perCell.csv")), row.names=FALSE)
  
  # ====== 自定义UMAP+条形图（聚类颜色一致，细胞类型全名显示） ======
  # 读取注释表
  anno_path <- file.path(output_dir, "03_Cell_Type_Annotation", "CellCluster_AutoAnno_SingleR.csv")
  if(file.exists(anno_path)) {
    anno_df <- read.csv(anno_path, stringsAsFactors = FALSE)
    cluster_col <- grep("cluster", colnames(anno_df), ignore.case = TRUE, value = TRUE)[1]
    celltype_col <- grep("cell.*type", colnames(anno_df), ignore.case = TRUE, value = TRUE)[1]
    
    if(!is.na(cluster_col) && !is.na(celltype_col)) {
      anno_df[[cluster_col]] <- as.character(anno_df[[cluster_col]])
      
      # UMAP数据
      umap_df <- as.data.frame(Embeddings(scObject, "umap"))
      umap_cols <- grep("UMAP|Dim", colnames(umap_df), value = TRUE, ignore.case = TRUE)
      if(length(umap_cols) >= 2) {
        colnames(umap_df)[match(umap_cols[1:2], colnames(umap_df))] <- c("UMAP_1", "UMAP_2")
      } else if(ncol(umap_df) == 2) {
        colnames(umap_df) <- c("UMAP_1", "UMAP_2")
      }
      
      umap_df$Cluster <- as.character(scObject$seurat_clusters)
      umap_df$CellBarcode <- rownames(umap_df)
      umap_df$Expression <- FetchData(scObject, vars = target_gene)[,1]
      
      # 合并注释
      umap_df <- left_join(umap_df, anno_df, by = setNames(cluster_col, "Cluster"))
      
      if(celltype_col %in% colnames(umap_df)) {
        # 统计每个聚类的平均表达
        cluster_stats <- umap_df %>%
          group_by(Cluster, .data[[celltype_col]]) %>%
          summarise(AvgExpr = mean(Expression), .groups = "drop") %>%
          arrange(Cluster)
        colnames(cluster_stats)[2] <- "CellType"
        
        # 计算聚类中心用于标注
        label_df <- umap_df %>%
          group_by(Cluster, .data[[celltype_col]]) %>%
          summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2), .groups = "drop")
        colnames(label_df)[2] <- "CellType"
        
        # 颜色配色
        cluster_levels <- sort(unique(umap_df$Cluster))
        cluster_colors <- setNames(colorRampPalette(brewer.pal(12, "Set3"))(length(cluster_levels)), cluster_levels)
        
        # UMAP主图
        p_umap <- ggplot(umap_df, aes(x=UMAP_1, y=UMAP_2)) +
          geom_point(aes(color=factor(Cluster)), size=1, alpha=0.7) +
          scale_color_manual(values=cluster_colors, name="Cluster") +
          geom_text(data=label_df,
                    aes(label=paste0("c-", Cluster, "\n", CellType)),
                    color="black", size=4, fontface="bold") +
          theme_classic(base_size=16) +
          ggtitle(target_gene) +
          theme(
            plot.title = element_text(size=22, face="bold", hjust=0.5),
            legend.position = "none"
          )
        
        # 条形图
        cluster_stats$ClusterLabel <- paste0("c-", cluster_stats$Cluster)
        cluster_stats$ClusterLabel <- factor(cluster_stats$ClusterLabel, levels = paste0("c-", cluster_levels))
        cluster_stats$CellTypeLabel <- cluster_stats$CellType
        max_expr <- max(cluster_stats$AvgExpr)
        offset <- max_expr * 0.03  # 右移3%用于显示全名
        
        cluster_stats$Cluster <- as.character(cluster_stats$Cluster)
        cluster_stats$ClusterLabel <- paste0("c-", cluster_stats$Cluster)
        cluster_stats$ClusterLabel <- factor(cluster_stats$ClusterLabel, levels = paste0("c-", cluster_levels))
        cluster_stats$CellTypeLabel <- cluster_stats$CellType
        
        p_legend <- ggplot(cluster_stats, aes(x=ClusterLabel, y=AvgExpr, fill=Cluster)) +
          geom_bar(stat="identity", width=0.7) +
          scale_fill_manual(values=cluster_colors, guide="none") +
          geom_text(aes(x=ClusterLabel, y=AvgExpr + offset, label=CellTypeLabel),
                    hjust=0, vjust=0.5, size=5, fontface="plain", color="black") +
          coord_flip(clip = "off") +
          labs(x=NULL, y="Avg Expression", title="") +
          theme_minimal(base_size=14) +
          theme(
            axis.text.y = element_text(size=13, face="bold"),
            axis.text.x = element_text(size=12),
            axis.title.x = element_text(size=14),
            plot.margin = margin(5, 80, 5, 5)
          )
        
        # 拼图输出
        p_final <- p_umap + p_legend + patchwork::plot_layout(widths = c(2.5, 1.2))
        
        pdf(file.path(gene_dir, paste0(target_gene, "_UMAP_withLegend.pdf")), width=12, height=6)
        print(p_final)
        dev.off()
        cat("已输出自定义UMAP表达图：", file.path(gene_dir, paste0(target_gene, "_UMAP_withLegend.pdf")), "\n")
      }
    }
  }
}

# ===============================================================================
# 第十部分：轨迹分析
# ===============================================================================

if (!requireNamespace("monocle3", quietly = TRUE)) {
  cat("警告：monocle3 包未安装，跳过轨迹分析步骤\n")
} else {
  cat("Step 14: 开始 monocle3 轨迹分析...\n")

  # 创建轨迹分析目录
  trajectory_dir <- file.path(output_dir, "06_Trajectory_Analysis")

  # 提取表达矩阵和meta信息
  expr_matrix <- GetAssayData(scObject, assay = "RNA", layer = "counts")
  cell_metadata <- scObject@meta.data
  gene_metadata <- data.frame(gene_short_name = rownames(expr_matrix))
  rownames(gene_metadata) <- rownames(expr_matrix)

  # 构建monocle3的cds对象
  cds <- new_cell_data_set(expr_matrix,
                           cell_metadata = cell_metadata,
                           gene_metadata = gene_metadata)

  # 预处理
  cds <- preprocess_cds(cds, num_dim = 20)

  # ========== 关键修改：使用Seurat对象中已有的UMAP坐标 ==========
  # 这样可以确保轨迹图的细胞分布与前面的聚类注释图一致
  cat("使用Seurat对象中已有的UMAP坐标，确保与聚类注释图一致...\n")

  # 从Seurat对象提取UMAP坐标
  seurat_umap <- Embeddings(scObject, "umap")

  # 确保细胞顺序一致
  common_cells <- intersect(colnames(cds), rownames(seurat_umap))
  cat(sprintf("共有 %d 个细胞匹配\n", length(common_cells)))

  # 将Seurat的UMAP坐标转移到monocle3的cds对象中
  umap_coords_for_cds <- seurat_umap[common_cells, , drop = FALSE]
  colnames(umap_coords_for_cds) <- c("UMAP_1", "UMAP_2")

  # 设置UMAP降维结果到cds对象
  reducedDims(cds)$UMAP <- umap_coords_for_cds[colnames(cds), ]

  cat("已成功将Seurat UMAP坐标导入到monocle3 cds对象\n")

  # 聚类（使用较低分辨率）
  cds <- cluster_cells(cds, resolution = 0.0001)

  # 强制所有细胞在同一分区（确保连接所有簇）
  colData(cds)$partition <- as.factor(1)

  # 学习轨迹图
  cds <- learn_graph(cds,
                     use_partition = FALSE,
                     close_loop = FALSE,
                     learn_graph_control = list(
                       minimal_branch_len = 3,
                       ncenter = 2000,
                       geodesic_distance_ratio = 0.75,
                       euclidean_distance_ratio = 0.25
                     ))

  # 轨迹推断
  # 轨迹排序：自动选择root principal node，避免order_cells弹出交互窗口后报“No root node was chosen”
  get_root_pr_node <- function(cds, root_cells) {
    closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
    closest_vertex <- as.matrix(closest_vertex)

    if (is.null(rownames(closest_vertex))) {
      rownames(closest_vertex) <- colnames(cds)
    }

    root_cells <- intersect(root_cells, rownames(closest_vertex))
    if (length(root_cells) == 0) {
      root_cells <- colnames(cds)
    }

    root_vertex <- names(which.max(table(closest_vertex[root_cells, 1])))
    graph_node_names <- igraph::V(principal_graph(cds)[["UMAP"]])$name

    if (root_vertex %in% graph_node_names) {
      return(root_vertex)
    }

    root_vertex_index <- suppressWarnings(as.numeric(root_vertex))
    if (!is.na(root_vertex_index) && root_vertex_index >= 1 && root_vertex_index <= length(graph_node_names)) {
      return(graph_node_names[root_vertex_index])
    }

    graph_node_names[1]
  }

  # 可按需要在脚本前面自定义：
  # trajectory_root_cluster <- "0"
  # trajectory_root_celltype <- "Naive CD4 T cells"
  if (!exists("trajectory_root_cluster")) trajectory_root_cluster <- "0"
  if (!exists("trajectory_root_celltype")) trajectory_root_celltype <- NA_character_

  root_cells <- character(0)
  if (!is.na(trajectory_root_celltype) && "SingleR_celltype" %in% colnames(colData(cds))) {
    root_cells <- colnames(cds)[as.character(colData(cds)$SingleR_celltype) == trajectory_root_celltype]
  }
  if (length(root_cells) == 0 && "seurat_clusters" %in% colnames(colData(cds))) {
    root_cells <- colnames(cds)[as.character(colData(cds)$seurat_clusters) == as.character(trajectory_root_cluster)]
  }
  if (length(root_cells) == 0) {
    root_cells <- colnames(cds)
  }

  root_pr_node <- get_root_pr_node(cds, root_cells)
  cat(sprintf("monocle3 root principal node: %s (root cells: %d)\n", root_pr_node, length(root_cells)))

  cds <- tryCatch({
    order_cells(
      cds,
      reduction_method = "UMAP",
      root_pr_nodes = root_pr_node
    )
  }, error = function(e) {
    cat(sprintf("root_pr_nodes排序失败，改用root_cells排序。原因: %s\n", e$message))
    order_cells(
      cds,
      reduction_method = "UMAP",
      root_cells = root_cells[1]
    )
  })
  # 按pseudotime着色（去掉标签，只保留轨迹连线）
  pdf(file.path(trajectory_dir, "Trajectory_Pseudotime.pdf"), width=7, height=6)
  print(plot_cells(cds,
                   color_cells_by = "pseudotime",
                   label_groups_by_cluster = FALSE,
                   label_leaves = FALSE,
                   label_branch_points = FALSE,
                   label_roots = FALSE,
                   label_principal_points = FALSE,
                   show_trajectory_graph = TRUE,
                   trajectory_graph_color = "gray30",
                   trajectory_graph_segment_size = 1.0,
                   cell_size = 0.8))
  dev.off()

  # 按细胞类型着色（去掉标签，只保留轨迹连线）
  if("SingleR_celltype" %in% colnames(colData(cds))) {
    pdf(file.path(trajectory_dir, "Trajectory_CellType.pdf"), width=9, height=7)
    print(plot_cells(cds,
                     color_cells_by = "SingleR_celltype",
                     label_groups_by_cluster = FALSE,
                     label_leaves = FALSE,
                     label_branch_points = FALSE,
                     label_roots = FALSE,
                     label_principal_points = FALSE,
                     show_trajectory_graph = TRUE,
                     trajectory_graph_color = "gray30",
                     trajectory_graph_segment_size = 1.0,
                     cell_size = 0.8) +
            theme_classic(base_size = 14) +
            theme(
              legend.position = "right",
              legend.title = element_text(face = "bold", size = 12),
              legend.text = element_text(face = "bold", size = 10),
              axis.title = element_text(face = "bold", size = 14),
              axis.text = element_text(face = "bold", size = 12),
              plot.title = element_text(face = "bold", hjust = 0.5, size = 16)
            ) +
            ggtitle("Cellular Trajectory by Cell Type"))
    dev.off()
  }
  
  # 目标基因轨迹分析
  if (exists("cds") && !is.null(cds)) {
    for (gene in target_genes) {
      if (!(gene %in% rownames(cds))) {
        cat("警告：基因", gene, "不在cds对象中，跳过。\n")
        next
      }
      
      # 基因专属轨迹目录
      gene_traj_dir <- file.path(trajectory_dir, paste0(gene, "_Trajectory"))
      if (!dir.exists(gene_traj_dir)) dir.create(gene_traj_dir, recursive = TRUE)
      
      # 获取该基因在所有细胞的表达量
      expr_vec <- Matrix::t(assay(cds)[gene, ])
      expr_vec <- as.numeric(expr_vec)
      names(expr_vec) <- colnames(cds)
      
      # 按中位数分组
      med <- median(expr_vec)
      expr_group <- ifelse(expr_vec > med, "High", "Low")
      
      # 写入cds的colData
      colData(cds)[[paste0(gene, "_Group")]] <- factor(expr_group, levels = c("Low", "High"))
      
      # 轨迹图（去掉序号圈，只保留轨迹连线）
      pdf(file.path(gene_traj_dir, paste0("Trajectory_", gene, "_HighLow.pdf")), width=7, height=6)
      print(
        plot_cells(
          cds,
          color_cells_by = paste0(gene, "_Group"),
          show_trajectory_graph = TRUE,
          label_groups_by_cluster = FALSE,
          label_leaves = FALSE,
          label_branch_points = FALSE,
          label_roots = FALSE,
          label_principal_points = FALSE,
          trajectory_graph_color = "gray30",
          trajectory_graph_segment_size = 1.0,
          cell_size = 0.8,
          reduction_method = "UMAP"
        ) +
          scale_color_manual(values = c("Low" = "#079EDF", "High" = "#D377A9")) +
          ggtitle(paste("Trajectory -", gene, "(High/Low by median)"))
      )
      dev.off()
      
      # 伪时序密度分析
      pseudotime <- pseudotime(cds)
      pseudotime <- as.numeric(pseudotime)
      names(pseudotime) <- colnames(cds)
      
      plot_df <- data.frame(
        Pseudotime = pseudotime,
        Group = factor(expr_group, levels = c("Low", "High"))
      )
      plot_df <- plot_df[!is.na(plot_df$Pseudotime), ]
      
      p <- ggplot(plot_df, aes(x = Pseudotime, fill = Group, color = Group)) +
        geom_density(alpha = 0.35, size = 1.5, adjust = 1.1) +
        scale_fill_manual(values = c("Low" = "#079EDF", "High" = "#D377A9")) +
        scale_color_manual(values = c("Low" = "#079EDF", "High" = "#D377A9")) +
        labs(
          title = bquote(.(gene)~" density along pseudotime"),
          x = "Pseudotime",
          y = "Density",
          fill = NULL,
          color = NULL
        ) +
        theme_classic(base_size = 20) +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 22),
          axis.title = element_text(face = "bold", size = 20),
          axis.text = element_text(face = "bold", size = 16),
          legend.position = "top",
          legend.title = element_blank(),
          legend.text = element_text(face = "bold", size = 16)
        )
      
      pdf(file.path(gene_traj_dir, paste0("Pseudotime_Density_", gene, "_HighLow.pdf")), width=6, height=5)
      print(p)
      dev.off()
    }
  }
}

# ===============================================================================
# 第十一部分：数据导出
# ===============================================================================

cat("Step 15: 保存最终数据...\n")

# 保存Seurat对象
save(scObject, file = file.path(output_dir, "08_Data_Export", "scRNAseq_SeuratObject.RData"))

# 导出元数据
meta_data <- scObject@meta.data
write.csv(meta_data, file = file.path(output_dir, "08_Data_Export", "metadata.csv"))

# 导出UMAP坐标
umap_coords <- Embeddings(scObject, "umap")
write.csv(umap_coords, file = file.path(output_dir, "08_Data_Export", "umap_coordinates.csv"))

# 合并元数据和UMAP坐标
meta_umap <- cbind(meta_data, umap_coords)
write.csv(meta_umap, file = file.path(output_dir, "08_Data_Export", "MetaData_UMAP.csv"))

# ===============================================================================
# 分析完成
# ===============================================================================

# 保存分析参数到文件
analysis_parameters <- list(
  logFC_filter = logFC_filter,
  p_adj_filter = p_adj_filter,
  min_cells_gene = min_cells_gene,
  min_genes_per_cell = min_genes_per_cell,
  post_filter_cells = post_filter_cells,
  post_filter_mito = post_filter_mito,
  n_top_var_features = n_top_var_features,
  n_pcs = n_pcs,
  n_topmarker_heat = n_topmarker_heat,
  cluster_resolution = cluster_resolution,
  neighbor_dims = neighbor_dims,
  target_genes = target_genes,
  important_marker_genes = important_marker_genes,
  ko_n_hvg = ko_n_hvg,
  ko_qc_mtThreshold = ko_qc_mtThreshold,
  ko_qc_minLSize = ko_qc_minLSize,
  ko_nc_nNet = ko_nc_nNet,
  ko_nc_nCells = ko_nc_nCells,
  ko_pval_threshold = ko_pval_threshold
)

# 保存参数到JSON格式文件
library(jsonlite)
write_json(analysis_parameters, file.path(output_dir, "analysis_parameters.json"), pretty = TRUE)

cat("分析参数已保存至：analysis_parameters.json\n")

# ===============================================================================
# 第十二部分：单细胞基因敲除模拟分析 (scTenifoldKnk)
# Single Cell Gene Knockout Simulation Analysis
# 生成类似文献图4 E/F风格的散点图和柱状图
# ===============================================================================

cat("\n================================================================================\n")
cat("         Step 16: 单细胞基因敲除模拟分析 (scTenifoldKnk)                        \n")
cat("         Single Cell Gene Knockout Simulation Analysis                         \n")
cat("================================================================================\n")

# 加载基因敲除模拟分析所需的包
if (!requireNamespace("scTenifoldKnk", quietly = TRUE)) {
  cat("警告：scTenifoldKnk 包未安装，跳过基因敲除模拟分析\n")
  cat("安装方法：devtools::install_github('cailab-tamu/scTenifoldKnk')\n")
} else {

  library(scTenifoldKnk)
  library(ggrepel)
  set.seed(123)

  # =============================================
  # 敲除分析参数设置（使用前面定义的参数）
  # =============================================

  knockout_params <- list(
    # 敲除目标基因（使用target_genes）
    knockout_genes = target_genes,

    # 用于构建网络的高变基因数量（已优化）
    n_hvg = ko_n_hvg,

    # scTenifoldKnk特定参数（已优化以加速）
    qc_mtThreshold = ko_qc_mtThreshold,
    qc_minLSize = ko_qc_minLSize,
    nc_nNet = ko_nc_nNet,
    nc_nCells = ko_nc_nCells,

    # 显著性阈值
    pval_threshold = ko_pval_threshold
  )

  cat("\n--- 敲除分析参数设置 ---\n")
  cat(sprintf("敲除目标基因: %s\n", paste(knockout_params$knockout_genes, collapse = ", ")))
  cat(sprintf("用于网络构建的高变基因数: %d\n", knockout_params$n_hvg))

  # =============================================
  # 创建输出目录
  # =============================================

  knockout_output_dir <- file.path(output_dir, "10_Gene_Knockout_Analysis")
  if (!dir.exists(knockout_output_dir)) {
    dir.create(knockout_output_dir, recursive = TRUE)
    cat(sprintf("创建敲除分析目录: %s\n", knockout_output_dir))
  }

  # =============================================
  # 对每个诊断基因进行敲除分析
  # =============================================

  # 检查可用基因
  available_knockout_genes <- knockout_params$knockout_genes[knockout_params$knockout_genes %in% rownames(scObject)]
  missing_knockout_genes <- setdiff(knockout_params$knockout_genes, available_knockout_genes)

  if (length(missing_knockout_genes) > 0) {
    cat(sprintf("\n警告: 以下基因不在数据中: %s\n", paste(missing_knockout_genes, collapse = ", ")))
  }

  if (length(available_knockout_genes) > 0) {

    cat(sprintf("\n========================================================\n"))
    cat(sprintf("开始对 %d 个诊断基因进行敲除模拟分析\n", length(available_knockout_genes)))
    cat(sprintf("========================================================\n"))

    # 提取表达矩阵
    countMat <- GetAssayData(scObject, layer = "counts")

    # 提取高可变基因
    scObject <- FindVariableFeatures(object = scObject, selection.method = "vst", nfeatures = knockout_params$n_hvg)
    hvgs <- VariableFeatures(scObject)

    # 遍历每个敲除基因
    for (target_gene in available_knockout_genes) {

      cat(sprintf("\n########## 分析基因: %s ##########\n", target_gene))

      # 创建基因专属目录
      gene_output_dir <- file.path(knockout_output_dir, paste0("Gene_", target_gene))
      if (!dir.exists(gene_output_dir)) {
        dir.create(gene_output_dir, recursive = TRUE)
      }

      tryCatch({
        # 准备数据（包含目标基因和高变基因）
        data <- as.data.frame(countMat[unique(c(target_gene, hvgs)), ])

        cat(sprintf("  表达矩阵: %d 基因 x %d 细胞\n", nrow(data), ncol(data)))

        # 执行虚拟敲除（使用参考.R的参数）
        cat("  运行scTenifoldKnk分析...\n")

        result <- scTenifoldKnk(
          countMatrix = data,
          gKO = target_gene,
          qc_mtThreshold = knockout_params$qc_mtThreshold,
          qc_minLSize = knockout_params$qc_minLSize,
          nc_nNet = knockout_params$nc_nNet,
          nc_nCells = knockout_params$nc_nCells
        )

        cat("  scTenifoldKnk分析完成!\n")

        # 输出差异分析的结果
        df <- result$diffRegulation
        df <- df[df$gene != target_gene, ]

        # 保存显著差异基因（CSV格式）
        outTab <- df[df$p.adj < knockout_params$pval_threshold, ]
        write.csv(outTab, file = file.path(gene_output_dir, paste0(target_gene, "_sigDiff.csv")),
                  row.names = FALSE)

        cat(sprintf("  发现 %d 个显著差异调控基因 (p.adj < 0.05)\n", nrow(outTab)))

        # =============================================
        # 散点图（Scatter Plot）
        # =============================================
        cat("  生成散点图...\n")

        # 准备散点图数据
        df$log_p.adj <- -log10(df$p.adj)
        df$significant <- ifelse(df$p.adj < knockout_params$pval_threshold, "Significant", "Not significant")
        label_genes <- subset(df, p.adj < knockout_params$pval_threshold)

        # 设置y轴上限（避免极端值）
        y_upper <- quantile(df$log_p.adj, 0.999, na.rm = TRUE)

        # 绑制散点图（使用Set3配色）
        scatter_colors <- colorRampPalette(brewer.pal(12, "Set3"))(2)
        p_scatter <- ggplot(df, aes(x = Z, y = log_p.adj, color = significant)) +
          geom_point(alpha = 0.7, size = 1.5) +
          scale_color_manual(values = c("Significant" = scatter_colors[1], "Not significant" = "gray70")) +
          geom_hline(yintercept = -log10(knockout_params$pval_threshold), linetype = "dashed", color = scatter_colors[1]) +
          geom_text_repel(data = label_genes, aes(label = gene), size = 3, max.overlaps = 50,
                          color = "black", fontface = "italic") +
          labs(title = paste0(target_gene, " Knockout"),
               x = "Z-score",
               y = "-log10(p.adj)") +
          theme_classic(base_size = 14) +
          coord_cartesian(ylim = c(0, y_upper)) +
          theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
            axis.title = element_text(size = 14, face = "bold"),
            axis.text = element_text(size = 12),
            legend.position = "none"
          )

        # 保存散点图
        pdf(file.path(gene_output_dir, paste0(target_gene, "_KO_scatter.pdf")), width = 6, height = 5)
        print(p_scatter)
        dev.off()
        cat("  已保存: ", target_gene, "_KO_scatter.pdf\n")

        # =============================================
        # 柱状图（Bar Plot）- Top 20基因
        # =============================================
        cat("  生成柱状图...\n")

        # 获取Top 20差异调控基因（按FC排序）
        top_genes <- head(df[order(-df$FC), ], 20)

        # 使用Set3配色为每个基因分配颜色
        bar_colors <- colorRampPalette(brewer.pal(12, "Set3"))(nrow(top_genes))

        # 绑制柱状图（使用Set3配色）
        p_bar <- ggplot(top_genes, aes(x = reorder(gene, FC), y = FC, fill = reorder(gene, FC))) +
          geom_bar(stat = 'identity', alpha = 0.9) +
          scale_fill_manual(values = bar_colors) +
          coord_flip() +
          labs(title = paste0("Top 20 Differentially Regulated Genes\n(", target_gene, " Knockout)"),
               x = "Gene",
               y = "FC") +
          theme_classic(base_size = 14) +
          theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
            axis.title = element_text(size = 13, face = "bold"),
            axis.text.y = element_text(size = 11, face = "italic"),
            axis.text.x = element_text(size = 11),
            legend.position = "none"
          )

        # 保存柱状图
        pdf(file.path(gene_output_dir, paste0(target_gene, "_KO_barplot.pdf")), width = 6, height = 5)
        print(p_bar)
        dev.off()
        cat("  已保存: ", target_gene, "_KO_barplot.pdf\n")

        # =============================================
        # 保存完整结果
        # =============================================
        write.csv(df, file = file.path(gene_output_dir, paste0(target_gene, "_KO_allResults.csv")),
                  row.names = FALSE)

        cat(sprintf("  基因 %s 敲除分析完成!\n", target_gene))

      }, error = function(e) {
        cat(sprintf("  基因 %s 敲除分析失败: %s\n", target_gene, e$message))
      })
    }

  } else {
    cat("警告: 没有有效的敲除基因!\n")
  }

  cat("\n================================================================================\n")
  cat("         基因敲除模拟分析完成!                                                  \n")
  cat("================================================================================\n")
  cat(sprintf("\n结果保存至: %s\n", knockout_output_dir))

  cat("\n=== 敲除分析输出结构 ===\n")
  cat("10_Gene_Knockout_Analysis/\n")
  cat("  └── Gene_[GeneName]/\n")
  cat("      ├── [Gene]_KO_scatter.pdf      (散点图)\n")
  cat("      ├── [Gene]_KO_barplot.pdf      (柱状图)\n")
  cat("      ├── [Gene]_sigDiff.csv         (显著差异基因)\n")
  cat("      └── [Gene]_KO_allResults.csv   (完整结果)\n")

}

cat("\n================================================================================\n")
cat("                       所有分析已成功完成!                                       \n")
cat("================================================================================\n")
