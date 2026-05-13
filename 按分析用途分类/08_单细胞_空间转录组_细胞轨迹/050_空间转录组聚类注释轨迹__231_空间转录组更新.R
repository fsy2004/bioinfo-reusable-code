# ==========================================================================
# 脚本名     : 空间转录组聚类注释轨迹.R
# 分类       : 08_单细胞_空间转录组_细胞轨迹
# 项目来源   : 从压缩包 231.空间转录组更新.rar 整理
# 原始文件   : 231.空间转录组更新\空间转录组.R
# 用途       : 10x Visium/Seurat 空间转录组分析：质控、空间聚类、空间 marker、空间可变基因、cluster-based 细胞类型注释和 monocle3 拟时序。
# 结果图     : 空间QC图；空间聚类图；单cluster空间图；Top marker空间FeaturePlot；小提琴图；空间可变基因图；UMAP细胞类型图；空间细胞类型注释图；轨迹总览图；空间pseudotime图；轨迹一致性图
# 非肿瘤消化适配: 适合。可用于炎症性肠病、肝病、胰腺/胃肠组织等非肿瘤空间转录组。
# 主要 R 包  : Seurat; ggplot2; patchwork; dplyr; monocle3; RColorBrewer
# 整理日期   : 2026-05-13
# 备注       : 保留原始代码逻辑，仅添加统一说明头；运行前请把 workDir/setwd 和输入文件名改成当前项目路径。
# ==========================================================================


# 加载所需R包
library(Seurat)        # 单细胞分析主包
library(ggplot2)       # 绘图
library(patchwork)     # 多图拼接
library(dplyr)         # 数据处理
library(Rfast2)        # 快速数据处理
library(hdf5r)         # 读写h5文件
library(viridis)
library(glmGamPoi)
library(SingleR)
library(celldex)
library(RColorBrewer)
library(monocle3)
library(SeuratWrappers)
# 自动生成足够的颜色
library(viridis)
# 用户自定义参数
project_id <- "BRCA"    # 项目/样本名称
work_dir <- "H:\\常用分析生信\\231.空间转录组更新"  # 工作目录
output_dir <- "SpatialAnalysisResults"  # 输出文件夹

# 创建各类结果的子文件夹
qc_dir <- file.path(output_dir, "01_Quality_Control")
cluster_dir <- file.path(output_dir, "02_Clustering_Analysis")
de_dir <- file.path(output_dir, "03_Differential_Expression")
spatial_dir <- file.path(output_dir, "04_Spatial_Features")
annotation_dir <- file.path(output_dir, "05_Cell_Type_Annotation")
trajectory_dir <- file.path(output_dir, "06_Trajectory_Analysis")
data_dir <- file.path(output_dir, "07_Data_Export")
raw_h5_file <- "GSM6177599_NYU_BRCA0_Vis_processed_filtered_feature_bc_matrix.h5" # 原始表达矩阵文件
genes_to_plot <- c("CENPA", "KIF2C", "EZH2")  # 可添加更多基因
logfc_cutoff <- 1      # logFC阈值
adjp_cutoff <- 0.05    # 调整后p值阈值

# =======================
# 步骤1：设置工作目录和输出目录
# =======================
setwd(work_dir)  # 设置当前工作目录
if (!dir.exists(output_dir)) { dir.create(output_dir) }  # 若输出文件夹不存在则新建

# 创建所有子文件夹
for (subdir in c(qc_dir, cluster_dir, de_dir, spatial_dir, annotation_dir, trajectory_dir, data_dir)) {
  if (!dir.exists(subdir)) { dir.create(subdir, recursive = TRUE) }
}

cat("Step 1/10: 文件重命名为标准格式...\n")  # 进度提示

# 文件重命名，标准化10x空间转录组文件名
# 10x Genomics空间转录组原始文件，重命名为Seurat等分析工具默认识别的标准文件名
file.rename("GSM6177599_NYU_BRCA0_Vis_processed_spatial_tissue_hires_image.png", "tissue_hires_image.png")      # 高分辨率组织图片
file.rename("GSM6177599_NYU_BRCA0_Vis_processed_spatial_tissue_lowres_image.png", "tissue_lowres_image.png")    # 低分辨率组织图片
file.rename("GSM6177599_NYU_BRCA0_Vis_processed_spatial_scalefactors_json.json", "scalefactors_json.json")      # 空间缩放因子
file.rename("GSM6177599_NYU_BRCA0_Vis_processed_spatial_tissue_positions_list.csv", "tissue_positions_list.csv")# 组织空间坐标

cat("Step 2/10: 读取表达矩阵...\n")  # 进度提示
# 读取表达矩阵
expr_matrix <- Read10X_h5(raw_h5_file)  # 读取h5格式表达矩阵

cat("Step 3/10: 读取空间图片和空间信息...\n")  # 进度提示
# 读取空间图片和空间信息
spatial_img <- Read10X_Image(
  image.dir = ".",      # 当前目录
  filter.matrix = TRUE  # 过滤矩阵
)

cat("Step 4/10: 创建Seurat对象并挂载图片...\n")  # 进度提示
# 创建Seurat对象
sp_obj <- CreateSeuratObject(counts = expr_matrix, assay = "Spatial", project = project_id)  # 创建对象
spatial_img <- spatial_img[Cells(x = sp_obj)]  # 匹配细胞
DefaultAssay(spatial_img) <- "Spatial"         # 设置默认assay
sp_obj[["spatial_slice"]] <- spatial_img       # 挂载图片到对象

cat("Step 5/10: 预处理可视化...\n")  # 进度提示
# 获取分组数量（如按cluster分组，否则为1）
group_count <- length(unique(Idents(sp_obj)))

my_colors <- viridis::viridis(group_count, option = "C")  # 你可以换成"A"/"B"/"D"等

# 绘制nCount_Spatial的小提琴图（无点，窄宽度，自动多色）
vln_plot <- VlnPlot(
  sp_obj,
  features = "nCount_Spatial",
  pt.size = 0,                        # 不显示点
  cols = my_colors                    # 自动分配颜色
) +
  geom_boxplot(width = 0.08, outlier.shape = NA, fill = "white", color = "black", alpha = 0.3) +  # 窄箱线
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black", size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.position = "none"
  ) +
  labs(
    title = "nCount_Spatial Distribution",
    x = "Group",
    y = "nCount_Spatial"
  )

spatial_plot <- SpatialFeaturePlot(sp_obj, features = "nCount_Spatial") + theme(legend.position = "right")  # 空间分布

# 分别输出两个独立的图到质量控制文件夹
pdf(file = file.path(qc_dir, "FeatureDistribution_violin.pdf"), width = 5, height = 6)
print(vln_plot)
dev.off()

pdf(file = file.path(qc_dir, "FeatureDistribution_spatial.pdf"), width = 9, height = 6)
print(spatial_plot)
dev.off()

cat("Step 6/10: 数据标准化...\n")  # 进度提示
# SCTransform标准化
sp_obj <- SCTransform(sp_obj, assay = "Spatial", verbose = FALSE)  # 标准化

cat("Step 7/10: 降维与聚类...\n")  # 进度提示
# PCA降维
sp_obj <- RunPCA(sp_obj, assay = "SCT", verbose = FALSE)  # PCA
sp_obj <- FindNeighbors(sp_obj, reduction = "pca", dims = 1:30)  # 邻居查找

# 使用固定分辨率进行聚类
resolution <- 0.6
sp_obj <- FindClusters(sp_obj, resolution = resolution, verbose = FALSE)  # 聚类
sp_obj$seurat_clusters <- Idents(sp_obj)

n_clusters <- length(unique(sp_obj$seurat_clusters))
cat(sprintf("使用分辨率 %.1f，共有 %d 个聚类\n", resolution, n_clusters))

sp_obj <- RunUMAP(sp_obj, reduction = "pca", dims = 1:30)  # UMAP降维

cat("Step 8/10: 聚类与空间分布可视化...\n")  # 进度提示
# 聚类和空间分布可视化
umap_plot <- DimPlot(sp_obj, reduction = "umap", label = TRUE)  # UMAP聚类
spatial_cluster_plot <- SpatialDimPlot(sp_obj, label = TRUE, label.size = 3)  # 空间聚类

# 分别输出两个独立的图到聚类分析文件夹
pdf(file = file.path(cluster_dir, "ClusterDistribution_umap.pdf"), width = 8, height = 6)
print(umap_plot)
dev.off()

pdf(file = file.path(cluster_dir, "ClusterDistribution_spatial.pdf"), width = 8, height = 6)
print(spatial_cluster_plot)
dev.off()

# 每个聚类单独空间分布
pdf(file = file.path(cluster_dir, "ClusterSingleSpatial.pdf"), width = 12, height = 8)
SpatialDimPlot(
  sp_obj,
  cells.highlight = CellsByIdentities(object = sp_obj, idents = levels(sp_obj$seurat_clusters)),
  facet.highlight = TRUE,
  ncol = 4
)
dev.off()
# 新建用于保存每个聚类单独空间分布图的文件夹
single_cluster_dir <- file.path(cluster_dir, "Individual_Clusters")
if (!dir.exists(single_cluster_dir)) dir.create(single_cluster_dir)

# 获取所有聚类编号
cluster_ids <- levels(sp_obj$seurat_clusters)

# 循环输出每个聚类的空间分布图
for (cid in cluster_ids) {
  # 生成当前聚类的空间分布图
  p <- SpatialDimPlot(
    sp_obj,
    cells.highlight = CellsByIdentities(object = sp_obj, idents = cid),
    facet.highlight = FALSE,
    cols.highlight = c("orange", "grey90"),  # 橙色高亮，灰色为背景
    label = TRUE,
    label.size = 4
  ) +
    ggtitle(paste("Spatial Distribution - Cluster", cid)) +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  
  # 保存为PDF
  pdf(file = file.path(single_cluster_dir, paste0("Cluster_", cid, "_Spatial.pdf")), width = 8, height = 8)
  print(p)
  dev.off()
}


cat("Step 9/10: 差异表达分析...\n")  # 进度提示
# 差异表达分析
all_markers <- FindAllMarkers(
  object = sp_obj,
  only.pos = FALSE,
  min.pct = 0.1,
  logfc.threshold = logfc_cutoff
)
# 按阈值筛选marker
sig_markers <- all_markers[
  (abs(as.numeric(as.vector(all_markers$avg_log2FC))) > logfc_cutoff &
     as.numeric(as.vector(all_markers$p_val_adj)) < adjp_cutoff),
]
write.table(sig_markers, file = file.path(de_dir, "SignificantMarkers.csv"), sep = ",", row.names = FALSE, quote = FALSE)

# 可视化部分marker基因
pdf(file = file.path(de_dir, "TopMarkersSpatial.pdf"), width = 10, height = 8)
SpatialFeaturePlot(object = sp_obj, features = rownames(all_markers)[1:6], alpha = c(0.1, 1), ncol = 3)
dev.off()
cat("Step 5/10: 预处理可视化...\n")  # 进度提示

# 获取分组数量（如按cluster分组，否则为1）
group_count <- length(unique(Idents(sp_obj)))


my_colors <- viridis::viridis(group_count, option = "C")  # 你可以换成"A"/"B"/"D"等

# 绘制nCount_Spatial的小提琴图（无点，窄宽度，自动多色）
vln_plot <- VlnPlot(
  sp_obj,
  features = "nCount_Spatial",
  pt.size = 0,                        # 不显示点
  cols = my_colors                    # 自动分配颜色
) +
  geom_boxplot(width = 0.08, outlier.shape = NA, fill = "white", color = "black", alpha = 0.3) +  # 窄箱线
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black", size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.position = "none"
  ) +
  labs(
    title = "Spatial Transcript Count Distribution",
    x = "Group",
    y = "nCount_Spatial"
  )

pdf(file = file.path(qc_dir, "FeatureDistribution_nCount_Spatial.pdf"), width = 7, height = 5)  # 输出PDF
print(vln_plot)
dev.off()  # 关闭设备
cat("Step 10/10: 空间可变基因分析与保存...\n")  # 进度提示
# 空间可变基因分析
sp_obj <- NormalizeData(sp_obj, assay = "Spatial")  # 归一化
sp_obj <- FindVariableFeatures(sp_obj, assay = "Spatial", selection.method = "vst", nfeatures = 2000)  # 高变基因
top_var_genes <- VariableFeatures(sp_obj, assay = "Spatial")[1:100]  # 前100高变基因
sp_obj <- ScaleData(sp_obj, assay = "Spatial", features = top_var_genes)  # 标准化
sp_obj <- FindSpatiallyVariableFeatures(
  sp_obj,
  assay = "Spatial",
  features = top_var_genes,
  selection.method = "moransi"
)

# 提取空间可变基因（新版Seurat直接用VariableFeatures获取）
top_spatial_genes <- head(VariableFeatures(sp_obj, assay = "Spatial"), 6)

# 可视化空间可变基因
pdf(file = file.path(spatial_dir, "SpatiallyVariableGenes.pdf"), width = 10, height = 8)
SpatialFeaturePlot(sp_obj, features = top_spatial_genes, ncol = 3, alpha = c(0.1, 1))
dev.off()

for (gene in genes_to_plot) {
  # 空间聚类分布图
  p_spatial <- SpatialDimPlot(sp_obj, label = TRUE, label.size = 3)
  # 指定基因空间表达图
  p_gene <- SpatialFeaturePlot(sp_obj, features = gene, alpha = c(0.1, 1))
  # 输出PDF，每个基因一个文件
  pdf(file = file.path(spatial_dir, paste0(gene, "_Spatial.pdf")), width = 10, height = 8)
  print(p_spatial + p_gene)
  dev.off()
}
# 保存Seurat对象
save(sp_obj, file = file.path(data_dir, "SpatialSeuratObject.RData"))  # 保存对象
# 假设 sp_obj 是 Seurat 对象
expr_matrix <- as.data.frame(GetAssayData(sp_obj, slot = "counts"))  # 或 "data" 获取归一化数据
write.csv(expr_matrix, file = file.path(data_dir, "expression_matrix.csv"))
meta_data <- sp_obj@meta.data
write.csv(meta_data, file = file.path(data_dir, "metadata.csv"))
umap_coords <- Embeddings(sp_obj, "umap")
write.csv(umap_coords, file = file.path(data_dir, "umap_coordinates.csv"))

# 3. 空间可变基因分布图
sp_obj <- FindSpatiallyVariableFeatures(sp_obj, assay = "Spatial", selection.method = "moransi")
top_spatial_genes <- head(VariableFeatures(sp_obj, assay = "Spatial"), 6)
pdf(file = file.path(spatial_dir, "SpatialFeaturePlot_TopSpatialGenes.pdf"), width = 12, height = 8)
print(SpatialFeaturePlot(sp_obj, features = top_spatial_genes, ncol = 3))
dev.off()


# 4. 指定基因空间表达图

for (gene in genes_to_plot) {
  pdf(file = file.path(spatial_dir, paste0("SpatialFeaturePlot_", gene, ".pdf")), width = 8, height = 6)
  print(SpatialFeaturePlot(sp_obj, features = gene))
  dev.off()
}


# 2. 获取参考数据集（人类用HumanPrimaryCellAtlasData，小鼠用MouseRNAseqData）
ref <- celldex::HumanPrimaryCellAtlasData()
# ref <- celldex::MouseRNAseqData()  # 如果是小鼠数据，取消注释本行

# 3. 提取表达矩阵（用归一化数据，适合SingleR）
expr <- GetAssayData(sp_obj, slot = "data")

# 4. 运行SingleR自动注释
singleR_result <- SingleR(
  test = expr,
  ref = ref,
  labels = ref$label.main
)

# 5. 将注释结果加入Seurat对象
sp_obj$celltype_original <- singleR_result$labels

# ========== 基于聚类的细胞类型注释 ==========
cat("基于聚类进行细胞类型注释...\n")

# 创建聚类与细胞类型的交叉表（直接使用原始SingleR注释）
cluster_celltype_table <- table(sp_obj$seurat_clusters, sp_obj$celltype_original)
cat("聚类与细胞类型交叉表：\n")
print(cluster_celltype_table)

# 为每个聚类找到主要的细胞类型
cluster_annotation <- data.frame(
  cluster = rownames(cluster_celltype_table),
  dominant_celltype = apply(cluster_celltype_table, 1, function(x) names(which.max(x))),
  max_count = apply(cluster_celltype_table, 1, max),
  total_count = rowSums(cluster_celltype_table),
  stringsAsFactors = FALSE
)

# 计算每个聚类的主导细胞类型占比
cluster_annotation$proportion <- cluster_annotation$max_count / cluster_annotation$total_count

# 打印聚类注释结果
cat("聚类注释结果：\n")
print(cluster_annotation)

# 检查每个聚类的主导性
cat("\n聚类质量评估：\n")
for (i in 1:nrow(cluster_annotation)) {
  cluster_id <- cluster_annotation$cluster[i]
  celltype <- cluster_annotation$dominant_celltype[i]
  prop <- round(cluster_annotation$proportion[i], 3)
  
  quality <- if (prop >= 0.7) "优秀" else if (prop >= 0.5) "良好" else if (prop >= 0.3) "一般" else "差"
  
  cat(sprintf("聚类 %s: %s (%.1f%%) - 质量: %s\n", 
              cluster_id, celltype, prop*100, quality))
}

# 基于聚类的细胞类型注释
sp_obj$celltype_cluster_based <- cluster_annotation$dominant_celltype[match(sp_obj$seurat_clusters, cluster_annotation$cluster)]

# 保存详细的聚类分析结果
cluster_analysis_dir <- file.path(annotation_dir, "Cluster_Analysis")
if (!dir.exists(cluster_analysis_dir)) dir.create(cluster_analysis_dir)

# 保存聚类注释映射表
write.csv(cluster_annotation, file = file.path(cluster_analysis_dir, "Cluster_CellType_Mapping.csv"), row.names = FALSE)

# 保存交叉表
write.csv(as.data.frame.matrix(cluster_celltype_table), file = file.path(cluster_analysis_dir, "Cluster_CellType_CrossTable.csv"))

# 6. 绘制基于聚类的细胞类型分布图

# 获取基于聚类的细胞类型
celltype_levels <- unique(sp_obj$celltype_cluster_based)
n_celltypes <- length(celltype_levels)

# 使用固定的颜色方案，确保可重复性
my_colors <- RColorBrewer::brewer.pal(min(n_celltypes, 8), "Set2")
if (n_celltypes > 8) {
  my_colors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(n_celltypes)
}

# UMAP图：基于聚类的细胞类型注释
umap_celltype_plot <- DimPlot(
  sp_obj,
  reduction = "umap",
  group.by = "celltype_cluster_based",
  label = TRUE,
  label.size = 5,
  repel = TRUE,
  cols = my_colors
) +
  ggtitle("Cell Type UMAP (Cluster-based Annotation)") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  )

pdf(file = file.path(annotation_dir, "UMAP_celltype_cluster_based.pdf"), width = 12, height = 9)
print(umap_celltype_plot)
dev.off()

# 比较图：原始SingleR注释 vs 基于聚类的注释
comparison_plot1 <- DimPlot(sp_obj, reduction = "umap", group.by = "celltype_original", label = TRUE, label.size = 3) + 
  ggtitle("Original SingleR Annotation") + theme(legend.position = "none")

comparison_plot2 <- DimPlot(sp_obj, reduction = "umap", group.by = "celltype_cluster_based", label = TRUE, label.size = 3) + 
  ggtitle("Cluster-based Annotation") + theme(legend.position = "none")

pdf(file = file.path(annotation_dir, "CellType_Annotation_Comparison.pdf"), width = 16, height = 8)
print(comparison_plot1 + comparison_plot2)
dev.off()

# 保存基于聚类的细胞类型统计
celltype_count_cluster_based <- as.data.frame(table(sp_obj$celltype_cluster_based))
colnames(celltype_count_cluster_based) <- c("CellType", "Count")
write.csv(celltype_count_cluster_based, file = file.path(annotation_dir, "CellType_Count_Cluster_Based.csv"), row.names = FALSE)

# 保存原始细胞类型统计
celltype_count_original <- as.data.frame(table(sp_obj$celltype_original))
colnames(celltype_count_original) <- c("CellType", "Count")
write.csv(celltype_count_original, file = file.path(annotation_dir, "CellType_Count_Original.csv"), row.names = FALSE)



# 7. 空间分布图：基于聚类的细胞类型注释

# 1. 确保celltype_cluster_based在meta.data中且无NA/空
sp_obj@meta.data$celltype_cluster_based <- as.character(sp_obj@meta.data$celltype_cluster_based)
sp_obj@meta.data$celltype_cluster_based[is.na(sp_obj@meta.data$celltype_cluster_based) | sp_obj@meta.data$celltype_cluster_based == ""] <- "Unknown"
sp_obj@meta.data$celltype_cluster_based <- factor(sp_obj@meta.data$celltype_cluster_based)

# 2. 设置基于聚类的细胞类型为当前分组
Idents(sp_obj) <- sp_obj@meta.data$celltype_cluster_based

# 3. 绘制空间分布图
spatial_celltype_plot <- SpatialDimPlot(
  sp_obj,
  label = TRUE,
  label.size = 4,
  repel = TRUE
) +
  ggtitle("Spatial Cell Type Annotation (Cluster-based)") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  )

# 4. 导出为PDF
pdf(file = file.path(annotation_dir, "Spatial_CellType_Annotation_Cluster_Based.pdf"), width = 14, height = 11)
print(spatial_celltype_plot)
dev.off()

# 保存原始细胞类型统计（用于比较）
# 已在上面保存，此处移除重复代码

# 保存合并后的细胞类型统计
# 已移除celltype_consolidated，此处删除相关代码
# 移除重复的细胞类型统计保存代码（已在上面保存过）
split_markers <- split(all_markers, all_markers$cluster)
for (cid in names(split_markers)) {
  write.csv(split_markers[[cid]], file = file.path(de_dir, paste0("Markers_Cluster_", cid, ".csv")), row.names = FALSE)
}
spatial_var_genes <- VariableFeatures(sp_obj, assay = "Spatial")
write.csv(spatial_var_genes, file = file.path(spatial_dir, "SpatiallyVariableGenes.csv"), row.names = FALSE)
meta_data <- sp_obj@meta.data
umap_coords <- Embeddings(sp_obj, "umap")
meta_umap <- cbind(meta_data, umap_coords)
write.csv(meta_umap, file = file.path(data_dir, "MetaData_UMAP.csv"))

cat("分析流程全部完成，结果已保存至：", output_dir, "\n")  # 结束提示

# ========== 空间拟时序轨迹分析完整流程 ==========
# ========== 空间拟时序轨迹分析完整流程 ==========

cat("开始空间轨迹分析...\n")

# 1. 转换Seurat对象为CDS对象（确保聚类信息正确传递）
cds <- as.cell_data_set(sp_obj)

# 重要：手动添加Seurat聚类信息到CDS对象的colData中
# 替代方案：直接在colData中检查和
# 确保seurat_clusters是因子类型
colData(cds)$seurat_clusters <- as.factor(sp_obj@meta.data$seurat_clusters)

# 或者使用不同的列名
colData(cds)$seurat_cluster_labels <- paste0("Cluster_", sp_obj@meta.data$seurat_clusters)

# 然后在plot_cells中使用
p2 <- plot_cells(cds, 
                 color_cells_by = "seurat_cluster_labels",
                 label_cell_groups = TRUE,
                 label_groups_by_cluster = TRUE)


# 2. Monocle3聚类和轨迹构建
cds <- cluster_cells(cds)
cds <- learn_graph(cds)

# 3. 选择起始点并计算拟时序
cds <- order_cells(cds)

# 4. 将pseudotime信息添加回Seurat对象
sp_obj$pseudotime <- pseudotime(cds)
sp_obj$monocle3_clusters <- clusters(cds)
# ========== 轨迹可视化面板 ==========

pdf(file = file.path(trajectory_dir, "Trajectory_Overview_Panel_Fixed.pdf"), width = 16, height = 12)

# p1: 按伪时序着色的轨迹
p1 <- plot_cells(cds, 
                 color_cells_by = "pseudotime", 
                 label_cell_groups = FALSE,
                 show_trajectory_graph = TRUE) + 
  ggtitle("Trajectory by Pseudotime") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

# p2: 直接使用Trajectory_Seurat_Clusters_Only的成功配置
p2 <- plot_cells(cds,
                 color_cells_by = "seurat_clusters",
                 label_cell_groups = TRUE,           # 显示标签
                 label_groups_by_cluster = TRUE,     # 按聚类显示标签
                 show_trajectory_graph = TRUE,       # 显示轨迹图
                 graph_label_size = 4,               # 标签大小
                 cell_size = 0.5) +                 # 细胞点大小
  ggtitle("Trajectory by Seurat Clusters") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

# p3: 按Monocle3 partition着色的轨迹
p3 <- plot_cells(cds, 
                 color_cells_by = "partition", 
                 label_groups_by_cluster = TRUE,
                 show_trajectory_graph = TRUE) + 
  ggtitle("Trajectory by Partition") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

# p4: 空间伪时序分布
p4 <- SpatialFeaturePlot(sp_obj, features = "pseudotime") + 
  ggtitle("Spatial Pseudotime") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

print((p1 + p2) / (p3 + p4))
dev.off()

cat("轨迹Overview面板已生成！\n")


# ========== 额外：验证Seurat聚类是否正确传递 ==========

# 检查CDS对象中是否包含seurat_clusters信息
cat("CDS对象中的colData列名：", colnames(colData(cds)), "\n")
cat("seurat_clusters的唯一值：", unique(colData(cds)$seurat_clusters), "\n")

# ========== 单独生成带Seurat cluster标签的轨迹图 ==========

pdf(file = file.path(trajectory_dir, "Trajectory_Seurat_Clusters_Only.pdf"), width = 10, height = 8)
plot_cells(cds,
           color_cells_by = "seurat_clusters",
           label_cell_groups = TRUE,           # 显示标签
           label_groups_by_cluster = TRUE,     # 按聚类显示标签
           show_trajectory_graph = TRUE,       # 显示轨迹图
           graph_label_size = 4,               # 标签大小
           cell_size = 0.5) +                 # 细胞点大小
  ggtitle("Trajectory with Seurat Cluster Labels") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10))
dev.off()

cat("轨迹分析完成！\n")


# ========== 基础轨迹可视化 ==========

# 1.1 多种轨迹图组合
pdf(file = file.path(trajectory_dir, "Trajectory_Overview_Panel.pdf"), width = 16, height = 12)
p1 <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE) + 
  ggtitle("Trajectory by Pseudotime")
p2 <- plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = TRUE) + 
  ggtitle("Trajectory by Cluster")

p3 <- plot_cells(cds, color_cells_by = "partition", label_groups_by_cluster = TRUE) + 
  ggtitle("Trajectory by Partition")
p4 <- SpatialFeaturePlot(sp_obj, features = "pseudotime") + 
  ggtitle("Spatial Pseudotime")

print((p1 + p2) / (p3 + p4))
dev.off()

# 1.2 轨迹主干图（带分支点标注）
pdf(file = file.path(trajectory_dir, "Trajectory_Backbone.pdf"), width = 10, height = 8)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups = FALSE,
           label_leaves = TRUE,           # 标注叶节点
           label_branch_points = TRUE,    # 标注分支点
           label_roots = TRUE,           # 标注根节点
           graph_label_size = 3) +
  ggtitle("Trajectory Backbone with Branch Points") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
dev.off()

## ========== 轨迹比较与验证 ==========

# 7.1 不同方法的轨迹比较（如果有多种方法）
pdf(file = file.path(trajectory_dir, "Trajectory_Method_Comparison.pdf"), width = 14, height = 6)
p1 <- DimPlot(sp_obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("Seurat Clustering")
p2 <- DimPlot(sp_obj, reduction = "umap", group.by = "monocle3_clusters", label = TRUE) +
  ggtitle("Monocle3 Clustering")
print(p1 + p2)
dev.off()

# 7.2 轨迹一致性验证
pdf(file = file.path(trajectory_dir, "Trajectory_Consistency.pdf"), width = 10, height = 8)
consistency_data <- data.frame(
  pseudotime = sp_obj$pseudotime,
  seurat_cluster = sp_obj$seurat_clusters,
  monocle_cluster = sp_obj$monocle3_clusters
)

ggplot(consistency_data, aes(x = seurat_cluster, y = pseudotime)) +
  geom_violin(aes(fill = seurat_cluster), alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white") +
  labs(title = "Pseudotime Consistency Across Clustering Methods",
       x = "Seurat Clusters", y = "Pseudotime") +
  theme_minimal() +
  theme(legend.position = "none")
dev.off()
#
