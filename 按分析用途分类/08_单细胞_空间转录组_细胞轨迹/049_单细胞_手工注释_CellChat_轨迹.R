# ==========================================================================
# 脚本名     : 单细胞手工注释Marker验证_CellChat_轨迹.R
# 分类       : 08_单细胞_空间转录组_细胞轨迹
# 项目来源   : 从压缩包 199，单细胞手工注释.rar 整理
# 原始文件   : 199，单细胞手工注释\单细胞细胞注释和细胞通讯umap聚类 .R
# 用途       : 单细胞 Seurat 流程，重点加入 marker 文件辅助的手工注释、聚类-细胞类型富集评分，并继续输出细胞组成、目标基因和轨迹分析。
# 结果图     : QC小提琴图；QC散点图；高变基因图；PCA图/热图；UMAP聚类图；marker热图；Marker DotPlot；Marker小提琴图；细胞类型-聚类富集热图；手工注释UMAP；细胞组成柱状图；目标基因图；monocle3轨迹图
# 非肿瘤消化适配: 很适合。非肿瘤消化系统单细胞最常需要手工注释，这个模块可作为主力模板。
# 主要 R 包  : Seurat; CellChat; celldex; SingleR; monocle3; ggplot2; dplyr; tidyr; patchwork
# 整理日期   : 2026-05-13
# 备注       : 保留原始代码逻辑，仅添加统一说明头；运行前请把 workDir/setwd 和输入文件名改成当前项目路径。
# ==========================================================================
# --先运行548----------------- 预加载R包及依赖 --------------------
library(Seurat)       # 主要的单细胞分析包
library(dplyr)        # 数据整理、管道操作
library(ggplot2)      # 画图
library(magrittr)     # 管道符 %>%
library(RColorBrewer) # 颜色方案
library(limma)        # 差异分析
library(tidyr)        # 数据整理（长宽表互转等）
library(NMF)          # 非负矩阵分解相关（部分包依赖）
library(CellChat)     # 细胞通讯分析
library(ggalluvial)   # Alluvial 图
library(svglite)      # 导出svg格式图
library(celldex)      # SingleR注释参考数据库
library(SingleR)      # 自动细胞类型注释
#library(monocle)      # 单细胞轨迹分析
library(patchwork)
library(Matrix)
library(monocle3)
# ----------------- 设置工作目录 -----------------
workDir <- "H:\\常用分析生信\\199，单细胞手工注释"  # 工作目录
setwd(workDir)  # 切换到工作目录

# ----------- 1. 参数设置区 ------------
logFC_filter       <- 1        # 差异分析logFC阈值
p_adj_filter       <- 0.05     # 差异分析校正p值阈值
min_cells_gene     <- 5        # 至少在5个细胞出现的基因参与分析
min_genes_per_cell <- 100      # 每个细胞至少检测到的基因数
post_filter_cells  <- 300      # 二次过滤时，细胞内基因数大于此值
post_filter_mito   <- 20       # 线粒体比例不得高于此值（%）
n_top_var_features <- 2500     # 选取前2500个高度变异基因
n_pcs              <- 22        # PCA降维时使用的主成分数
n_topmarker_heat   <- 10       # 每簇热图展示top10 marker
cluster_resolution <- 0.6      # 聚类分辨率
neighbor_dims      <- 15       # 聚类与UMAP时用多少PC
qc_vln_width       <- 15       # 质控小提琴图宽度
qc_vln_height      <- 7        # 质控小提琴图高度
input_expr_file    <- "single_cell_data.csv"  # 输入表达数据文件名
output_dir         <- "analysis_results"      # 结果输出文件夹
# ----- 可自定义关注的基因（支持多基因） -----
target_genes <- c("CSF1R","CYP3A4","KDR","MAPK3")   # 可多个基因

# ----------- 2. 文件夹确认 -------
if(!dir.exists(output_dir)) dir.create(output_dir, recursive=TRUE)  # 如果输出目录不存在则创建
cat("输出目录设置为：", output_dir, "\n")  # 输出目录提示

# ----------- 3. 数据读取与初处理 -------
cat("Step 1: 读取表达数据...\n")  # 流程进度提示
if(!file.exists(input_expr_file)) stop("表达数据文件不存在！")  # 检查输入文件是否存在
expr_table <- read.csv(input_expr_file, header=TRUE, sep=",", stringsAsFactors=FALSE, check.names=FALSE) # 读取表达矩阵
if(nrow(expr_table)==0) stop("表达矩阵无内容！")  # 检查表格非空
rownames(expr_table) <- expr_table[,1]   # 首列为基因名，设为行名
gene_names <- expr_table[,1]             # 保存基因名向量
expr_values <- as.matrix(expr_table[,-1]) # 去掉第一列后的原始表达值
expr_matrix <- expr_values                # 表达矩阵后续用于Seurat

# ----------- 4. 创建Seurat对象及QC小提琴图 ---------
cat("Step 2: 创建Seurat对象...\n") # 流程提示
scObject <- CreateSeuratObject(
  counts = as.matrix(expr_matrix),         # 原始表达量
  project = "scRNAseqProject",             # 项目名
  min.cells = min_cells_gene,              # 参与分析的基因必须在多少细胞中存在
  min.features = min_genes_per_cell,       # 细胞需检测多少基因
  names.delim = "_"                        # 分隔符
)
scObject[["percent.mito"]] <- PercentageFeatureSet(scObject, pattern = "^MT-") # 线粒体基因比例
pdf(file.path(output_dir, "QC_violin_basicMetrics_pubstyle.pdf"), width=12, height=6)
p <- VlnPlot(
  scObject,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mito"),
  pt.size = 0,
  group.by = "orig.ident"
) +
  theme_classic(base_size = 16) +  # 高级美化
  theme(
    legend.position = "right",
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.text.x  = element_text(size = 12, face = "bold", angle = 45, vjust=1, hjust=1),
    axis.text.y  = element_text(size = 12, face = "bold"),
    strip.text   = element_text(size = 15, face = "bold", color = "black")
  )
print(p)  # 必须print
dev.off()

# ----------- 5. 二次过滤与QC散点 -------
cat("Step 3: 二次细胞质控...\n") # 流程提示
cell_num_before <- ncol(scObject) # 过滤前细胞数
scObject <- subset(scObject, subset = nFeature_RNA > post_filter_cells & percent.mito < post_filter_mito) # 基于基因数和线粒体比过滤
cell_num_after <- ncol(scObject) # 过滤后细胞数
cat(sprintf("细胞过滤: 原%d，剩%d。\n", cell_num_before, cell_num_after)) # 输出细胞过滤信息
if(cell_num_after < 10) stop("过滤后过少细胞，检查阈值！") # 过少细胞则终止

pdf(file.path(output_dir, "QC_scatter_metrics.pdf"), width=13, height=7) # PDF输出
plot1 <- FeatureScatter(scObject, feature1 = "nCount_RNA", feature2 = "percent.mito", pt.size = 1.5)   # UMI-线粒体比例
plot2 <- FeatureScatter(scObject, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 1.5)   # UMI-基因数
CombinePlots(plots = list(plot1, plot2)) # 合并图
dev.off()

# ----------- 6. 归一化及变异基因 -------
cat("Step 4: 归一化+高变基因...\n") # 流程提示
scObject <- NormalizeData(scObject, normalization.method="LogNormalize", scale.factor=10000) # 归一化
scObject <- FindVariableFeatures(scObject, selection.method="vst", nfeatures=n_top_var_features) # 高变基因
pdf(file.path(output_dir, "VarGenes_overview.pdf"), width=10, height=6) # 输出高变基因分布
VariableFeaturePlot(scObject)
dev.off()

# ----------- 7. PCA降维与可视化 ----------
cat("Step 5: PCA降维...\n") # 流程提示
scObject <- ScaleData(scObject) # 标准化数据（均值0方差1）
scObject <- RunPCA(scObject, npcs=n_pcs, features=VariableFeatures(scObject)) # PCA降维，n_pcs主成分
pdf(file.path(output_dir, "PCA_geneLoadings.pdf")) # 主成分主要基因可视化
VizDimLoadings(scObject, dims=1:3, reduction="pca", nfeatures=25)
dev.off()
pdf(file.path(output_dir, "PCA_heatmap_topGenes.pdf"), width=8, height=7) # 前三主成分主要基因热图
DimHeatmap(scObject, dims=1:3, cells=400, balanced=TRUE)
dev.off()




# ----------- 8. 聚类与UMAP降维 ----------
cat("Step 6: 聚类与UMAP降维...\n") # 流程提示
scObject <- FindNeighbors(scObject, dims=1:neighbor_dims)      # 构建K近邻图
scObject <- FindClusters(scObject, resolution=cluster_resolution) # 聚类
scObject <- RunUMAP(scObject, dims=1:neighbor_dims)            # UMAP降维
nColors <- length(unique(scObject$seurat_clusters))             # 统计聚类类别数量
myColors <- colorRampPalette(brewer.pal(12, "Set3"))(nColors)  # 自定义颜色
pdf(file.path(output_dir, "UMAP_clustered_samples.pdf"), width=7, height=5)
DimPlot(scObject, reduction="umap", label=TRUE, cols=myColors) # 绘制UMAP并标注聚类
dev.off()
write.csv(data.frame(Cell=colnames(scObject), Cluster=scObject$seurat_clusters), 
          file=file.path(output_dir, "CellCluster_UMAP_assignments.csv"), row.names=FALSE) # 保存每个细胞的聚类信息

# ----------- 9. 差异基因查找/热图 -----------
cat("Step 7: 差异表达分析...\n") # 流程提示
marker_df <- tryCatch({
  FindAllMarkers(scObject, only.pos=TRUE, min.pct=0.2, logfc.threshold=logFC_filter) # 所有聚类差异表达基因
}, error=function(e){cat("Marker分析异常\n"); NULL})
if(!is.null(marker_df) && nrow(marker_df)>0){
  sigMarkers <- marker_df[abs(marker_df$avg_log2FC)>logFC_filter & marker_df$p_val_adj<p_adj_filter,]
  write.csv(sigMarkers, file=file.path(output_dir, "Cluster_Markers_DEGs.csv"), row.names=FALSE)
  topmarker <- marker_df %>% group_by(cluster) %>% top_n(n_topmarker_heat, avg_log2FC)
  # 检查基因是否都在对象中
  valid_genes <- topmarker$gene[topmarker$gene %in% rownames(scObject)]
  if(length(valid_genes) > 0){
    pdf(file.path(output_dir, "Markers_DoHeatmap.pdf"))
    print(DoHeatmap(scObject, features=valid_genes, size=4) + NoLegend())
    dev.off()
  } else {
    cat("未找到可用于热图的marker基因。\n")
  }
} else {
  cat('未检出显著marker。\n')
}

###########
# ----------- 7.5. 细胞类型标记基因与聚类相关性分析（在手工注释前） -----------
cat("Step 7.5: 细胞类型标记基因与聚类相关性分析...\n")

# 首先检查聚类是否存在
if(!"seurat_clusters" %in% colnames(scObject@meta.data)) {
  cat("警告：未找到聚类结果，请先完成聚类分析\n")
} else {
  
  # 1. 读取细胞类型标记基因注释文件
  marker_file <- "Cell Annotation.txt"
  if(!file.exists(marker_file)) {
    cat("警告：标记基因文件", marker_file, "不存在，跳过标记基因-聚类相关性分析\n")
  } else {
    cat("读取标记基因文件：", marker_file, "\n")
    
    # 读取标记基因文件，跳过第一行（标题行）
    refMarker <- read.table(marker_file, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
    
    # 解析基因列表
    genes <- list()
    for(i in 1:nrow(refMarker)){
      celltype_name <- refMarker[i, 1]
      # 跳过标题行
      if(celltype_name == "Cell Type (English)" || celltype_name == "Cell Type") next
      
      gene_string <- refMarker[i, 2]
      if(!is.na(gene_string) && gene_string != "Marker Genes") {
        genes[[celltype_name]] <- unlist(strsplit(gene_string, "\\,|\\;|\\|"))  # 支持多种分隔符
        genes[[celltype_name]] <- trimws(genes[[celltype_name]])  # 去除空格
        genes[[celltype_name]] <- genes[[celltype_name]][genes[[celltype_name]] != ""]  # 去除空字符
      }
    }
    
    # 移除空的基因列表
    genes <- genes[sapply(genes, length) > 0]
    
    cat("解析的细胞类型和基因：\n")
    for(ct in names(genes)) {
      cat(sprintf("  %s: %d个基因 (%s)\n", ct, length(genes[[ct]]), 
                  paste(head(genes[[ct]], 3), collapse = ", ")))
    }
    
    if(length(genes) == 0) {
      cat("警告：没有找到有效的基因列表，请检查文件格式\n")
    } else {
      
      # 创建标记基因分析结果文件夹
      marker_analysis_dir <- file.path(output_dir, "Marker_Gene_Cluster_Analysis")
      if(!dir.exists(marker_analysis_dir)) dir.create(marker_analysis_dir, recursive = TRUE)
      
      # 2. 验证基因是否在数据中存在
      all_marker_genes <- unique(unlist(genes))
      valid_marker_genes <- all_marker_genes[all_marker_genes %in% rownames(scObject)]
      missing_genes <- setdiff(all_marker_genes, valid_marker_genes)
      
      cat(sprintf("总标记基因: %d, 在数据中找到: %d, 缺失: %d\n", 
                  length(all_marker_genes), length(valid_marker_genes), length(missing_genes)))
      
      if(length(missing_genes) > 0) {
        write.csv(data.frame(Missing_Genes = missing_genes), 
                  file = file.path(marker_analysis_dir, "Missing_Marker_Genes.csv"), 
                  row.names = FALSE)
        cat("缺失的基因：", paste(missing_genes, collapse = ", "), "\n")
      }
      
      # 过滤基因列表，只保留存在的基因
      genes_filtered <- lapply(genes, function(x) x[x %in% rownames(scObject)])
      genes_filtered <- genes_filtered[sapply(genes_filtered, length) > 0]  # 移除空列表
      
      if(length(genes_filtered) == 0) {
        cat("警告：没有找到有效的标记基因，跳过分析\n")
      } else {
        
        # 获取聚类信息
        cluster_levels <- sort(unique(scObject$seurat_clusters))
        
        # 3. 绘制标记基因的气泡图
        cat("绘制标记基因气泡图...\n")
        
        # 检查是否有足够的基因用于分析
        total_genes_for_plot <- sum(sapply(genes_filtered, length))
        if(total_genes_for_plot > 0) {
          
          # 使用Seurat的DotPlot函数
          tryCatch({
            # 如果基因过多，限制每个细胞类型最多显示前10个基因
            genes_for_plot <- lapply(genes_filtered, function(x) head(x, 10))
            
            # 创建带有细胞类型标签的基因名列表
            labeled_genes_for_plot <- list()
            for(celltype in names(genes_for_plot)) {
              # 简化细胞类型名称（取前15个字符，避免标签过长）
              celltype_short <- substr(celltype, 1, 15)
              if(nchar(celltype) > 15) celltype_short <- paste0(celltype_short, "...")
              
              # 为每个基因添加细胞类型前缀
              labeled_genes <- paste0(celltype_short, "_", genes_for_plot[[celltype]])
              names(labeled_genes) <- genes_for_plot[[celltype]]  # 保持原基因名用于数据提取
              labeled_genes_for_plot[[celltype]] <- labeled_genes
            }
            
            # 为了绘图，我们需要创建一个基因映射表
            gene_mapping <- data.frame()
            plot_gene_list <- list()
            
            for(celltype in names(genes_for_plot)) {
              celltype_short <- substr(celltype, 1, 15)
              if(nchar(celltype) > 15) celltype_short <- paste0(celltype_short, "...")
              
              original_genes <- genes_for_plot[[celltype]]
              labeled_genes <- paste0(celltype_short, "_", original_genes)
              
              # 创建映射关系
              temp_mapping <- data.frame(
                original_gene = original_genes,
                labeled_gene = labeled_genes,
                celltype = celltype,
                stringsAsFactors = FALSE
              )
              gene_mapping <- rbind(gene_mapping, temp_mapping)
              
              # 为绘图准备基因列表（使用原基因名）
              plot_gene_list[[celltype]] <- original_genes
            }
            
      #      # 绘制基础的DotPlot
          #  pdf(file.path(marker_analysis_dir, "Marker_Gene_DotPlot.pdf"), width = 12, height = 16)
            
            p_dot <- DotPlot(scObject, 
                             features = plot_gene_list, 
                             group.by = "seurat_clusters",
                             dot.scale = 8,
                             cols = c("lightgrey", "red")) +
              theme_classic(base_size = 12) +
              theme(
                axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", size = 10),
                axis.text.y = element_text(face = "bold", size = 10),
                plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
                legend.title = element_text(face = "bold"),
                strip.text = element_text(face = "bold", size = 10)
              ) +
              labs(
                title = "Marker Genes Expression Across Clusters",
                x = "Clusters", 
                y = "Genes"
              )
            
            # 修改Y轴标签，添加细胞类型信息
            plot_data <- p_dot$data
            
            # 创建新的基因标签
            new_gene_labels <- c()
            for(gene in levels(plot_data$features.plot)) {
              matching_row <- gene_mapping[gene_mapping$original_gene == gene, ]
              if(nrow(matching_row) > 0) {
                celltype_short <- substr(matching_row$celltype[1], 1, 15)
                if(nchar(matching_row$celltype[1]) > 15) celltype_short <- paste0(celltype_short, "...")
                new_label <- paste0(celltype_short, "_", gene)
              } else {
                new_label <- gene
              }
              new_gene_labels <- c(new_gene_labels, new_label)
            }
            
            # 应用新标签
            p_dot <- p_dot + 
              scale_y_discrete(labels = new_gene_labels) +
              scale_x_discrete(labels = paste("Cluster", levels(plot_data$id)))
            
            print(p_dot)
            dev.off()
          
            
            
            pdf(file.path(marker_analysis_dir, "Marker_Gene_DotPlot_Flipped.pdf"), width = 30, height = 12)
            
            p_dot_flipped <- DotPlot(scObject, 
                                     features = plot_gene_list, 
                                     group.by = "seurat_clusters",
                                     dot.scale = 8,
                                     cols = c("lightgrey", "red")) +
              theme_classic(base_size = 12) +
              theme(
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "bold", size = 9),
                axis.text.y = element_text(face = "bold", size = 10),
                plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
                legend.title = element_text(face = "bold"),
                strip.text = element_text(face = "bold", size = 10)
              ) +
              labs(
                title = "Marker Genes Expression Across Clusters (Flipped)",
                x = "Genes", 
                y = "Clusters"
              ) +
              coord_flip()  # 这里进行坐标翻转
            
            # 对翻转后的图应用标签
            p_dot_flipped <- p_dot_flipped + 
              scale_x_discrete(labels = new_gene_labels) +
              scale_y_discrete(labels = paste("Cluster", levels(plot_data$id)))
            
            print(p_dot_flipped)
            dev.off()
            cat("已输出: Marker_Gene_DotPlot_Flipped.pdf\n")
            
          }, error = function(e) {
            cat("绘制气泡图时出错：", e$message, "\n")
            
            # 尝试使用VlnPlot作为替代
            tryCatch({
              pdf(file.path(marker_analysis_dir, "Marker_Gene_ViolinPlot.pdf"), width = 16, height = 12)
              
              # 选择前20个基因进行展示
              plot_genes <- head(unlist(genes_for_plot), 20)
              
              p_vln <- VlnPlot(scObject, 
                               features = plot_genes, 
                               group.by = "seurat_clusters",
                               ncol = 4,
                               pt.size = 0) +
                theme(plot.title = element_text(hjust = 0.5, face = "bold"))
              
              print(p_vln)
              dev.off()
              cat("已输出备选图: Marker_Gene_ViolinPlot.pdf\n")
            }, error = function(e2) {
              cat("备选绘图也失败：", e2$message, "\n")
            })
          })
          
          # 4. 计算每个细胞类型标记基因的聚类富集分数
          cat("计算细胞类型-聚类富集分数...\n")
          
          enrichment_results <- data.frame()
          
          for(celltype in names(genes_filtered)) {
            marker_genes <- genes_filtered[[celltype]]
            cat(sprintf("  分析细胞类型: %s (%d个基因)\n", celltype, length(marker_genes)))
            
            # 计算每个聚类中这些标记基因的平均表达
            cluster_scores <- data.frame()
            
            for(cluster in cluster_levels) {
              cluster_cells <- colnames(scObject)[scObject$seurat_clusters == cluster]
              
              if(length(cluster_cells) == 0) next
              
              # 计算该聚类中标记基因的表达
              tryCatch({
                if(length(marker_genes) == 1) {
                  gene_expr <- FetchData(scObject, vars = marker_genes, cells = cluster_cells)[,1]
                } else {
                  gene_expr_matrix <- FetchData(scObject, vars = marker_genes, cells = cluster_cells)
                  gene_expr <- rowMeans(gene_expr_matrix, na.rm = TRUE)  # 取平均表达
                }
                
                cluster_score <- data.frame(
                  CellType = celltype,
                  Cluster = as.character(cluster),
                  Mean_Expression = mean(gene_expr, na.rm = TRUE),
                  Median_Expression = median(gene_expr, na.rm = TRUE),
                  Positive_Cells = sum(gene_expr > 0, na.rm = TRUE),
                  Total_Cells = length(gene_expr),
                  Positive_Percentage = round(sum(gene_expr > 0, na.rm = TRUE) / length(gene_expr) * 100, 2),
                  Gene_Count = length(marker_genes),
                  stringsAsFactors = FALSE
                )
                
                cluster_scores <- rbind(cluster_scores, cluster_score)
              }, error = function(e) {
                cat(sprintf("    聚类 %s 分析出错: %s\n", cluster, e$message))
              })
            }
            
            enrichment_results <- rbind(enrichment_results, cluster_scores)
          }
          
          if(nrow(enrichment_results) > 0) {
            # 输出富集结果
            write.csv(enrichment_results, 
                      file = file.path(marker_analysis_dir, "CellType_Cluster_Enrichment_Scores.csv"), 
                      row.names = FALSE)
            
            # 5. 生成细胞类型-聚类富集热图
            cat("生成细胞类型-聚类富集热图...\n")
            
            tryCatch({
              # 创建热图矩阵（平均表达）
              heatmap_matrix_mean <- enrichment_results %>%
                select(CellType, Cluster, Mean_Expression) %>%
                pivot_wider(names_from = Cluster, values_from = Mean_Expression, values_fill = 0) %>%
                column_to_rownames("CellType") %>%
                as.matrix()
              
              # 确保聚类按顺序排列
              cluster_order <- sort(as.numeric(colnames(heatmap_matrix_mean)))
              heatmap_matrix_mean <- heatmap_matrix_mean[, as.character(cluster_order), drop = FALSE]
              
              pdf(file.path(marker_analysis_dir, "CellType_Cluster_Enrichment_Heatmap_MeanExpression.pdf"), 
                  width = max(8, ncol(heatmap_matrix_mean) * 0.8), 
                  height = max(6, nrow(heatmap_matrix_mean) * 0.8))
              pheatmap::pheatmap(
                heatmap_matrix_mean,
                cluster_rows = TRUE,
                cluster_cols = FALSE,
                scale = "row",
                color = colorRampPalette(c("blue", "white", "red"))(100),
                main = "Cell Type Marker Genes Mean Expression Across Clusters",
                fontsize = 12,
                angle_col = 0
              )
              dev.off()
              
              # 创建热图矩阵（阳性比例）
              heatmap_matrix_pct <- enrichment_results %>%
                select(CellType, Cluster, Positive_Percentage) %>%
                pivot_wider(names_from = Cluster, values_from = Positive_Percentage, values_fill = 0) %>%
                column_to_rownames("CellType") %>%
                as.matrix()
              
              heatmap_matrix_pct <- heatmap_matrix_pct[, as.character(cluster_order), drop = FALSE]
              
              pdf(file.path(marker_analysis_dir, "CellType_Cluster_Enrichment_Heatmap_PositivePercentage.pdf"), 
                  width = max(8, ncol(heatmap_matrix_pct) * 0.8), 
                  height = max(6, nrow(heatmap_matrix_pct) * 0.8))
              pheatmap::pheatmap(
                heatmap_matrix_pct,
                cluster_rows = TRUE,
                cluster_cols = FALSE,
                scale = "none",
                color = colorRampPalette(c("white", "orange", "red"))(100),
                main = "Cell Type Marker Genes Positive Percentage Across Clusters",
                fontsize = 12,
                angle_col = 0
              )
              dev.off()
            }, error = function(e) {
              cat("生成热图时出错：", e$message, "\n")
            })
            
            # 7. 生成推荐的聚类-细胞类型对应表
            cat("生成聚类-细胞类型推荐对应表...\n")
            
            tryCatch({
              # 为每个聚类找到最高得分的细胞类型
              cluster_celltype_prediction <- enrichment_results %>%
                group_by(Cluster) %>%
                slice_max(Mean_Expression, n = 1, with_ties = FALSE) %>%
                select(Cluster, Predicted_CellType = CellType, Max_Score = Mean_Expression, 
                       Positive_Percentage, Gene_Count) %>%
                ungroup() %>%
                arrange(as.numeric(Cluster))
              
              write.csv(cluster_celltype_prediction, 
                        file = file.path(marker_analysis_dir, "Cluster_CellType_Prediction.csv"), 
                        row.names = FALSE)
              
              # 8. 输出分析总结
              cat("\n=== 细胞类型标记基因-聚类分析完成 ===\n")
              cat("输出文件夹：", marker_analysis_dir, "\n")
              cat("主要输出文件：\n")
        
              cat("- Marker_Gene_DotPlot_Flipped.pdf: 标记基因气泡图（坐标翻转版）\n")
              cat("- CellType_Cluster_Enrichment_Scores.csv: 细胞类型-聚类富集分数\n")
              cat("- CellType_Cluster_Enrichment_Heatmap_*.pdf: 富集热图\n")
              cat("- Cluster_CellType_Prediction.csv: 聚类-细胞类型预测结果\n\n")
              
              # 输出推荐的聚类注释
              cat("=== 基于标记基因的聚类注释推荐 ===\n")
              for(i in 1:nrow(cluster_celltype_prediction)) {
                cat(sprintf("聚类 %s: %s (得分: %.3f, 阳性比例: %.1f%%)\n",
                            cluster_celltype_prediction$Cluster[i],
                            cluster_celltype_prediction$Predicted_CellType[i],
                            cluster_celltype_prediction$Max_Score[i],
                            cluster_celltype_prediction$Positive_Percentage[i]))
              }
            }, error = function(e) {
              cat("生成预测结果时出错：", e$message, "\n")
            })
          } else {
            cat("警告：没有生成有效的富集结果\n")
          }
        }
      }
    }
  }
}
#
# 清理内存
gc()

# 先运行到这里-----------# -----------# -----------# -----------# -----------# -----------# -----------# -----------# -----------# -----------# -----------# -----------# -----------# -----------# -----------# -----------# -----------# -----------# -----------# -----------# -----------



#再运行后面， -----# -----------# -----------# -----------# -----------# -----------# -----------# -----------# -----------# -----------# -----------# -----------# -----------# -----------# -----------
cat("Step 手动细胞类型注释...\n")

# 首先检查当前聚类信息
cat("=== 调试信息 ===\n")
cat("聚类数据类型：", class(scObject$seurat_clusters), "\n")
cat("聚类唯一值：", paste(sort(unique(scObject$seurat_clusters)), collapse = ", "), "\n")
cat("聚类唯一值类型：", class(unique(scObject$seurat_clusters)), "\n")

# 您想要的注释
cluster_names <- c(
  "0" = "T Cells",                    
  "1" = "Monocytes/Macrophages",     
  "2" = "B Cells",                    
  "3" = "Monocytes/Macrophages",     
  "4" = "Monocytes/Macrophages",      
  "5" = "B Cells",                    
  "6" = "Monocytes/Macrophages",     
  "7" = "Monocytes/Macrophages",     
  "8" = "Neurons",                  
  "9" = "Monocytes/Macrophages",     
  "10" = "Monocytes/Macrophages",    
  "11" = "Neurons",                  
  "12" = "Monocytes/Macrophages",    
  "13" = "Monocytes/Macrophages",     
  "14" = "Monocytes/Macrophages",    
  "15" = "B Cells",                   
  "16" = "Neurons",                   
  "17" = "Monocytes/Macrophages"      
)

# 方法1：最直接的方法 - 创建映射向量
# 将聚类转换为字符串进行匹配
cell_clusters <- as.character(scObject$seurat_clusters)
cat("转换后的聚类示例：", head(cell_clusters), "\n")

# 创建细胞类型向量
Manual_celltype <- character(length(cell_clusters))
names(Manual_celltype) <- colnames(scObject)

# 逐个匹配
for(i in 1:length(cell_clusters)) {
  cluster_id <- cell_clusters[i]
  if(cluster_id %in% names(cluster_names)) {
    Manual_celltype[i] <- cluster_names[cluster_id]
  } else {
    Manual_celltype[i] <- "Unknown"
  }
}

# 检查映射结果
cat("映射前10个细胞的结果：\n")
for(i in 1:min(10, length(cell_clusters))) {
  cat("细胞", i, "聚类", cell_clusters[i], "->", Manual_celltype[i], "\n")
}

# 统计映射结果
cat("手动注释统计（添加前）：\n")
print(table(Manual_celltype))

# 方法2：如果方法1不工作，尝试直接使用数字索引
if(all(Manual_celltype == "Unknown")) {
  cat("方法1失败，尝试方法2...\n")
  
  # 将聚类转换为数字再转回字符串
  numeric_clusters <- as.numeric(as.character(scObject$seurat_clusters))
  string_clusters <- as.character(numeric_clusters)
  
  Manual_celltype <- cluster_names[string_clusters]
  names(Manual_celltype) <- colnames(scObject)
  
  # 处理NA值
  Manual_celltype[is.na(Manual_celltype)] <- "Unknown"
  
  cat("方法2映射结果：\n")
  print(table(Manual_celltype))
}

# 方法3：如果还是不行，使用switch语句
if(all(Manual_celltype == "Unknown")) {
  cat("方法2失败，尝试方法3...\n")
  
  Manual_celltype <- sapply(as.character(scObject$seurat_clusters), function(x) {
    switch(x,
           "0" = "T Cells",
           "1" = "Monocytes/Macrophages",
           "2" = "B Cells",
           "3" = "Monocytes/Macrophages",
           "4" = "Monocytes/Macrophages",
           "5" = "B Cells",
           "6" = "Monocytes/Macrophages",
           "7" = "Monocytes/Macrophages",
           "8" = "Neurons",
           "9" = "Monocytes/Macrophages",
           "10" = "Monocytes/Macrophages",
           "11" = "Neurons",
           "12" = "Monocytes/Macrophages",
           "13" = "Monocytes/Macrophages",
           "14" = "Monocytes/Macrophages",
           "15" = "B Cells",
           "16" = "Neurons",
           "17" = "Monocytes/Macrophages",
           "Unknown")  # 默认值
  })
  
  names(Manual_celltype) <- colnames(scObject)
  
  cat("方法3映射结果：\n")
  print(table(Manual_celltype))
}

# 确保没有NA值
if(any(is.na(Manual_celltype))) {
  cat("发现NA值，替换为Unknown\n")
  Manual_celltype[is.na(Manual_celltype)] <- "Unknown"
}

# 添加到Seurat对象
scObject$Manual_celltype <- Manual_celltype

# 验证结果
cat("=== 最终验证 ===\n")
cat("细胞类型注释统计：\n")
final_table <- table(scObject$Manual_celltype)
print(final_table)

# 显示每个聚类对应的细胞类型
cat("聚类-细胞类型对应关系验证：\n")
cluster_celltype_check <- table(scObject$seurat_clusters, scObject$Manual_celltype)
print(cluster_celltype_check)

# 检查是否成功
if(length(unique(scObject$Manual_celltype)) > 1 && !all(scObject$Manual_celltype == "Unknown")) {
  cat("✓ 注释成功！\n")
} else {
  cat("✗ 注释失败，所有细胞仍为Unknown\n")
  
  # 进一步调试
  cat("进一步调试信息：\n")
  cat("scObject$seurat_clusters的前10个值：\n")
  print(head(scObject$seurat_clusters, 10))
  cat("cluster_names的名称：\n")
  print(names(cluster_names))
  cat("匹配测试：\n")
  test_cluster <- as.character(scObject$seurat_clusters[1])
  cat("第一个聚类值：'", test_cluster, "'\n", sep="")
  cat("是否在cluster_names中：", test_cluster %in% names(cluster_names), "\n")
}

# 创建manual_anno数据框
unique_clusters <- sort(unique(as.character(scObject$seurat_clusters)))
manual_anno <- data.frame(
  Cluster = unique_clusters,
  CellType = cluster_names[unique_clusters],
  stringsAsFactors = FALSE
)

# 处理NA值
manual_anno$CellType[is.na(manual_anno$CellType)] <- "Unknown"

# 输出文件
write.csv(manual_anno, file = file.path(output_dir, "CellCluster_ManualAnno.csv"), row.names=FALSE)

# 合并marker和手动注释
if(exists("marker_df") && !is.null(marker_df) && nrow(marker_df) > 0) {
  marker_df$cluster <- as.character(marker_df$cluster)
  marker_df_anno <- marker_df %>% left_join(manual_anno, by = c("cluster" = "Cluster"))
  write.csv(marker_df_anno, file=file.path(output_dir,"Cluster_Markers_DEGs_with_ManualCellType.csv"), row.names=FALSE)
}

# 绘制UMAP图
unique_celltypes <- unique(scObject$Manual_celltype)
nCellTypes <- length(unique_celltypes)
celltypeColors <- colorRampPalette(brewer.pal(min(12, max(3, nCellTypes)), "Set3"))(nCellTypes)

pdf(file.path(output_dir, "UMAP_celltype_manualAnnot_Set3.pdf"), width = 10, height = 8)
p <- DimPlot(scObject, group.by = "Manual_celltype", reduction = "umap", 
             label = TRUE, label.size = 4, cols = celltypeColors) +
  ggtitle("Manual Cell Type Annotation") +
  theme(plot.title = element_text(hjust = 0.5))
print(p)
dev.off()

cat("手动细胞类型注释完成！\n")


# ----------- 11. 细胞类型统计可视化及marker数 Alluvial等 -----------
cat("Step 9: 各种统计可视化输出...\n") # 流程提示
meta_df <- scObject@meta.data
celltype_stats <- as.data.frame(table(
  Sample = meta_df$orig.ident, 
  CellType = meta_df$Manual_celltype
))         # 按样本和细胞类型做数量统计
celltype_stats <- celltype_stats %>%
  group_by(Sample) %>%
  mutate(Ratio = Freq / sum(Freq)) %>%
  ungroup()  # 计算各样本细胞类型的比例
celltype_stats$Sample <- factor(celltype_stats$Sample, levels = unique(celltype_stats$Sample)) # 因子化
celltype_stats$CellType <- factor(celltype_stats$CellType, levels = unique(celltype_stats$CellType)) # 因子化
cell_types <- levels(celltype_stats$CellType) # 获得全部细胞类型
nColors <- length(cell_types) # 类型数
myColors <- colorRampPalette(brewer.pal(12, "Set3"))(nColors) # 颜色
names(myColors) <- cell_types # 命名

pdf(file.path(output_dir, "Sample_CellType_Composition_barplot.pdf"), width=7, height=6) # 打开pdf
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
  ) # 绘制每个样本各细胞类型比例的条形图
dev.off()

cell_counts <- as.data.frame(table(Sample = meta_df$orig.ident, CellType = meta_df$Manual_celltype)) # 统计数量
write.csv(cell_counts, file.path(output_dir, "CellNumber_perSample_CellType.csv"), row.names = FALSE) # 输出分组计数表
celltype_matrix <- celltype_stats %>% select(Sample, CellType, Ratio) %>%
  tidyr::spread(CellType, Ratio, fill = 0)
write.csv(celltype_matrix, file.path(output_dir, "CellType_Proportion_Matrix.csv"), row.names=FALSE) # 输出矩阵
cell_annot <- data.frame(CellBarcode = rownames(meta_df), Sample = meta_df$orig.ident,
                         Cluster = meta_df$seurat_clusters, CellType = meta_df$Manual_celltype)
write.csv(cell_annot, file.path(output_dir, "Cell_Full_Annotation.csv"), row.names=FALSE) # 输出所有细胞信息注释表
# marker计数及绘图
if (exists("marker_df_anno") && nrow(marker_df_anno) > 0) {
  gene_count_by_celltype <- marker_df_anno %>%
    group_by(CellType) %>%
    summarise(gene_number = n_distinct(gene)) %>%
    arrange(desc(gene_number))
  write.csv(gene_count_by_celltype, file = file.path(output_dir, "GeneNumber_per_CellType.csv"), row.names = FALSE)
  pdf(file.path(output_dir, "MarkerGeneCount_perCellType_bar.pdf"), width=8, height=5)
  p_bar <- ggplot(gene_count_by_celltype, aes(x = reorder(CellType, gene_number), y = gene_number, fill=CellType)) +
    geom_bar(stat="identity", width=0.7, color="grey40") +
    scale_fill_manual(values=colorRampPalette(brewer.pal(8,"Set2"))(nrow(gene_count_by_celltype))) +
    coord_flip() +
    labs(x="Cell Type", y="Number of Marker Genes", title="Number of Marker Genes per Cell Type") +
    theme_classic(base_size=16) +
    theme(legend.position = "none")
  print(p_bar)
  dev.off()
  
  # Alluvial流图 - 移到这里面
  top_genes_per_CT <- marker_df_anno %>%
    group_by(CellType) %>% top_n(5, abs(avg_log2FC)) %>% ungroup()
  pdf(file.path(output_dir, "MarkerGene_CellType_Alluvial.pdf"), width=9, height=6)
  p_alluvial <- ggplot(top_genes_per_CT, aes(axis1 = CellType, axis2 = gene, y = abs(avg_log2FC))) +
    geom_alluvium(aes(fill=CellType), width=1/12) +
    geom_stratum(width=1/6, fill="grey90", color="black") +
    geom_text(stat="stratum", aes(label=after_stat(stratum)), size=3) +
    scale_x_discrete(limits = c("CellType", "Gene"), expand = c(.05, .05)) +
    theme_minimal(base_size=14) +
    ggtitle("Top Marker Genes Alluvial Plot")
  print(p_alluvial)
  dev.off()
} else {
  cat("未找到marker注释数据，跳过marker相关图表\n")
}


# ----------- PCA方差解释与推荐主成分数输出 -----------
pca_stdev <- scObject[["pca"]]@stdev
pca_var_explained <- (pca_stdev^2) / sum(pca_stdev^2)
cumulative_var <- cumsum(pca_var_explained)
pca_table <- data.frame(
  PC = seq_along(pca_var_explained),
  Variance_Explained = pca_var_explained,
  Cumulative_Variance = cumulative_var
)
write.csv(pca_table, file = file.path(output_dir, "PCA_Variance_Explained.csv"), row.names=FALSE)

opt_cutoff <- 0.8  # 推荐主成分累计解释率阈值80%
best_n_pcs <- which(cumulative_var >= opt_cutoff)[1]
write.csv(
  data.frame(Recommended_n_PCs = best_n_pcs, Cumulative_Variance = cumulative_var[best_n_pcs]),
  file = file.path(output_dir, "PCA_Recommended_nPCs.csv"),
  row.names=FALSE
)

# 


for (target_gene in target_genes) {
  # ---- 检查基因是否在数据中 ----
  if(! target_gene %in% rownames(scObject)) {
    cat("警告：指定基因", target_gene, "不在表达矩阵中！\n")
    next
  }
  
  # --- 1. 绘制UMAP表达分布 ---
  pdf(file.path(output_dir, paste0(target_gene, "_UMAP_FeaturePlot.pdf")), width=7, height=6)
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
  
  # 按聚类
  pdf(file.path(output_dir, paste0(target_gene, "_VlnPlot_byCluster.pdf")), width=8, height=6)
  p2 <- VlnPlot(
    scObject, 
    features = target_gene, 
    group.by = "seurat_clusters", 
    pt.size = 0, # 不显示小黑点
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
  
  # --- 3. 导出分组均值（celltype/cluster） ---
  avg_by_celltype <- AverageExpression(scObject, features = target_gene, group.by = "Manual_celltype")$RNA
  avg_by_cluster  <- AverageExpression(scObject, features = target_gene, group.by = "seurat_clusters")$RNA
  
  write.csv(avg_by_celltype, file.path(output_dir, paste0(target_gene, "_Average_By_CellType.csv")))
  write.csv(avg_by_cluster,  file.path(output_dir, paste0(target_gene, "_Average_By_Cluster.csv")))
  
  # --- 4. 导出每个细胞表达量 ---
  gene_expr_vec <- FetchData(scObject, vars = target_gene)
  cell_anno <- data.frame(
    CellBarcode = colnames(scObject),
    CellType = scObject$Manual_celltype,
    Cluster = scObject$seurat_clusters,
    Gene_Expression = as.vector(gene_expr_vec[,1])
  )
  write.csv(cell_anno, file.path(output_dir, paste0(target_gene, "_Expression_perCell.csv")), row.names=FALSE)
  
  # ====== 自定义UMAP+条形图（聚类颜色一致，细胞类型全名显示） ======
  # 读取注释表
  anno_path <- file.path(output_dir, "CellCluster_ManualAnno.csv")
  anno_df <- read.csv(anno_path, stringsAsFactors = FALSE)
  cluster_col <- grep("cluster", colnames(anno_df), ignore.case = TRUE, value = TRUE)[1]
  celltype_col <- grep("cell.*type", colnames(anno_df), ignore.case = TRUE, value = TRUE)[1]
  if(is.na(cluster_col) | is.na(celltype_col)) stop("注释表中未找到cluster或celltype相关的列，请检查！")
  anno_df[[cluster_col]] <- as.character(anno_df[[cluster_col]])
  
  # UMAP数据
  umap_df <- as.data.frame(Embeddings(scObject, "umap"))
  umap_cols <- grep("UMAP|Dim", colnames(umap_df), value = TRUE, ignore.case = TRUE)
  if(length(umap_cols) >= 2) {
    colnames(umap_df)[match(umap_cols[1:2], colnames(umap_df))] <- c("UMAP_1", "UMAP_2")
  } else if(ncol(umap_df) == 2) {
    colnames(umap_df) <- c("UMAP_1", "UMAP_2")
  } else {
    stop("找不到UMAP坐标列，请检查umap_df的列名！")
  }
  umap_df$Cluster <- as.character(scObject$seurat_clusters)
  umap_df$CellBarcode <- rownames(umap_df)
  gene_found <- grep(paste0("^", target_gene, "$"), rownames(scObject), value = TRUE, ignore.case = TRUE)
  if(length(gene_found) == 0) stop(paste("指定基因", target_gene, "不在表达矩阵中！"))
  gene_to_use <- gene_found[1]
  umap_df$Expression <- FetchData(scObject, vars = gene_to_use)[,1]
  
  # 合并注释
  umap_df <- left_join(umap_df, anno_df, by = setNames(cluster_col, "Cluster"))
  if(!(celltype_col %in% colnames(umap_df))) stop(paste("合并后未找到细胞类型列", celltype_col, "！"))
  
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
  
  pdf(file.path(output_dir, paste0(target_gene, "_UMAP_withLegend.pdf")), width=12, height=6)
  print(p_final)
  dev.off()
  cat("已输出自定义UMAP表达图：", file.path(output_dir, paste0(target_gene, "_UMAP_withLegend.pdf")), "\n")
}

cat("所有结果/表格/图形均已输出到文件夹：", output_dir, "\n") # 流程结束提示

# 


# 1. 提取表达矩阵和meta信息
expr_matrix <- GetAssayData(scObject, assay = "RNA", slot = "counts")
cell_metadata <- scObject@meta.data
gene_metadata <- data.frame(gene_short_name = rownames(expr_matrix))
rownames(gene_metadata) <- rownames(expr_matrix)

# 2. 构建monocle3的cds对象
cds <- new_cell_data_set(expr_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_metadata)

# 3. 预处理和降维
cds <- preprocess_cds(cds, num_dim = 30)
cds <- reduce_dimension(cds, reduction_method = "UMAP")
cds <- cluster_cells(cds)
cds <- learn_graph(cds)

# 4. 轨迹推断（指定 reduction_method = "UMAP"）
cds <- order_cells(cds, reduction_method = "UMAP")

# 5. 绘图并保存为PDF

# 5.1 按pseudotime着色
pdf("Trajectory_Pseudotime.pdf", width=7, height=6)
print(plot_cells(cds, color_cells_by = "pseudotime"))
dev.off()

# 5.2 按细胞类型着色
if("Manual_celltype" %in% colnames(colData(cds))) {
  pdf("Trajectory_CellType.pdf", width=7, height=6)
  print(plot_cells(cds, color_cells_by = "Manual_celltype"))
  dev.off()
}
##

# 假设cds已经完成了 preprocess_cds, reduce_dimension, cluster_cells, learn_graph, order_cells
# target_genes <- c("EZH2", "CENPA", "KIF2C")

for (gene in target_genes) {
  # 检查基因是否在cds对象中
  if (!(gene %in% rownames(cds))) {
    cat("警告：基因", gene, "不在cds对象中，跳过。\n")
    next
  }
  
  # 1. 获取该基因在所有细胞的表达量
  expr_vec <- Matrix::t(assay(cds)[gene, ])  # 取该基因的表达量，转为列向量
  expr_vec <- as.numeric(expr_vec)
  names(expr_vec) <- colnames(cds)
  
  # 2. 按中位数分组
  med <- median(expr_vec)
  expr_group <- ifelse(expr_vec > med, "High", "Low")
  
  # 3. 写入cds的colData，临时新列
  colData(cds)[[paste0(gene, "_Group")]] <- factor(expr_group, levels = c("Low", "High"))
  
  # 4. 绘图
  pdf(paste0("Trajectory_", gene, "_HighLow.pdf"), width=7, height=6)
  print(
    plot_cells(
      cds,
      color_cells_by = paste0(gene, "_Group"),
      show_trajectory_graph = TRUE,
      label_groups_by_cluster = FALSE,
      label_leaves = TRUE,
      label_branch_points = TRUE,
      graph_label_size = 4,
      cell_size = 1.5,
      reduction_method = "UMAP"
    ) +
      scale_color_manual(values = c("Low" = "#079EDF", "High" = "#D377A9")) +
      ggtitle(paste("Trajectory -", gene, "(High/Low by median)"))
  )
  dev.off()
  cat("已输出：Trajectory_", gene, "_HighLow.pdf\n")
}





# 获取伪时序向量
pseudotime <- pseudotime(cds)
pseudotime <- as.numeric(pseudotime)
names(pseudotime) <- colnames(cds)

for (gene in target_genes) {
  # 检查基因是否在cds对象中
  if (!(gene %in% rownames(cds))) {
    cat("警告：基因", gene, "不在cds对象中，跳过。\n")
    next
  }
  
  # 获取该基因在所有细胞的表达量
  expr_vec <- Matrix::t(assay(cds)[gene, ])
  expr_vec <- as.numeric(expr_vec)
  names(expr_vec) <- colnames(cds)
  
  # 按中位数分组
  med <- median(expr_vec)
  expr_group <- ifelse(expr_vec > med, "High", "Low")
  
  # 构建数据框
  plot_df <- data.frame(
    Pseudotime = pseudotime,
    Group = factor(expr_group, levels = c("Low", "High"))
  )
  plot_df <- plot_df[!is.na(plot_df$Pseudotime), ]  # 去除NA
  
  # 绘制SCI风格密度图
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
      axis.line = element_line(size = 1.2),
      axis.ticks = element_line(size = 1.2),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(face = "bold", size = 16),
      legend.key = element_blank(),
      panel.grid = element_blank()
    )
  
  # 输出PDF
  pdf(paste0("Pseudotime_Density_", gene, "_HighLow.pdf"), width=6, height=5)
  print(p)
  dev.off()
  cat("已输出：Pseudotime_Density_", gene, "_HighLow.pdf\n")
}



# 假设你已经有scObject和output_dir
# target_genes <- c("EZH2", "CENPA", "KIF2C")
# output_dir <- "analysis_results"

# 1. 提取UMAP坐标和注释
umap_df <- as.data.frame(Embeddings(scObject, "umap"))
colnames(umap_df)[1:2] <- c("UMAP_1", "UMAP_2")
umap_df$Cell <- rownames(umap_df)
umap_df$CellType <- scObject$Manual_celltype[umap_df$Cell]

for (gene in target_genes) {
  cat("当前基因：", gene, "\n")
  if (!(gene %in% rownames(scObject))) {
    cat("警告：基因", gene, "不在scObject中，跳过。\n")
    next
  }
  # 获取表达量并分组
  expr_vec <- FetchData(scObject, vars = gene)[,1]
  med <- median(expr_vec)
  expr_group <- ifelse(expr_vec > med, "High", "Low")
  
  plot_df <- umap_df
  plot_df$Group <- factor(expr_group, levels = c("Low", "High"))
  plot_df$Expression <- expr_vec
  
  # 创建基因专属文件夹
  gene_dir <- file.path(output_dir, paste0(gene, "_Celltype_UMAP"))
  if (!dir.exists(gene_dir)) dir.create(gene_dir, recursive = TRUE)
  
  for (ct in sort(unique(na.omit(plot_df$CellType)))) {
    sub_df <- subset(plot_df, CellType == ct)
    cat("  注释类型：", ct, "，细胞数：", nrow(sub_df), "\n")
    if (nrow(sub_df) == 0) next
    
    # 颜色与图例标签（高=红，低=蓝）
    sub_df$GroupLabel <- ifelse(sub_df$Group == "High",
                                paste0("High ", gene, " ", ct),
                                paste0("Low ", gene, " ", ct))
    color_map <- setNames(
      c("#4682B4", "#CD2626"),
      c(paste0("Low ", gene, " ", ct), paste0("High ", gene, " ", ct))
    )
    
    p <- ggplot(sub_df, aes(x = UMAP_1, y = UMAP_2, color = GroupLabel)) +
      geom_point(size = 1.2, alpha = 0.8) +
      scale_color_manual(values = color_map) +
      labs(
        title = paste0(gene, " expression in ", ct, " (by median)"),
        color = NULL
      ) +
      theme_classic(base_size = 18) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        axis.title = element_blank(),
        axis.text = element_text(face = "bold", size = 14),
        axis.line = element_line(size = 1.1),
        axis.ticks = element_line(size = 1.1),
        legend.position = "top",
        legend.text = element_text(face = "bold", size = 14),
        legend.key = element_blank(),
        panel.grid = element_blank()
      )
    
    # 输出PDF，文件名去除特殊字符
    pdf_path <- file.path(gene_dir, paste0(gene, "_", gsub("[/\\:*?\"<>| ]", "_", ct), "_UMAP.pdf"))
    cat("  输出：", pdf_path, "\n")
    pdf(pdf_path, width=6, height=6)
    print(p)
    dev.off()
  }
  cat("已输出：", gene_dir, "下所有注释细胞类型UMAP图\n")
}










# 读取聚类-主细胞类型注释表
anno_path <- file.path(output_dir, "CellCluster_ManualAnno.csv")
anno_df <- read.csv(anno_path, stringsAsFactors = FALSE)
colnames(anno_df)[1:2] <- c("Cluster", "CellType") # 确保列名一致

# 自动查找所有 *_Expression_perCell.csv 文件
csv_files <- list.files(path = output_dir, pattern = "_Expression_perCell\\.csv$", full.names = TRUE)
if (length(csv_files) == 0) stop("未找到 *_Expression_perCell.csv 文件！")

for (csv_file in csv_files) {
  gene_name <- sub("_Expression_perCell\\.csv$", "", basename(csv_file))
  cat("正在处理基因：", gene_name, "\n")
  expr_df <- read.csv(csv_file, stringsAsFactors = FALSE)
  if (!all(c("Cluster", "CellType", "Gene_Expression") %in% colnames(expr_df))) {
    cat("文件", csv_file, "缺少 Cluster/CellType/Gene_Expression 列，跳过。\n")
    next
  }
  # 1. 生成聚类标签和顺序
  anno_df <- anno_df[order(as.numeric(as.character(anno_df$Cluster))), ]
  cluster_levels <- as.character(anno_df$Cluster)
  cluster_labels <- paste0("c-", anno_df$Cluster, " (", anno_df$CellType, ")")
  names(cluster_labels) <- cluster_levels
  nClusters <- length(cluster_levels)
  cluster_colors <- colorRampPalette(brewer.pal(12, "Set3"))(nClusters)
  names(cluster_colors) <- cluster_levels
  
  # 2. 统计每个聚类的表达均值和阳性比例
  bubble_df <- expr_df %>%
    group_by(Cluster) %>%
    summarise(
      avg_expr = mean(Gene_Expression, na.rm = TRUE),
      pct_pos = mean(Gene_Expression > 0, na.rm = TRUE) * 100,
      .groups = "drop"
    )
  bubble_df$Cluster <- factor(bubble_df$Cluster, levels = cluster_levels)
  
  # 3. 绘制气泡图，左边增加页边距和y轴留白
  bubble_df$Y <- " "  # 新增一列，y轴为一个空格字符串
  p <- ggplot(bubble_df, aes(x = Cluster, y = Y, size = pct_pos, color = avg_expr)) +
    geom_point() +
    scale_x_discrete(labels = cluster_labels, expand = expansion(add = 0.5)) +  # x轴两端留白
    scale_y_discrete(expand = expansion(add = 0.5)) +  # y轴两端留白
    scale_size(range = c(5, 20), name = "Percent Positive") +
    scale_color_gradientn(colors = c("blue", "white", "red"), name = "Avg Expression") +
    theme_classic(base_size = 16) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold", size = 16),  # 显示y轴
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 40)
    ) +
    labs(title = paste0(gene_name, " DotPlot by Cluster (Main Cell Type)"))
  
  # 4. 输出PDF
  pdf_name <- file.path(output_dir, paste0(gene_name, "_DotPlot_byCluster_withMainCellType.pdf"))
  pdf(pdf_name, width=14, height=7)
  print(p)
  dev.off()
  cat("已输出：", pdf_name, "\n")
}# 



# ----------- 设置输出目录和注释表路径 -----------
output_dir <- "analysis_results"  # 请根据实际情况修改
anno_path <- file.path(output_dir, "CellCluster_ManualAnno.csv")

# ----------- 读取聚类-主细胞类型注释表 -----------
anno_df <- read.csv(anno_path, stringsAsFactors = FALSE)
colnames(anno_df)[1:2] <- c("Cluster", "CellType") # 确保列名一致
anno_df <- anno_df[order(as.numeric(as.character(anno_df$Cluster))), ]
cluster_levels <- as.character(anno_df$Cluster)
cluster_labels <- paste0("c-", anno_df$Cluster, " (", anno_df$CellType, ")")
names(cluster_labels) <- cluster_levels
nClusters <- length(cluster_levels)
cluster_colors <- colorRampPalette(brewer.pal(12, "Set3"))(nClusters)
names(cluster_colors) <- cluster_levels

# ----------- 合并所有基因的表达数据 -----------
csv_files <- list.files(path = output_dir, pattern = "_Expression_perCell\\.csv$", full.names = TRUE)
if (length(csv_files) == 0) stop("未找到 *_Expression_perCell.csv 文件！")

all_bubble_df <- data.frame()
for (csv_file in csv_files) {
  gene_name <- sub("_Expression_perCell\\.csv$", "", basename(csv_file))
  expr_df <- read.csv(csv_file, stringsAsFactors = FALSE)
  if (!all(c("Cluster", "CellType", "Gene_Expression") %in% colnames(expr_df))) next
  expr_df$Gene <- gene_name
  all_bubble_df <- rbind(all_bubble_df, expr_df[, c("Cluster", "Gene_Expression", "Gene")])
}

if (nrow(all_bubble_df) == 0) stop("没有有效的表达数据可用于绘图！")

# ----------- 统计每个聚类-基因的均值和阳性比例 -----------
bubble_df <- all_bubble_df %>%
  group_by(Cluster, Gene) %>%
  summarise(
    avg_expr = mean(Gene_Expression, na.rm = TRUE),
    pct_pos = mean(Gene_Expression > 0, na.rm = TRUE) * 100,
    .groups = "drop"
  )

bubble_df$Cluster <- factor(bubble_df$Cluster, levels = cluster_levels)
gene_order <- unique(bubble_df$Gene)
bubble_df$Gene <- factor(bubble_df$Gene, levels = gene_order)

# ----------- 绘制合成气泡图，左边增加页边距 -----------
bubble_df$Y <- bubble_df$Gene  # y轴为基因名
p <- ggplot(bubble_df, aes(x = Cluster, y = Y, size = pct_pos, color = avg_expr)) +
  geom_point() +
  scale_x_discrete(labels = cluster_labels, expand = expansion(add = 0.5)) +  # x轴显示主细胞类型并留白
  scale_y_discrete(expand = expansion(add = 0.5)) +  # y轴上下留白
  scale_size(range = c(5, 20), name = "Percent Positive") +
  scale_color_gradientn(colors = c("blue", "white", "red"), name = "Avg Expression") +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold", size = 16),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 40)  # 左边页边距加大
  ) +
  labs(title = "Gene Expression DotPlot by Cluster (Main Cell Type)")

# ----------- 输出PDF -----------
pdf_name <- file.path(output_dir, "AllGenes_DotPlot_byCluster_withMainCellType.pdf")
pdf(pdf_name, width = 15, height = 4 + 0.7 * length(gene_order))
print(p)
dev.off()
cat("已输出：", pdf_name, "\n")

# 


#// 

# ----------- VlnPlot_byCluster_withMainCellType：每个基因一张，聚类顺序和标签与DotPlot一致，标签分两行防重叠 -----------

# 设置输出目录和注释表路径
output_dir <- "analysis_results"  # 请根据实际情况修改
anno_path <- file.path(output_dir, "CellCluster_ManualAnno.csv")

# 读取聚类-主细胞类型注释表
anno_df <- read.csv(anno_path, stringsAsFactors = FALSE)
colnames(anno_df)[1:2] <- c("Cluster", "CellType") # 确保列名一致
anno_df <- anno_df[order(as.numeric(as.character(anno_df$Cluster))), ]
cluster_levels <- as.character(anno_df$Cluster)
# 标签分两行，防止重叠
cluster_labels <- paste0("c-", anno_df$Cluster, "\n", anno_df$CellType)
names(cluster_labels) <- cluster_levels
nClusters <- length(cluster_levels)
cluster_colors <- colorRampPalette(brewer.pal(12, "Set3"))(nClusters)
names(cluster_colors) <- cluster_levels

# 自动查找所有 *_Expression_perCell.csv 文件
csv_files <- list.files(path = output_dir, pattern = "_Expression_perCell\\.csv$", full.names = TRUE)
if (length(csv_files) == 0) stop("未找到 *_Expression_perCell.csv 文件！")

for (csv_file in csv_files) {
  gene_name <- sub("_Expression_perCell\\.csv$", "", basename(csv_file))
  cat("正在处理基因：", gene_name, "\n")
  expr_df <- read.csv(csv_file, stringsAsFactors = FALSE)
  if (!all(c("Cluster", "CellType", "Gene_Expression") %in% colnames(expr_df))) {
    cat("文件", csv_file, "缺少 Cluster/CellType/Gene_Expression 列，跳过。\n")
    next
  }
  # 聚类顺序和标签
  expr_df$Cluster <- factor(expr_df$Cluster, levels = cluster_levels)
  expr_df$ClusterLabel <- factor(expr_df$Cluster, levels = cluster_levels, labels = cluster_labels)
  
  # 绘制小提琴图
  p_vln <- ggplot(expr_df, aes(x = ClusterLabel, y = Gene_Expression, fill = Cluster)) +
    geom_violin(scale = "width", trim = TRUE, adjust = 1, width = 0.9) +
    stat_summary(fun = median, geom = "point", shape = 23, size = 2, fill = "white", color = "black") +
    scale_fill_manual(values = cluster_colors) +
    theme_classic(base_size = 16) +
    theme(
      axis.text.x = element_text(
        angle = 45,      # 斜着显示
        hjust = 1,       # 右对齐
        vjust = 1,       # 靠下
        face = "bold",
        size = 12        # 字体适中
      ),
      axis.text.y = element_text(face = "bold"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(face = "bold"),
      strip.text = element_text(face = "bold", size = 16),
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.margin = margin(t = 5, r = 5, b = 60, l = 40)  # 底部页边距加大
    ) +
    labs(y = "Expression", title = paste0(gene_name, " Expression by Cluster (Main Cell Type)"))
  
  
  # 输出PDF
  pdf_name_vln <- file.path(output_dir, paste0(gene_name, "_VlnPlot_byCluster_withMainCellType.pdf"))
  pdf(pdf_name_vln, width = 14, height = 7)
  print(p_vln)
  dev.off()
  cat("已输出：", pdf_name_vln, "\n")
}
#########
# ----------- 合并所有基因的小提琴图（VlnPlot_byCluster_withMainCellType），标签分两行且斜着显示 -----------


# 设置输出目录和注释表路径
output_dir <- "analysis_results"  # 请根据实际情况修改
anno_path <- file.path(output_dir, "CellCluster_ManualAnno.csv")

# 读取聚类-主细胞类型注释表
anno_df <- read.csv(anno_path, stringsAsFactors = FALSE)
colnames(anno_df)[1:2] <- c("Cluster", "CellType") # 确保列名一致
anno_df <- anno_df[order(as.numeric(as.character(anno_df$Cluster))), ]
cluster_levels <- as.character(anno_df$Cluster)
# 标签分两行
cluster_labels <- paste0("c-", anno_df$Cluster, "\n", anno_df$CellType)
names(cluster_labels) <- cluster_levels
nClusters <- length(cluster_levels)
cluster_colors <- colorRampPalette(brewer.pal(12, "Set3"))(nClusters)
names(cluster_colors) <- cluster_levels

# 合并所有 *_Expression_perCell.csv
csv_files <- list.files(path = output_dir, pattern = "_Expression_perCell\\.csv$", full.names = TRUE)
if (length(csv_files) == 0) stop("未找到 *_Expression_perCell.csv 文件！")

all_vln_df <- data.frame()
for (csv_file in csv_files) {
  gene_name <- sub("_Expression_perCell\\.csv$", "", basename(csv_file))
  expr_df <- read.csv(csv_file, stringsAsFactors = FALSE)
  if (!all(c("Cluster", "CellType", "Gene_Expression") %in% colnames(expr_df))) next
  expr_df$Gene <- gene_name
  all_vln_df <- rbind(all_vln_df, expr_df[, c("Cluster", "Gene_Expression", "Gene")])
}

if (nrow(all_vln_df) == 0) stop("没有有效的表达数据可用于绘图！")

# 聚类顺序和标签
all_vln_df$Cluster <- factor(all_vln_df$Cluster, levels = cluster_levels)
all_vln_df$ClusterLabel <- factor(all_vln_df$Cluster, levels = cluster_levels, labels = cluster_labels)
gene_order <- unique(all_vln_df$Gene)
all_vln_df$Gene <- factor(all_vln_df$Gene, levels = gene_order)

# 合成小提琴图
p_vln_merge <- ggplot(all_vln_df, aes(x = ClusterLabel, y = Gene_Expression, fill = Cluster)) +
  geom_violin(scale = "width", trim = TRUE, adjust = 1, width = 0.9) +
  stat_summary(fun = median, geom = "point", shape = 23, size = 2, fill = "white", color = "black") +
  scale_fill_manual(values = cluster_colors) +
  facet_wrap(~Gene, ncol = 1, scales = "free_y") +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(
      angle = 45,      # 斜着显示
      hjust = 1,
      vjust = 1,
      face = "bold",
      size = 10
    ),
    axis.text.y = element_text(face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 16),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.margin = margin(t = 5, r = 5, b = 60, l = 40)
  ) +
  labs(y = "Expression", title = "All Genes Expression by Cluster (Main Cell Type)")

# 输出PDF
pdf_name_vln_merge <- file.path(output_dir, "AllGenes_VlnPlot_byCluster_withMainCellType.pdf")
pdf(pdf_name_vln_merge, width = 14, height = 4 + 2 * length(gene_order))
print(p_vln_merge)
dev.off()
cat("已输出：", pdf_name_vln_merge, "\n")

#// 
# ... existing code ...细胞轨迹
output_dir <- "analysis_results"
anno_path <- file.path(output_dir, "CellCluster_ManualAnno.csv")
anno_df <- read.csv(anno_path, stringsAsFactors = FALSE)
colnames(anno_df)[1:2] <- c("Cluster", "CellType")
anno_df <- anno_df[order(as.numeric(as.character(anno_df$Cluster))), ]
cluster_levels <- as.character(anno_df$Cluster)

# 获取伪时序向量
pseudotime_vec <- as.numeric(monocle3::pseudotime(cds))
names(pseudotime_vec) <- colnames(cds)

for (gene in target_genes) {
  if (!(gene %in% rownames(cds))) {
    cat("警告：基因", gene, "不在cds对象中，跳过。\n")
    next
  }
  # 合并输出到同一个文件夹
  gene_dir <- file.path(output_dir, paste0(gene, "_Trajectory_ClusterGroup_HighLow"))
  if (!dir.exists(gene_dir)) dir.create(gene_dir, recursive = TRUE)
  
  cell_meta <- as.data.frame(colData(cds))
  if (!"seurat_clusters" %in% colnames(cell_meta)) {
    stop("cds对象中没有seurat_clusters列，请检查聚类信息是否写入。")
  }
  cell_meta$Cluster <- as.character(cell_meta$seurat_clusters)
  
  for (i in seq_len(nrow(anno_df))) {
    cluster_id <- as.character(anno_df$Cluster[i])
    celltype <- as.character(anno_df$CellType[i])
    cells_in_cluster <- rownames(cell_meta)[cell_meta$Cluster == cluster_id]
    n_cells <- length(cells_in_cluster)
    if (n_cells == 0) next
    
    # 取表达量
    expr_vals <- assay(cds)[gene, cells_in_cluster]
    med <- median(expr_vals, na.rm = TRUE)
    expr_group <- ifelse(expr_vals > med, "High", "Low")
    names(expr_group) <- cells_in_cluster
    
    # ----------- 轨迹图 -----------
    highlight_vec <- rep("Other", nrow(cell_meta))
    names(highlight_vec) <- rownames(cell_meta)
    highlight_vec[cells_in_cluster] <- expr_group
    colData(cds)$ExprGroup_HighLow <- factor(highlight_vec, levels = c("Other", "Low", "High"))
    
    p_traj <- plot_cells(
      cds,
      color_cells_by = "ExprGroup_HighLow",
      show_trajectory_graph = TRUE,
      label_groups_by_cluster = FALSE,
      label_leaves = FALSE,
      label_branch_points = FALSE,
      cell_size = 1.5,
      reduction_method = "UMAP"
    ) +
      scale_color_manual(
        values = c("Other" = "grey80", "Low" = "#079EDF", "High" = "#D377A9"),
        name = "Expression Group",
        labels = c("Other", "Low (Blue)", "High (Red)")
      ) +
      labs(
        title = paste0("Gene: ", gene, " | Cluster: ", cluster_id),
        subtitle = paste0("CellType: ", celltype, " | n=", n_cells)
      ) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 14),
        legend.position = "top"
      )
    
    pdf_path_traj <- file.path(
      gene_dir,
      paste0(gene, "_Cluster", cluster_id, "_", gsub("[/\\:*?\"<>| ]", "_", celltype), "_Trajectory_HighLow.pdf")
    )
    pdf(pdf_path_traj, width = 7, height = 6)
    print(p_traj)
    dev.off()
    
    # ----------- 表达量密度图 -----------
    plot_df_expr <- data.frame(
      Expression = expr_vals,
      ExprGroup = factor(expr_group, levels = c("Low", "High"))
    )
    p_expr_density <- ggplot(plot_df_expr, aes(x = Expression, fill = ExprGroup, color = ExprGroup)) +
      geom_density(alpha = 0.35, size = 1.5, adjust = 1.1) +
      scale_fill_manual(values = c("Low" = "#079EDF", "High" = "#D377A9")) +
      scale_color_manual(values = c("Low" = "#079EDF", "High" = "#D377A9")) +
      labs(
        title = paste0("Gene: ", gene, " | Cluster: ", cluster_id),
        subtitle = paste0("CellType: ", celltype, " | n=", n_cells, " | AvgExpr=", signif(mean(expr_vals, na.rm=TRUE), 3)),
        x = "Expression",
        y = "Density",
        fill = NULL,
        color = NULL
      ) +
      theme_classic(base_size = 18) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 14),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(face = "bold", size = 14)
      )
    pdf_path_expr_density <- file.path(
      gene_dir,
      paste0(gene, "_Cluster", cluster_id, "_", gsub("[/\\:*?\"<>| ]", "_", celltype), "_ExpressionDensity_HighLow.pdf")
    )
    pdf(pdf_path_expr_density, width = 7, height = 5)
    print(p_expr_density)
    dev.off()
    
    # ----------- 伪时序密度图（如有伪时序）并输出对应表达 -----------
    pseudotime_vec_cluster <- pseudotime_vec[cells_in_cluster]
    valid_idx <- which(!is.na(pseudotime_vec_cluster) & !is.na(expr_vals))
    if (length(valid_idx) > 0) {
      pt_valid <- pseudotime_vec_cluster[valid_idx]
      expr_group_valid <- expr_group[valid_idx]
      expr_vals_valid <- expr_vals[valid_idx]
      plot_df_pt <- data.frame(
        Pseudotime = pt_valid,
        ExprGroup = factor(expr_group_valid, levels = c("Low", "High")),
        Expression = expr_vals_valid
      )
      p_pt_density <- ggplot(plot_df_pt, aes(x = Pseudotime, fill = ExprGroup, color = ExprGroup)) +
        geom_density(alpha = 0.35, size = 1.5, adjust = 1.1) +
        scale_fill_manual(values = c("Low" = "#079EDF", "High" = "#D377A9")) +
        scale_color_manual(values = c("Low" = "#079EDF", "High" = "#D377A9")) +
        labs(
          title = paste0("Gene: ", gene, " | Cluster: ", cluster_id),
          subtitle = paste0("CellType: ", celltype, " | n=", n_cells, " | Pseudotime cells: ", length(valid_idx)),
          x = "Pseudotime",
          y = "Density",
          fill = NULL,
          color = NULL
        ) +
        theme_classic(base_size = 18) +
        theme(
          plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
          plot.subtitle = element_text(hjust = 0.5, size = 14),
          legend.position = "top",
          legend.title = element_blank(),
          legend.text = element_text(face = "bold", size = 14)
        )
      pdf_path_pt_density <- file.path(
        gene_dir,
        paste0(gene, "_Cluster", cluster_id, "_", gsub("[/\\:*?\"<>| ]", "_", celltype), "_PseudotimeDensity_HighLow.pdf")
      )
      pdf(pdf_path_pt_density, width = 7, height = 5)
      print(p_pt_density)
      dev.off()
      
      # 输出伪时序-表达量对应表
      expr_table_path <- file.path(
        gene_dir,
        paste0(gene, "_Cluster", cluster_id, "_", gsub("[/\\:*?\"<>| ]", "_", celltype), "_Pseudotime_Expression.csv")
      )
      write.csv(
        data.frame(Cell = cells_in_cluster[valid_idx],
                   Pseudotime = pt_valid,
                   Expression = expr_vals_valid,
                   ExprGroup = expr_group_valid),
        file = expr_table_path,
        row.names = FALSE
      )
    }
  }
  cat("已输出：", gene_dir, "下所有聚类主注释的高低表达轨迹图、表达量密度图、伪时序密度图（如有）及伪时序-表达量表\n")
}

