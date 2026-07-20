# ==========================================================================
# 脚本名     : 单细胞CellChat细胞通讯.R
# 分类       : 08_singlecell_spatial_trajectory
# 项目来源   : 从压缩包 245.单细胞分析的细胞通讯分析.rar 整理
# 原始文件   : 245.单细胞分析的细胞通讯分析\细胞通讯-训练集.R
# 用途       : 基于已注释单细胞表达矩阵和细胞类型标签运行 CellChat，分析细胞间配体-受体通讯网络和通路强度。
# 结果图     : CellChat通讯数量circle图；通讯强度circle图；单细胞类型发送/接收网络图；LR bubble图；所有细胞类型LR气泡图
# 非肿瘤消化适配: 很适合。非肿瘤消化系统炎症、纤维化、免疫微环境文章可作为新意图模块。
# 主要 R 包  : Seurat; CellChat; dplyr; ggplot2; Matrix; patchwork; celldex; SingleR
# 整理日期   : 2026-05-13
# 备注       : 保留bioinfo-reusable-code逻辑，仅添加统一说明头；运行前请把 workDir/setwd 和输入文件名改成当前项目路径。
# ==========================================================================
# ===============================================================================
# 单细胞RNA测序 - 细胞通讯分析 (强度和数量分析专版)
# 参考细胞通讯.R的输出格式
# ===============================================================================

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

# 细胞通讯分析包
library(CellChat)         # 细胞通讯分析
library(celldex)          # SingleR注释参考数据库
library(SingleR)          # 自动细胞类型注释

message("✓ 所有必需的R包已成功加载完成！")

# ------------------- 2. 设置工作目录和参数 --------------------
workDir <- "H:\\常用分析生信\\245.单细胞分析的细胞通讯分析"
setwd(workDir)

# 分析参数
min_cells_gene     <- 5       # 至少在5个细胞出现的基因参与分析
min_genes_per_cell <- 200     # 每个细胞至少检测到的基因数
post_filter_cells  <- 300     # 二次过滤时，细胞内基因数大于此值
post_filter_mito   <- 10      # 线粒体比例不得高于此值（%）
n_top_var_features <- 2500    # 选取前2500个高度变异基因
neighbor_dims      <- 15      # 聚类与UMAP时用多少PC
cluster_resolution <- 0.3     # 聚类分辨率

# 输出目录设置
output_dir         <- "cellchat_analysis_results"
input_expr_file    <- "single_cell_data.rds.rds"

# ------------------- 3. 创建输出文件夹结构 --------------------
# 创建主输出目录
if(!dir.exists(output_dir)) dir.create(output_dir, recursive=TRUE)

cat("细胞通讯分析目录结构已创建：", output_dir, "\n")

# ===============================================================================
# 第一部分：数据读取与预处理
# ===============================================================================

cat("Step 1: 读取表达矩阵...\n")
counts <- readRDS(input_expr_file)

# 检查读取的数据结构并提取表达矩阵
cat("读取的数据类型:", class(counts), "\n")

if (inherits(counts, "Seurat")) {
  cat("检测到Seurat对象，提取counts矩阵\n")
  if (packageVersion("Seurat") >= "5.0.0") {
    cat("检测到Seurat v5，使用新的API\n")
    tryCatch({
      counts <- JoinLayers(counts)
      expr_matrix <- GetAssayData(counts, assay = "RNA", layer = "counts")
    }, error = function(e) {
      cat("无法合并layers，尝试直接提取第一个layer\n")
      expr_matrix <- LayerData(counts, assay = "RNA", layer = "counts")
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

message(sprintf("✓ 成功读取数据: %d 个基因, %d 个细胞", nrow(expr_matrix), ncol(expr_matrix)))

# ===============================================================================
# 第二部分：快速预处理用于细胞通讯分析
# ===============================================================================

cat("Step 2: 创建Seurat对象与质量控制...\n")

# 处理矩阵格式并检查NA值
if (any(is.na(expr_matrix))) {
  cat("警告：数据中包含NA值，将替换为0\n")
  expr_matrix[is.na(expr_matrix)] <- 0
}

# 创建Seurat对象
scObject <- CreateSeuratObject(
  counts = expr_matrix,
  project = "CellChatProject",
  min.cells = min_cells_gene,
  min.features = min_genes_per_cell
)

# 计算线粒体基因比例
scObject[["percent.mito"]] <- PercentageFeatureSet(scObject, pattern = "^MT-")

# 二次过滤
cat("Step 3: 二次细胞质控...\n")
cell_num_before <- ncol(scObject)
scObject <- subset(scObject, subset = nFeature_RNA > post_filter_cells & percent.mito < post_filter_mito)
cell_num_after <- ncol(scObject)
cat(sprintf("细胞过滤: 原%d，剩%d。\n", cell_num_before, cell_num_after))

if(cell_num_after < 10) stop("过滤后过少细胞，检查阈值！")

# 快速标准化与聚类
cat("Step 4: 归一化+高变基因+聚类...\n")
scObject <- NormalizeData(scObject, normalization.method="LogNormalize", scale.factor=10000)
scObject <- FindVariableFeatures(scObject, selection.method="vst", nfeatures=n_top_var_features)
scObject <- ScaleData(scObject)
scObject <- RunPCA(scObject, features=VariableFeatures(scObject))
scObject <- FindNeighbors(scObject, dims=1:neighbor_dims)
scObject <- FindClusters(scObject, resolution=cluster_resolution) 
scObject <- RunUMAP(scObject, dims=1:neighbor_dims)

# 快速细胞类型注释
cat("Step 5: SingleR细胞类型自动注释...\n")

if (!requireNamespace("SingleR", quietly = TRUE) || !requireNamespace("celldex", quietly = TRUE)) {
  cat("警告：SingleR 或 celldex 包未安装，使用聚类编号作为细胞类型\n")
  scObject$CellType <- paste0("CellType_", scObject$seurat_clusters)
} else {
  expr4singler <- GetAssayData(scObject, slot = "data")
  clusters4singler <- scObject$seurat_clusters
  
  tryCatch({
    ref_hpa <- celldex::HumanPrimaryCellAtlasData()
    singler_result <- SingleR(
      test = expr4singler,
      ref = ref_hpa,
      labels = ref_hpa$label.main,
      clusters = clusters4singler
    )
    
    scObject$CellType <- singler_result$labels[as.numeric(scObject$seurat_clusters)+1]
    
    # 保存注释结果
    write.csv(
      data.frame(Cluster=rownames(singler_result), CellType=singler_result$labels),
      file = file.path(output_dir, "CellCluster_AutoAnno_SingleR.csv"), row.names=FALSE
    )
    
  }, error = function(e) {
    cat("警告：SingleR注释失败，使用默认注释\n")
    scObject$CellType <- paste0("CellType_", scObject$seurat_clusters)
    singler_result <- NULL
  })
}

message(sprintf("✓ 细胞类型注释完成: %d 个细胞获得注释", ncol(scObject)))

# ===============================================================================
# 第三部分：细胞通讯分析 (CellChat)
# ===============================================================================

message("✓ 开始细胞通讯分析，CellChat分析")

# ----------- 1. 准备CellChat所需数据 -----------
cat("  准备CellChat分析数据...\n")

# 提取归一化的表达数据
normalized_data <- GetAssayData(scObject, assay = "RNA", slot = "data")

# 获取细胞类型信息
cell_types <- scObject$CellType
names(cell_types) <- colnames(scObject)

# 统计每种细胞类型的细胞数量
celltype_counts <- table(cell_types)
cat("  细胞类型统计：\n")
for(ct in names(celltype_counts)) {
  cat("    ", ct, ":", celltype_counts[ct], "个细胞\n")
}

# 过滤掉细胞数量过少的细胞类型（少于10个细胞）
min_cells_per_type <- 10
valid_celltypes <- names(celltype_counts)[celltype_counts >= min_cells_per_type]
valid_cells <- names(cell_types)[cell_types %in% valid_celltypes]

cat("  过滤后保留", length(valid_celltypes), "种细胞类型，", length(valid_cells), "个细胞\n")

if(length(valid_celltypes) < 2) {
  stop("  警告：有效细胞类型少于2种，无法进行细胞通讯分析\n")
} 

# 筛选数据
normalized_data_filtered <- normalized_data[, valid_cells]
cell_types_filtered <- cell_types[valid_cells]

# ----------- 2. 创建CellChat对象 -----------
cat("  创建CellChat对象...\n")

# 创建CellChat对象
meta <- data.frame(labels = cell_types_filtered, row.names = names(cell_types_filtered))
cellchat <- createCellChat(object = normalized_data_filtered, meta = meta, group.by = "labels")

# 设置数据库
cellchat@DB <- CellChatDB.human

# 预处理
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# 添加PPI网络投影（如果可用）
tryCatch({
  cellchat <- projectData(cellchat, PPI.human)
}, error = function(e) {
  cat("    PPI网络投影跳过（可能未安装相关数据）\n")
})

# 计算通讯概率和聚合网络
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

cat("  细胞通讯网络推断完成\n")

# ===============================================================================
# 第四部分：生成与参考代码相同的输出图表
# ===============================================================================

cat("  生成细胞通讯分析图表...\n")

groupSize <- as.numeric(table(cellchat@idents))

# ----------- 1. 通讯网络概览 - 数量图 -----------
cat("    生成通讯网络数量图...\n")

pdf(file.path(output_dir, "CellChat_Network_Count.pdf"), width = 8, height = 8)
par(mfrow = c(1,1), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge = T, title.name = "Number of interactions")
dev.off()

# ----------- 2. 通讯网络概览 - 权重图 -----------
cat("    生成通讯网络权重图...\n")

pdf(file.path(output_dir, "CellChat_Network_Weight.pdf"), width = 8, height = 8)
par(mfrow = c(1,1), xpd=TRUE)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge = T, title.name = "Interaction weights/strength")
dev.off()

# ----------- 3. 每个细胞类型的交互权重 - 单独保存 -----------
cat("    生成每个细胞类型的交互权重图...\n")

mat_weight <- cellchat@net$weight
for (i in 1:nrow(mat_weight)) {
  cell_type_name <- rownames(mat_weight)[i]
  safe_name <- gsub("[^A-Za-z0-9]", "_", cell_type_name)
  
  pdf(file.path(output_dir, paste0("CellChat_", safe_name, "_Individual_Weight.pdf")), width = 8, height = 8)
  par(mfrow = c(1,1), xpd=TRUE)
  mat2_weight <- matrix(0, nrow = nrow(mat_weight), ncol = ncol(mat_weight), dimnames = dimnames(mat_weight))
  mat2_weight[i, ] <- mat_weight[i, ]
  netVisual_circle(mat2_weight, vertex.weight = groupSize, weight.scale = T, 
                   edge.weight.max = max(mat_weight), label.edge = T, 
                   title.name = paste("Interaction weights -", cell_type_name))
  dev.off()
}

# ----------- 4. 每个细胞类型的交互数量 - 单独保存 -----------
cat("    生成每个细胞类型的交互数量图...\n")

mat_count <- cellchat@net$count
for (i in 1:nrow(mat_count)) {
  cell_type_name <- rownames(mat_count)[i]
  safe_name <- gsub("[^A-Za-z0-9]", "_", cell_type_name)
  
  pdf(file.path(output_dir, paste0("CellChat_", safe_name, "_Individual_Count.pdf")), width = 8, height = 8)
  par(mfrow = c(1,1), xpd=TRUE)
  mat2_count <- matrix(0, nrow = nrow(mat_count), ncol = ncol(mat_count), dimnames = dimnames(mat_count))
  mat2_count[i, ] <- mat_count[i, ]
  netVisual_circle(mat2_count, vertex.weight = groupSize, weight.scale = T, 
                   edge.weight.max = max(mat_count), label.edge = T, 
                   title.name = paste("Number of interactions -", cell_type_name))
  dev.off()
}

# ===============================================================================
# 第五部分：配体受体相互作用气泡图分析
# ===============================================================================

cat("  生成配体受体相互作用气泡图...\n")

# 获取所有细胞类型
cell_types_unique <- levels(cellchat@idents)

# 为每个细胞类型生成气泡图
for (cell_type in cell_types_unique) {
  tryCatch({
    # 提取该细胞类型相关的配体受体相互作用
    df.net <- subsetCommunication(cellchat, sources.use = cell_type)
    
    if (nrow(df.net) == 0) {
      message(sprintf("  细胞类型 %s 无相互作用数据", cell_type))
      next
    }
    
    # 准备气泡图数据
    df.net$interaction_name_2 <- paste(df.net$ligand, df.net$receptor, sep = " - ")
    df.net$prob.original <- df.net$prob
    
    # 创建细胞通信方向标签（Source→Target格式）
    df.net$cell_communication <- paste(df.net$source, df.net$target, sep = "→")
    
    # 选择前20个最强的相互作用
    top_interactions <- df.net %>%
      arrange(desc(prob)) %>%
      head(20)
    
    # 创建气泡图（横纵坐标互换，纵坐标显示通信方向）
    p <- ggplot(top_interactions, aes(x = interaction_name_2, y = cell_communication)) +
      geom_point(aes(size = prob, color = prob)) +
      scale_size_continuous(name = "Communication\nProbability", 
                            range = c(1, 8),
                            breaks = c(0.25, 0.5, 0.75, 1.0),
                            labels = c("0.25", "0.5", "0.75", "1.0")) +
      scale_color_gradient(low = "lightblue", high = "red", 
                           name = "Communication\nProbability") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
            axis.text.y = element_text(size = 9),
            axis.title = element_text(size = 12, face = "bold"),
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
            legend.title = element_text(size = 10, face = "bold"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(color = "black", size = 0.5)) +
      labs(title = paste("Ligand-Receptor Interactions\n", cell_type),
           x = "Ligand - Receptor",
           y = "Cell Communication Direction (Source→Target)") +
      guides(size = guide_legend(override.aes = list(color = "red")))
    
    # 保存PDF
    safe_name <- gsub("[^A-Za-z0-9]", "_", cell_type)
    ggsave(filename = file.path(output_dir, paste0("LR_Bubble_", safe_name, ".pdf")),
           plot = p, width = 14, height = 10, device = "pdf")
    
  }, error = function(e) {
    message(sprintf("  细胞类型 %s 生成气泡图时出错: %s", cell_type, e$message))
  })
}

# 生成所有细胞类型的汇总气泡图
tryCatch({
  df.net.all <- subsetCommunication(cellchat)
  
  if (nrow(df.net.all) > 0) {
    df.net.all$interaction_name_2 <- paste(df.net.all$ligand, df.net.all$receptor, sep = " - ")
    
    # 选择每个细胞类型组合的前5个相互作用
    top_interactions_all <- df.net.all %>%
      group_by(source, target) %>%
      top_n(5, prob) %>%
      ungroup() %>%
      arrange(desc(prob)) %>%
      head(50)  # 总体限制在50个相互作用
    
    # 创建细胞类型组合标签
    top_interactions_all$cell_pair <- paste(top_interactions_all$source, 
                                            top_interactions_all$target, sep = "→")
    
    # 汇总气泡图（横纵坐标互换）
    p_all <- ggplot(top_interactions_all, aes(x = interaction_name_2, y = cell_pair)) +
      geom_point(aes(size = prob, color = prob)) +
      scale_size_continuous(name = "Communication\nProbability", 
                            range = c(1, 6),
                            breaks = c(0.25, 0.5, 0.75, 1.0),
                            labels = c("0.25", "0.5", "0.75", "1.0")) +
      scale_color_gradient(low = "lightblue", high = "red", 
                           name = "Communication\nProbability") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
            axis.text.y = element_text(size = 8),
            axis.title = element_text(size = 12, face = "bold"),
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
            legend.title = element_text(size = 10, face = "bold"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(color = "black", size = 0.5)) +
      labs(title = "Top Ligand-Receptor Interactions - All Cell Types",
           x = "Ligand - Receptor",
           y = "Cell Communication Direction (Source→Target)")
    
    ggsave(filename = file.path(output_dir, "LR_Bubble_ALL_CellTypes.pdf"),
           plot = p_all, width = 16, height = 10, device = "pdf")
  }
  
}, error = function(e) {
  message(sprintf("  生成汇总气泡图时出错: %s", e$message))
})

# ===============================================================================
# 第六部分：生成配体受体相互作用汇总表
# ===============================================================================

cat("  生成配体受体相互作用汇总表...\n")

# 提取所有相互作用
df.net <- subsetCommunication(cellchat)

if (nrow(df.net) > 0) {
  # 添加相互作用名称和通信方向
  df.net$interaction_name <- paste(df.net$ligand, df.net$receptor, sep = " - ")
  df.net$cell_communication <- paste(df.net$source, df.net$target, sep = "→")
  
  # 选择主要列
  interaction_table <- df.net %>%
    select(source, target, cell_communication, ligand, receptor, interaction_name, 
           prob, pval, pathway_name) %>%
    arrange(desc(prob))
  
  # 保存表格
  write.csv(interaction_table, 
            file.path(output_dir, "LR_Interactions_detailed.csv"), 
            row.names = FALSE)
  
  # 创建按通路分组的汇总
  pathway_summary <- df.net %>%
    group_by(pathway_name) %>%
    summarise(
      n_interactions = n(),
      avg_prob = mean(prob),
      max_prob = max(prob),
      n_cell_pairs = n_distinct(paste(source, target)),
      .groups = "drop"
    ) %>%
    arrange(desc(avg_prob))
  
  write.csv(pathway_summary, 
            file.path(output_dir, "LR_Pathway_Summary.csv"), 
            row.names = FALSE)
}

# ===============================================================================
# 第七部分：数据导出和汇总
# ===============================================================================

cat("  导出细胞通讯分析数据...\n")

# 保存CellChat对象
saveRDS(cellchat, file.path(output_dir, "CellChat_object.rds"))

# 创建汇总报告
summary_data <- data.frame(
  Analysis = "CellChat_Analysis",
  Cell_Types = length(levels(cellchat@idents)),
  Total_Cells = length(cell_types_filtered),
  Total_Interactions = sum(cellchat@net$count),
  Signaling_Pathways = length(cellchat@netP$pathways),
  Analysis_Date = as.character(Sys.Date()),
  stringsAsFactors = FALSE
)

write.csv(summary_data, file.path(output_dir, "CellChat_Summary.csv"), row.names = FALSE)

# 保存通讯矩阵
write.csv(cellchat@net$count, file.path(output_dir, "Communication_Count_Matrix.csv"))
write.csv(cellchat@net$weight, file.path(output_dir, "Communication_Weight_Matrix.csv"))

# 计算每个细胞类型的信号统计
outgoing_signals <- apply(cellchat@net$weight, 1, sum)
incoming_signals <- apply(cellchat@net$weight, 2, sum) 

signal_stats <- data.frame(
  CellType = names(outgoing_signals),
  Outgoing_Signal_Strength = outgoing_signals,
  Incoming_Signal_Strength = incoming_signals,
  Total_Signal_Strength = outgoing_signals + incoming_signals,
  stringsAsFactors = FALSE
) %>%
arrange(desc(Total_Signal_Strength))

write.csv(signal_stats, file.path(output_dir, "Signal_Strength_Statistics.csv"), row.names = FALSE)

message("✓ 细胞通讯分析完成")

# ===============================================================================
# 分析完成
# ===============================================================================

cat("细胞通讯分析流程全部完成，结果已保存至：", output_dir, "\n")
cat("========================================\n")
