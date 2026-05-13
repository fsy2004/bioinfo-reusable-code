# =============================================================================
# 编号       : R025
# 脚本名     : 代码.R
# 分类       : 08_单细胞_空间转录组_细胞轨迹
# 项目来源   : 网络毒理学_scTenifoldKnk_单细胞_空间转录组_免疫_诊断
# 用途       : 读取并整理单细胞数据，构建/保存 Seurat 或 RDS 分析对象。
# 结果图     : 未检测到明确作图输出
# 主要 R 包  : dplyr; Matrix; Seurat; stringr
# 整理时间   : 2026-05-10
# =============================================================================
# =============================================================================
# 1. 环境准备
# =============================================================================

# 加载必要的库
library(Seurat)
library(Matrix)
library(dplyr)
library(stringr)

# 设置工作目录
work_dir <- "C:\\Users\\wo\\Desktop\\充电视频\\397.肺腺癌单细胞数据整理\\GSE123902_RAW"
setwd(work_dir)
cat("工作目录设置为:", getwd(), "\n")

# =============================================================================
# 2. 文件检查和样本信息提取
# =============================================================================

# 获取所有dense.csv.gz文件
csv_files <- list.files(pattern = "_dense\\.csv\\.gz$", full.names = FALSE)

if(length(csv_files) == 0) {
  stop("未找到_dense.csv.gz文件，请检查路径")
}

cat("找到", length(csv_files), "个CSV文件:\n")
print(csv_files)

# 提取样本信息（从文件名中提取GSM ID）
sample_info <- data.frame(
  file_name = csv_files,
  gsm_id = gsub("_.*", "", csv_files),  # 提取GSM ID
  stringsAsFactors = FALSE
)

print(sample_info)

# =============================================================================
# 3. 批量读取和处理CSV文件
# =============================================================================

cat("\n开始批量处理CSV文件...\n")

# 初始化存储列表
seurat_list <- list()
processing_log <- list()

# 批量处理函数
process_csv_file <- function(file_name, gsm_id) {
  cat("正在处理:", gsm_id, "...\n")

  start_time <- Sys.time()

  tryCatch({
    # 读取CSV.GZ文件
    data <- read.csv(gzfile(file_name), header = TRUE, check.names = FALSE, row.names = 1)

    cat("  - 原始细胞数:", nrow(data), "\n")
    cat("  - 原始基因数:", ncol(data), "\n")

    # 转置矩阵：行变为基因，列变为细胞
    counts <- t(as.matrix(data))

    # 处理重复的基因名
    gene_names <- rownames(counts)
    if(any(duplicated(gene_names))) {
      cat("  - 发现重复基因名，正在处理...\n")
      # 使用make.unique处理重复基因名
      gene_names <- make.unique(gene_names, sep = "_")
      rownames(counts) <- gene_names
      cat("  - 重复基因名已处理\n")
    }

    # 获取原始细胞barcode（列名）
    original_barcodes <- colnames(counts)

    # 重命名细胞：GSM_barcode格式
    new_cell_names <- paste0(gsm_id, "_", sprintf("%.0f", as.numeric(original_barcodes)))
    colnames(counts) <- new_cell_names

    # 转换为稀疏矩阵
    counts_sparse <- as(counts, "dgCMatrix")

    # 创建Seurat对象
    seurat_obj <- CreateSeuratObject(
      counts = counts_sparse,
      project = gsm_id,
      min.cells = 0,
      min.features = 0
    )

    # 计算线粒体基因比例
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

    # 计算核糖体基因比例
    seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")

    # 添加样本信息到metadata
    seurat_obj$sample_id <- gsm_id
    seurat_obj$orig.file <- file_name

    end_time <- Sys.time()
    processing_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

    # 记录处理信息
    log_entry <- list(
      sample_id = gsm_id,
      genes = nrow(seurat_obj),
      cells = ncol(seurat_obj),
      processing_time = round(processing_time, 2),
      status = "Success"
    )

    cat("  完成! 基因数:", nrow(seurat_obj),
        ", 细胞数:", ncol(seurat_obj), "\n")
    cat("  处理时间:", round(processing_time, 2), "秒\n\n")

    return(list(seurat = seurat_obj, log = log_entry))

  }, error = function(e) {
    cat("  处理失败:", e$message, "\n\n")

    log_entry <- list(
      sample_id = gsm_id,
      error_message = e$message,
      status = "Failed"
    )

    return(list(seurat = NULL, log = log_entry))
  })
}

# 执行批量处理
for(i in 1:nrow(sample_info)) {
  result <- process_csv_file(sample_info$file_name[i], sample_info$gsm_id[i])

  if(!is.null(result$seurat)) {
    seurat_list[[sample_info$gsm_id[i]]] <- result$seurat
  }

  processing_log[[i]] <- result$log
}

# =============================================================================
# 4. 数据合并
# =============================================================================

if(length(seurat_list) == 0) {
  stop("没有成功处理任何样本数据")
}

cat("开始合并", length(seurat_list), "个样本的数据...\n")

# 合并所有成功处理的样本
if(length(seurat_list) == 1) {
  merged_seurat <- seurat_list[[1]]
} else {
  merged_seurat <- merge(
    seurat_list[[1]],
    y = seurat_list[-1],
    add.cell.ids = NULL,
    project = "GSE123902"
  )
}

# =============================================================================
# 5. 合并后的数据统计
# =============================================================================

cat("数据合并完成！\n")
cat("数据统计:\n")
cat("  - 总样本数:", length(unique(merged_seurat$sample_id)), "\n")
cat("  - 总细胞数:", ncol(merged_seurat), "\n")
cat("  - 总基因数:", nrow(merged_seurat), "\n")

# 样本统计
sample_stats <- table(merged_seurat$sample_id)
cat("\n各样本细胞数分布:\n")
print(sample_stats)

# 基础统计
meta_data <- merged_seurat@meta.data
basic_stats <- aggregate(
  cbind(nFeature_RNA, nCount_RNA, percent.mt, percent.ribo) ~ sample_id,
  data = meta_data,
  FUN = function(x) c(mean = mean(x), median = median(x))
)

# 计算每个样本的细胞数
cell_counts <- as.data.frame(table(meta_data$sample_id))
colnames(cell_counts) <- c("sample_id", "n_cells")

# 整理统计结果
basic_stats_df <- data.frame(
  sample_id = cell_counts$sample_id,
  n_cells = cell_counts$n_cells,
  mean_genes = round(basic_stats$nFeature_RNA[, "mean"], 0),
  median_genes = round(basic_stats$nFeature_RNA[, "median"], 0),
  mean_UMI = round(basic_stats$nCount_RNA[, "mean"], 0),
  mean_mt_pct = round(basic_stats$percent.mt[, "mean"], 2),
  mean_ribo_pct = round(basic_stats$percent.ribo[, "mean"], 2)
)

print(basic_stats_df)

# =============================================================================
# 6. 数据保存
# =============================================================================

cat("\n保存数据文件...\n")

# 保存主要RDS文件
output_file <- "scRNA_matrix.rds"
saveRDS(merged_seurat, file = output_file)
cat("主数据文件已保存:", output_file, "\n")

# 保存处理日志
log_df <- do.call(rbind, lapply(processing_log, function(x) {
  missing_cols <- setdiff(c("sample_id", "genes", "cells", "processing_time", "status", "error_message"),
                          names(x))
  for(col in missing_cols) {
    x[[col]] <- NA
  }
  return(as.data.frame(x, stringsAsFactors = FALSE))
}))
write.csv(log_df, "processing_log.csv", row.names = FALSE)

# 保存基础统计
write.csv(basic_stats_df, "basic_statistics.csv", row.names = FALSE)

# 保存样本统计
sample_summary <- data.frame(
  sample_id = names(sample_stats),
  cell_count = as.numeric(sample_stats),
  percentage = round(as.numeric(sample_stats)/sum(sample_stats)*100, 2)
)
write.csv(sample_summary, "sample_summary.csv", row.names = FALSE)

# =============================================================================
# 7. 提取并保存前10个基因的表达矩阵
# =============================================================================

cat("\n提取前10个基因的表达数据...\n")

tryCatch({
  # 获取counts矩阵
  seurat_version <- packageVersion("Seurat")

  if(seurat_version >= "5.0.0") {
    if("LayerData" %in% ls("package:Seurat")) {
      counts_matrix <- LayerData(merged_seurat, assay = "RNA", layer = "counts")
    } else {
      counts_matrix <- merged_seurat[["RNA"]]@layers$counts
      rownames(counts_matrix) <- rownames(merged_seurat[["RNA"]])
      colnames(counts_matrix) <- colnames(merged_seurat[["RNA"]])
    }
  } else {
    counts_matrix <- GetAssayData(merged_seurat, assay = "RNA", slot = "counts")
  }

  # 提取前10个基因
  n_genes <- min(10, nrow(counts_matrix))
  first_genes_matrix <- counts_matrix[1:n_genes, ]

  # 转换为数据框
  first_genes_df <- as.data.frame(as.matrix(first_genes_matrix))
  first_genes_df <- data.frame(
    Gene = rownames(first_genes_df),
    first_genes_df,
    stringsAsFactors = FALSE
  )

  write.csv(first_genes_df, "first_10_genes_expression.csv", row.names = FALSE)
  cat("前10个基因的表达矩阵已保存\n")

}, error = function(e) {
  cat("提取表达矩阵时出错:", e$message, "\n")
})

# 保存前10个细胞的metadata
first_cells_meta <- merged_seurat@meta.data[1:min(10, ncol(merged_seurat)), ]
first_cells_meta <- data.frame(
  Cell_ID = rownames(first_cells_meta),
  first_cells_meta,
  stringsAsFactors = FALSE
)
write.csv(first_cells_meta, "first_10_cells_metadata.csv", row.names = FALSE)

# =============================================================================
# 8. 完成
# =============================================================================

cat("\n所有处理完成！\n")
cat("\n各样本详细信息:\n")
for(i in 1:nrow(basic_stats_df)) {
  cat(sprintf("样本 %s: %d细胞, %d平均基因数, %.2f%%线粒体基因\n",
              basic_stats_df$sample_id[i],
              basic_stats_df$n_cells[i],
              basic_stats_df$mean_genes[i],
              basic_stats_df$mean_mt_pct[i]))
}

cat("\n已生成的文件列表:\n")
cat("1. RDS文件:", output_file, "\n")
cat("2. CSV文件:\n")
cat("   - processing_log.csv\n")
cat("   - basic_statistics.csv\n")
cat("   - sample_summary.csv\n")
cat("   - first_10_genes_expression.csv\n")
cat("   - first_10_cells_metadata.csv\n")
