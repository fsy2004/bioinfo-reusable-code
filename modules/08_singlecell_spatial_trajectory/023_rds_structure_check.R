# =============================================================================
# 编号       : R023
# 脚本名     : 查看RDS前十行.R
# 分类       : 08_singlecell_spatial_trajectory
# 项目来源   : 网络毒理学_scTenifoldKnk_单细胞_空间转录组_免疫_诊断
# 用途       : 读取 RDS 文件并查看对象结构/前几行，作为单细胞数据检查脚本。
# 结果图     : 未检测到明确作图输出
# 主要 R 包  : 
# 整理时间   : 2026-05-10
# =============================================================================
# ============================================================================
# 查看RDS文件的前十行内容并输出为CSV
# ============================================================================

# 设置工作目录
setwd("C:/Users/fsy/Desktop/网络毒理学+单细胞+scTenifoldKnk模拟基因敲除+空间转录组+细胞轨迹+机器学习+分子对接+免疫分析+诊断模型/32.单细胞数据整理/常见的单细胞数据整理/单细胞数据整理")

# 读取RDS文件
rds_file <- "single_cell_expression_matrix.rds"
cat("正在读取文件：", rds_file, "\n\n")

data <- readRDS(rds_file)

# 检查数据类型和结构
cat("=== 数据基本信息 ===\n")
cat("数据类型：", class(data)[1], "\n")
cat("数据维度：", nrow(data), " × ", ncol(data), "\n\n")

# 根据数据类型提取前10行并保存为CSV
if(class(data)[1] == "Seurat") {
  cat("=== Seurat 对象 ===\n")
  cat("细胞数：", ncol(data), "\n")
  cat("基因数：", nrow(data), "\n\n")

  # 提取表达矩阵的前10行
  expr_matrix <- GetAssayData(data, assay = "RNA", slot = "counts")
  data_10rows <- as.data.frame(expr_matrix[1:min(10, nrow(expr_matrix)), ])

  # 添加基因名称列
  data_10rows <- cbind(Gene = rownames(data_10rows), data_10rows)

  # 保存为CSV
  output_file <- "single_cell_expression_matrix_head10.csv"
  write.csv(data_10rows, output_file, row.names = FALSE)
  cat("✓ 前10行数据已保存为：", output_file, "\n")

} else if(class(data)[1] %in% c("dgCMatrix", "matrix")) {
  cat("=== 矩阵 ===\n")
  cat("基因数：", nrow(data), "\n")
  cat("细胞数：", ncol(data), "\n\n")

  # 提取前10行
  data_10rows <- as.data.frame(data[1:min(10, nrow(data)), ])

  # 添加基因名称列
  data_10rows <- cbind(Gene = rownames(data_10rows), data_10rows)

  # 保存为CSV
  output_file <- "single_cell_expression_matrix_head10.csv"
  write.csv(data_10rows, output_file, row.names = FALSE)
  cat("✓ 前10行数据已保存为：", output_file, "\n")

} else if(is.list(data)) {
  cat("=== 列表结构 ===\n")
  cat("列表元素数：", length(data), "\n")
  cat("元素名称：", paste(names(data), collapse = ", "), "\n\n")

  # 查找矩阵元素
  possible_names <- c("expr_matrix", "counts", "expression", "data", "matrix", "expr")
  matrix_found <- FALSE

  for(name in possible_names) {
    if(name %in% names(data) && class(data[[name]])[1] %in% c("dgCMatrix", "matrix")) {
      cat("找到表达矩阵：", name, "\n\n")
      expr_matrix <- data[[name]]

      # 提取前10行
      data_10rows <- as.data.frame(expr_matrix[1:min(10, nrow(expr_matrix)), ])

      # 添加基因名称列
      data_10rows <- cbind(Gene = rownames(data_10rows), data_10rows)

      # 保存为CSV
      output_file <- "single_cell_expression_matrix_head10.csv"
      write.csv(data_10rows, output_file, row.names = FALSE)
      cat("✓ 前10行数据已保存为：", output_file, "\n")

      matrix_found <- TRUE
      break
    }
  }

  if(!matrix_found) {
    # 如果没找到特定名称，尝试第一个矩阵对象
    for(i in seq_along(data)) {
      if(class(data[[i]])[1] %in% c("dgCMatrix", "matrix")) {
        cat("找到表达矩阵（第", i, "个元素）\n\n")
        expr_matrix <- data[[i]]

        # 提取前10行
        data_10rows <- as.data.frame(expr_matrix[1:min(10, nrow(expr_matrix)), ])

        # 添加基因名称列
        data_10rows <- cbind(Gene = rownames(data_10rows), data_10rows)

        # 保存为CSV
        output_file <- "single_cell_expression_matrix_head10.csv"
        write.csv(data_10rows, output_file, row.names = FALSE)
        cat("✓ 前10行数据已保存为：", output_file, "\n")

        matrix_found <- TRUE
        break
      }
    }
  }

  if(!matrix_found) {
    cat("❌ 无法找到表达矩阵\n")
  }

} else {
  cat("=== 其他数据类型 ===\n")
  # 尝试转换为数据框并提取前10行
  data_10rows <- as.data.frame(data)
  data_10rows <- head(data_10rows, 10)

  output_file <- "single_cell_expression_matrix_head10.csv"
  write.csv(data_10rows, output_file, row.names = TRUE)
  cat("✓ 前10行数据已保存为：", output_file, "\n")
}

cat("\n=== 输出完成 ===\n")
cat("CSV文件位置：", normalizePath("single_cell_expression_matrix_head10.csv"), "\n")
