# =============================================================================
# 编号       : R024
# 脚本名     : 单细胞数据整理rds文件.R
# 分类       : 08_singlecell_spatial_trajectory
# 项目来源   : 网络毒理学_scTenifoldKnk_单细胞_空间转录组_免疫_诊断
# 用途       : 读取并整理单细胞数据，构建/保存 Seurat 或 RDS 分析对象。
# 结果图     : 未检测到明确作图输出
# 主要 R 包  : celldex; dplyr; ggplot2; limma; magrittr; progress; RColorBrewer; Seurat; SingleR
# 整理时间   : 2026-05-10
# =============================================================================
# 安装所需的包
#if (!requireNamespace("progress", quietly = TRUE)) {
#  install.packages("progress")
#}

# 导入必要的库
library(progress)      # 用于创建进度条
library(Seurat)        # 用于单细胞分析
library(limma)         # 用于差异分析
library(dplyr)         # 用于数据处理
library(magrittr)      # 用于管道操作符（%>%）
library(celldex)       # 用于细胞类型注释
library(SingleR)       # 用于单细胞注释
library(ggplot2)       # 用于数据可视化
library(RColorBrewer)  # 用于生成调色板

# 2. 设置工作目录
workDir <- "C:/Users/fsy/Desktop/网络毒理学+单细胞+scTenifoldKnk模拟基因敲除+空间转录组+细胞轨迹+机器学习+分子对接+免疫分析+诊断模型/32.单细胞数据整理/常见的单细胞数据整理/单细胞数据整理"  # 设置工作目录为存储单细胞数据的文件夹
setwd(workDir)  # 更改当前工作目录为上述指定的路径

# ===================== 自定义函数：读取Matrix Market格式数据 =====================
read_matrix_market <- function(data_dir, sample_name) {
  tryCatch({
    # 检查所需文件是否存在
    barcode_file <- file.path(data_dir, "barcodes.tsv.gz")
    gene_file <- file.path(data_dir, "genes.tsv.gz")
    matrix_file <- file.path(data_dir, "matrix.mtx.gz")

    if (!all(file.exists(barcode_file, gene_file, matrix_file))) {
      stop(paste("缺少必要的文件。需要：barcodes.tsv.gz, genes.tsv.gz, matrix.mtx.gz"))
    }

    # 读取数据
    barcodes <- read.csv(barcode_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    genes <- read.csv(gene_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    matrix <- Matrix::readMM(matrix_file)

    # 处理基因名称：如果有多列，使用第一列作为基因ID，第二列作为基因名
    if (ncol(genes) >= 2) {
      gene_names <- genes[, 2]  # 使用第二列作为基因名
    } else {
      gene_names <- genes[, 1]  # 如果只有一列，使用该列
    }

    # 设置行名和列名
    rownames(matrix) <- gene_names

    # 为细胞条形码添加样本前缀（使用下划线分隔）
    cell_barcodes <- paste0(sample_name, "_", barcodes[, 1])
    colnames(matrix) <- cell_barcodes

    # 检查重复的基因名
    if (anyDuplicated(rownames(matrix)) > 0) {
      duplicated_genes <- rownames(matrix)[duplicated(rownames(matrix))]
      message(paste("  ⚠ 样本内发现", length(unique(duplicated_genes)), "个重复基因，保留第一个"))
      # 只保留第一次出现的基因
      matrix <- matrix[!duplicated(rownames(matrix)), ]
    }

    message(paste("成功读取样本:", sample_name, "维度:", nrow(matrix), "x", ncol(matrix)))
    return(list(matrix = matrix, sample_name = sample_name))

  }, error = function(e) {
    message(paste("读取样本", sample_name, "时发生错误：", e$message))
    return(NULL)
  })
}

# 3. 自动检测输入方式：单个样本 vs 多个样本目录
cat("========== 数据输入方式自动检测 ==========\n")

# 检查当前目录是否有必要的文件
has_single_sample <- all(file.exists("barcodes.tsv.gz", "genes.tsv.gz", "matrix.mtx.gz"))

# 检查是否有子目录
dirs <- list.dirs(workDir)
dirs_sample <- dirs[-1]
has_subdirs <- length(dirs_sample) > 0

# 确定输入方式
if (has_single_sample) {
  input_type <- "1"
  cat("检测到：单个样本目录模式\n")
  cat("数据位置: ", workDir, "\n")
} else if (has_subdirs) {
  input_type <- "2"
  cat("检测到：多样本子目录模式\n")
  cat("找到", length(dirs_sample), "个样本子目录\n")
} else {
  stop("未找到数据！请检查:\n1. 当前目录是否包含 barcodes.tsv.gz, genes.tsv.gz, matrix.mtx.gz\n2. 或检查是否有样本子目录")
}

if (input_type == "1") {
  # 单样本模式
  cat("\n进入单样本处理模式\n")

  sample_result <- read_matrix_market(workDir, basename(workDir))

  if (!is.null(sample_result)) {
    counts <- list(sample_result$matrix)
    names(counts) <- sample_result$sample_name
  } else {
    stop("无法读取样本数据！")
  }

} else if (input_type == "2") {
  # 多样本模式
  cat("\n进入多样本处理模式\n")

  # 提取子目录名称作为样本的名称
  names(dirs_sample) <- gsub(".+\\/(.+)", "\\1", dirs_sample)

  # 创建进度条
  progress_bar <- progress::progress_bar$new(
    format = " 进度 [:bar] :percent 完成 时间: :elapsedfull 剩余时间: :eta",
    total = length(dirs_sample),
    clear = FALSE,
    width = 60
  )

  # 初始化列表和控制变量
  counts <- list()
  skip_to_next <- FALSE
  user_abort <- FALSE

  # 逐一读取每个样本目录的数据
  for (i in seq_along(dirs_sample)) {
    if (user_abort) break  # 如果用户选择中止，跳出循环

    sample_name <- names(dirs_sample)[i]
    message(paste("\n[", i, "/", length(dirs_sample), "] 正在读取样本：", sample_name))

    # 尝试读取数据
    sample_result <- read_matrix_market(dirs_sample[i], sample_name)

    if (!is.null(sample_result)) {
      counts[[length(counts) + 1]] <- sample_result$matrix
      names(counts)[length(counts)] <- sample_name
      message("✓ 成功读取")
    } else {
      message("✗ 读取失败")
      message("请选择:")
      message("  y  - 继续处理下一个样本")
      message("  n  - 跳过此样本")
      message("  quit - 中止数据处理")

      # 修复：使用更可靠的输入方法
      user_choice <- tolower(trimws(readline("你的选择: ")))

      if (user_choice == "quit") {
        user_abort <- TRUE
        message("中止数据处理")
        break
      } else if (user_choice == "n") {
        message("跳过此样本")
      }
      # 其他情况（包括 "y" 或其他输入）都继续
      skip_to_next <- TRUE
    }

    progress_bar$tick()  # 更新进度条
  }

}

# 6. 合并所有读取的数据
if (length(counts) == 0) {
  stop("没有成功读取任何数据，请检查输入路径和数据格式！")
}

message(paste("\n成功读取", length(counts), "个样本"))

# 尝试合并数据
counts_combined <- tryCatch({
  if (length(counts) == 1) {
    # 单样本情况，直接使用
    counts[[1]]
  } else {
    # 多样本情况，按列合并
    do.call(cbind, counts)
  }
}, error = function(e) {
  message(paste("合并数据时发生错误：", e$message))
  return(NULL)
})

if (is.null(counts_combined)) {
  stop("数据合并失败！")
}

# 7. 合并后处理：详细调试
message("\n=============== 数据质量检查与调试 ===============")

message("\n[1] 检查矩阵维度：")
message(paste("  行数（基因）：", nrow(counts_combined)))
message(paste("  列数（细胞）：", ncol(counts_combined)))

message("\n[2] 检查行名（基因名）重复：")
dup_rownames <- sum(duplicated(rownames(counts_combined)))
message(paste("  重复数：", dup_rownames))
if (dup_rownames > 0) {
  message("  重复的基因:")
  print(table(rownames(counts_combined))[table(rownames(counts_combined)) > 1])
}

message("\n[3] 检查列名（细胞条形码）重复：")
dup_colnames <- sum(duplicated(colnames(counts_combined)))
message(paste("  重复数：", dup_colnames))
if (dup_colnames > 0) {
  message("  ⚠ 发现重复的细胞条形码，正在处理...")
  message(paste("  重复的条形码数量：", dup_colnames))

  # 保留第一个，删除后续重复的列
  before_col <- ncol(counts_combined)
  counts_combined <- counts_combined[, !duplicated(colnames(counts_combined))]
  after_col <- ncol(counts_combined)
  message(paste("  ✓ 处理完成：", before_col, "→", after_col, "个细胞"))
}

message("\n[4] 统一处理特征名中的下划线：")
rownames(counts_combined) <- gsub("_", "-", rownames(counts_combined))
message("  ✓ 下划线已转换为破折号")

message("\n[5] 再次检查行名重复：")
dup_rownames_after <- sum(duplicated(rownames(counts_combined)))
message(paste("  重复数：", dup_rownames_after))
if (dup_rownames_after > 0) {
  message("  ⚠ 发现重复的基因名，正在删除...")
  before_count <- nrow(counts_combined)
  counts_combined <- counts_combined[!duplicated(rownames(counts_combined)), ]
  after_count <- nrow(counts_combined)
  message(paste("  ✓ 处理完成：", before_count, "→", after_count, "个基因"))
}

message("\n[6] 最终验证：")
message(paste("  行名唯一性：", if(anyDuplicated(rownames(counts_combined)) == 0) "✓ 通过" else "✗ 失败"))
message(paste("  列名唯一性：", if(anyDuplicated(colnames(counts_combined)) == 0) "✓ 通过" else "✗ 失败"))
message(paste("  矩阵类型：", class(counts_combined)))
message(paste("  最终维度：", nrow(counts_combined), "x", ncol(counts_combined)))

# 8. 创建Seurat对象，用于后续分析
message("\n创建Seurat对象...")
tryCatch({
  Peripheral_Blood_Mononuclear_Cells <- CreateSeuratObject(
    counts_combined,
    min.cells = 10,
    min.features = 40
  )
  message("✓ Seurat对象创建成功")
  message(paste("  包含", ncol(Peripheral_Blood_Mononuclear_Cells), "个细胞"))
  message(paste("  包含", nrow(Peripheral_Blood_Mononuclear_Cells), "个基因"))
}, error = function(e) {
  message(paste("✗ 创建Seurat对象时发生错误：", e$message))
  message("\n DEBUG信息：")
  message(paste("  行名唯一值数：", length(unique(rownames(counts_combined)))))
  message(paste("  实际行数：", nrow(counts_combined)))
  message(paste("  列名唯一值数：", length(unique(colnames(counts_combined)))))
  message(paste("  实际列数：", ncol(counts_combined)))

  # 检查是否有NA值
  message(paste("  行名中是否有NA：", any(is.na(rownames(counts_combined)))))
  message(paste("  列名中是否有NA：", any(is.na(colnames(counts_combined)))))

  stop(e$message)
})

# 9. 获取基因表达数据
message("提取基因表达数据...")
# ===================== 修复版：提取基因表达矩阵 =====================
tryCatch({
  # Seurat 5 新版写法 ✅
  single_cell_data <- GetAssayData(
    object = Peripheral_Blood_Mononuclear_Cells,
    layer = "counts"   # 把 slot 换成 layer 即可
  )
  message("✓ 表达数据提取成功！")
}, error = function(e) {
  message(paste("提取表达数据时发生错误：", e$message))
  stop(e$message)
})

# 9. 输出数据到RDS文件
message("\n正在导出数据为 RDS 文件...")

tryCatch({
  # 保存原始表达矩阵为RDS文件
  output_file <- "single_cell_expression_matrix.rds"
  saveRDS(single_cell_data, file = output_file)
  message(paste("✓ 表达矩阵 RDS 文件导出成功！"))
  message(paste("  文件名：", output_file))
}, error = function(e) {
  message(paste("导出RDS文件时发生错误：", e$message))
})

# 11. 输出处理完成的消息
message(paste(rep("=", 50), collapse = ""))
message("✓ 数据导出完成！生成的文件包括：")
message("1. single_cell_expression_matrix.rds - 表达矩阵RDS文件")
message(paste(rep("=", 50), collapse = ""))

# 12. 输出数据矩阵的维度信息（行数和列数）
cat("\n=== 数据矩阵基本信息 ===\n")
cat("数据矩阵的维度为：", dim(single_cell_data)[1], "行（基因数）", dim(single_cell_data)[2], "列（细胞数）\n")

# 13. 输出文件大小信息
rds_matrix_file <- "single_cell_expression_matrix.rds"

if (file.exists(rds_matrix_file)) {
  rds_matrix_size <- file.size(rds_matrix_file) / (1024^2)  # 转换为MB
  cat("表达矩阵RDS文件大小：", round(rds_matrix_size, 2), "MB\n")
}

# 14. 快速数据质量检查
message("\n=== 数据质量快速检查 ===")
cat("- 基因总数：", nrow(single_cell_data), "\n")
cat("- 细胞总数：", ncol(single_cell_data), "\n")
cat("- 非零表达值数量：", sum(single_cell_data > 0), "\n")
sparsity <- round((1 - sum(single_cell_data > 0) / (nrow(single_cell_data) * ncol(single_cell_data))) * 100, 2)
cat("- 稀疏度：", sparsity, "%\n")

# 15. 输出样本信息汇总
message("\n=== 处理的样本信息 ===")
if (length(counts) > 0) {
  for (i in seq_along(names(counts))) {
    cat(paste(i, ". ", names(counts)[i], "\n"))
  }
}

message("\n✓ 全部处理完成！")
