# =============================================================================
# 编号       : R015
# 脚本名     : venn交集.R
# 分类       : 04_机器学习筛选特征基因
# 项目来源   : 网络毒理学_scTenifoldKnk_单细胞_空间转录组_免疫_诊断
# 用途       : 计算多个基因集合交集/两两交集，并生成 Venn/集合图与交集基因表。
# 结果图     : Venn图；UpSet图
# 主要 R 包  : ComplexUpset; ggplot2; ggvenn; RColorBrewer; viridis
# 整理时间   : 2026-05-10
# =============================================================================
library(ggvenn)
library(RColorBrewer)
library(viridis)
library(ComplexUpset)
library(ggplot2)

# 1. 设置工作目录
setwd("C:/Users/fsy/Desktop/网络毒理学+单细胞+scTenifoldKnk模拟基因敲除+空间转录组+细胞轨迹+机器学习+分子对接+免疫分析+诊断模型/19.三个机器学习方法取交集")

# 2. 创建输出文件夹
output_folder <- "output_folder"
if (!dir.exists(output_folder)) dir.create(output_folder)

# 3. 获取所有TXT和CSV文件
all_files <- c(list.files(pattern = "\\.txt$", full.names = TRUE),
               list.files(pattern = "\\.csv$", full.names = TRUE))

# 4. 读取所有文件并存入列表（自动检测表头）
gene_list <- list()

for (file in all_files) {
  file_name <- tools::file_path_sans_ext(basename(file))
  file_ext <- tolower(tools::file_ext(file))

  if (file_ext == "csv") {
    # 读取CSV文件
    rt <- read.csv(file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    # 检查第一列列名是否为gene（不区分大小写）
    if (tolower(colnames(rt)[1]) == "gene") {
      geneNames <- trimws(rt[, 1])
    } else {
      geneNames <- trimws(rt[, 1])
    }
  } else {
    # 读取TXT文件，先读取第一行判断是否为表头
    firstLine <- trimws(readLines(file, n = 1, warn = FALSE))
    if (tolower(firstLine) == "gene") {
      # 有表头，跳过第一行
      rt <- readLines(file, warn = FALSE)[-1]
    } else {
      rt <- readLines(file, warn = FALSE)
    }
    geneNames <- trimws(rt)
  }

  # 去除空值和重复
  geneNames <- geneNames[geneNames != ""]
  geneNames <- unique(geneNames)

  gene_list[[file_name]] <- geneNames
  cat(file_name, "基因数：", length(geneNames), "\n")
}

# 5. 利用 viridis 调色板生成颜色
nColors <- length(gene_list)
myColors <- viridis(nColors)

# 6. 绘制Venn图并保存为PDF
pdf(file = file.path(output_folder, "venn.pdf"), width = 10, height = 10)

ggvenn(gene_list, show_percentage = TRUE,
       stroke_color = "white", stroke_size = 1.5,
       fill_color = myColors,
       set_name_color = "black",
       set_name_size = 8,
       text_size = 6,
       text_color = "black")

dev.off()

cat("Venn图保存至：", file.path(output_folder, "venn.pdf"), "\n")

# 7. 获取所有文件的全局交集基因并保存（CSV格式）
intersect_genes <- Reduce(intersect, gene_list)
cat("全局交集基因数：", length(intersect_genes), "\n")

# 将全局交集转为数据框保存
global_intersect_df <- data.frame(Gene = intersect_genes, stringsAsFactors = FALSE)
write.csv(global_intersect_df,
          file = file.path(output_folder, "IntersectionGenes.csv"),
          row.names = FALSE,
          quote = FALSE)

cat("全局交集基因已保存至CSV文件：", file.path(output_folder, "IntersectionGenes.csv"), "\n")

# 8. 生成每两个文件之间的交集基因列表，并保存为CSV表格
pairwise_intersections <- data.frame(
  File1 = character(),
  File2 = character(),
  Intersection_Count = integer(),
  Intersection_Genes = character(),
  stringsAsFactors = FALSE
)

files <- names(gene_list)
for (i in 1:(length(files) - 1)) {
  for (j in (i + 1):length(files)) {
    file1 <- files[i]
    file2 <- files[j]

    common_genes <- intersect(gene_list[[file1]], gene_list[[file2]])

    pairwise_intersections <- rbind(pairwise_intersections,
                                    data.frame(
                                      File1 = file1,
                                      File2 = file2,
                                      Intersection_Count = length(common_genes),
                                      Intersection_Genes = paste(common_genes, collapse = ";"),
                                      stringsAsFactors = FALSE
                                    ))
  }
}

write.csv(pairwise_intersections,
          file = file.path(output_folder, "PairwiseIntersectionGenes.csv"),
          row.names = FALSE,
          quote = FALSE)

cat("每两个文件间的交集基因表格已保存至CSV文件：", file.path(output_folder, "PairwiseIntersectionGenes.csv"), "\n")


cat("\n所有分析完成！\n")
