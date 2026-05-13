# =============================================================================
# 编号       : R005
# 脚本名     : OMIM和genecard数据库并集.R
# 分类       : 01_网络药理学与靶点数据库
# 项目来源   : 网络毒理学_scTenifoldKnk_单细胞_空间转录组_免疫_诊断
# 用途       : 合并多个来源的基因集合，输出并集基因并绘制集合关系图。
# 结果图     : Venn图；条形图/柱状图
# 主要 R 包  : ggplot2; ggvenn; grid; gridExtra
# 整理时间   : 2026-05-10
# =============================================================================
# ======================== OMIM和GeneCards数据库靶点取并集 ========================
# 目的：读取OMIM和GeneCards的基因列表，取并集并绘制韦恩图
# 输入：路径下的CSV和TXT文件（自动识别表头）
# 输出：并集基因列表、韦恩图、柱状图
# ==============================================================================

library(ggvenn)
library(gridExtra)
library(ggplot2)
library(grid)

set.seed(12345)

# 设置工作路径
work_dir <- "C:/Users/fsy/Desktop/网络毒理学+单细胞+scTenifoldKnk模拟基因敲除+空间转录组+细胞轨迹+机器学习+分子对接+免疫分析+诊断模型/08.OMIM和genecard数据库并集"
setwd(work_dir)

cat("================== OMIM和GeneCards靶点取并集开始 ==================\n\n")

# 创建输出文件夹
outputDir <- "output_folder"
if (!dir.exists(outputDir)) dir.create(outputDir)

# =============== 辅助函数：自动识别表头并读取基因 ===============
read_genes <- function(filepath) {
  # 读取第一行判断是否有表头
  first_line <- readLines(filepath, n = 1, warn = FALSE)
  # 取第一个字段（兼容csv和txt）
  first_field <- trimws(unlist(strsplit(first_line, "[,\t]"))[1])
  # 去掉引号
  first_field <- gsub("\"", "", first_field)

  # 判断：包含gene（不区分大小写）就认为是表头
  has_header <- grepl("gene", first_field, ignore.case = TRUE)

  cat("  文件：", basename(filepath), "\n")
  cat("  第一行：", first_line, "\n")
  cat("  识别为表头：", ifelse(has_header, "是", "否"), "\n")

  # 根据文件扩展名选择读取方式
  ext <- tolower(tools::file_ext(filepath))

  if (ext == "csv") {
    data <- read.csv(filepath, header = has_header, stringsAsFactors = FALSE,
                     fileEncoding = "UTF-8")
  } else {
    data <- read.table(filepath, header = has_header, sep = "\t",
                       stringsAsFactors = FALSE, fill = TRUE, quote = "")
  }

  # 取第一列作为基因
  genes <- as.character(data[, 1])
  genes <- trimws(genes)
  genes <- genes[!is.na(genes) & genes != ""]
  genes <- unique(genes)

  cat("  基因数量：", length(genes), "\n\n")
  return(genes)
}

# =============== 读取所有基因文件 ===============
cat("步骤1: 读取基因文件...\n\n")

# 获取所有csv和txt文件
all_files <- list.files(path = work_dir, pattern = "\\.(csv|txt)$",
                        full.names = TRUE, ignore.case = TRUE)

if (length(all_files) == 0) {
  stop("错误：未找到CSV或TXT文件！")
}

geneList <- list()

for (f in all_files) {
  fname <- basename(f)
  # 用文件名（去扩展名）作为集合名
  set_name <- tools::file_path_sans_ext(fname)
  geneList[[set_name]] <- read_genes(f)
}

cat("共读取", length(geneList), "个基因集合\n\n")

# =============== 计算并集 ===============
cat("步骤2: 计算并集...\n")

unionGenes <- Reduce(union, geneList)
unionCount <- length(unionGenes)

cat("并集基因数量：", unionCount, "\n\n")

# =============== 保存并集基因列表 ===============
cat("步骤3: 保存并集基因列表...\n")

write.csv(data.frame(Gene = unionGenes),
          file = file.path(outputDir, "Disease_Union_Genes.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

cat("已保存：", file.path(outputDir, "Disease_Union_Genes.csv"), "\n\n")

# =============== 绘制韦恩图 ===============
cat("步骤4: 绘制韦恩图...\n")

pdf(file = file.path(outputDir, "venn.pdf"), width = 6, height = 6)
venn_plot <- ggvenn(
  geneList,
  show_percentage = TRUE,
  stroke_color = "white",
  stroke_size = 0.5,
  fill_color = c("#FFA700", "#1E90FF")[seq_along(geneList)],
  set_name_color = c("#FFA700", "#1E90FF")[seq_along(geneList)],
  set_name_size = 6,
  text_size = 4.5
)
grid.arrange(venn_plot,
             textGrob(paste("Union Gene Count:", unionCount), gp = gpar(fontsize = 14)),
             ncol = 1, heights = c(5, 1))
dev.off()

cat("已保存：", file.path(outputDir, "venn.pdf"), "\n\n")

# =============== 绘制柱状图 ===============
cat("步骤5: 绘制集合大小柱状图...\n")

setSizes <- sapply(geneList, length)
pdf(file = file.path(outputDir, "set_barplot.pdf"), width = 6, height = 4)
print(
  ggplot(data.frame(Set = names(setSizes), Size = setSizes),
         aes(x = Set, y = Size, fill = Set)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    ylab("Gene Count") + xlab("") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
)
dev.off()

cat("已保存：", file.path(outputDir, "set_barplot.pdf"), "\n\n")

# =============== 汇总报告 ===============
cat("================== 分析汇总 ==================\n")
for (set_name in names(geneList)) {
  cat("  -", set_name, "：", length(geneList[[set_name]]), "个基因\n")
}
cat("并集基因数量：", unionCount, "\n")
cat("输出文件夹：", outputDir, "\n")
cat("==========================================\n\n")

cat("处理完成！\n")
