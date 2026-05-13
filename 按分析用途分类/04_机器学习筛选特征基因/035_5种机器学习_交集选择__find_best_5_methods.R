# =============================================================================
# 编号       : R035
# 脚本名     : find_best_5_methods.R
# 分类       : 04_机器学习筛选特征基因
# 项目来源   : 5种机器学习_交集选择
# 用途       : 遍历/比较 12 种方法组合，选择合适的 5 种方法并取交集。
# 结果图     : Venn图；UpSet图；特征重要性图
# 主要 R 包  : dplyr; ggvenn; RColorBrewer; tidyverse; UpSetR
# 整理时间   : 2026-05-10
# =============================================================================
# 从12种机器学习方法中选择方法并绘制Venn图

library(ggvenn)
library(RColorBrewer)
library(dplyr)
library(tidyverse)
library(UpSetR)

# 设置工作目录
setwd("H:\\常用分析生信\\344.12种学习方法选择合适的五种取交集")

# ============================================================================
# ★★★ 用户自定义设置区域 ★★★
# 可用方法: "Lasso", "RF", "SVM", "LDA", "GBM", "ElasticNet",
#          "NeuralNet", "PLS", "AdaBoost", "Logistic", "NaiveBayes", "C5.0"

custom_methods <- c("AdaBoost", "ElasticNet", "NeuralNet", "PLS", "RF")

# ============================================================================
# 1. 读取所有12种方法的基因列表
# ============================================================================

csv_files <- list.files(pattern = "importanceGene\\..*\\.csv$", full.names = TRUE)
cat("找到", length(csv_files), "个基因文件:\n")
print(basename(csv_files))

gene_lists <- list()
method_names <- c()

for (file in csv_files) {
  method_name <- gsub(".*importanceGene\\.(.*)\\.csv", "\\1", basename(file))
  method_names <- c(method_names, method_name)
  data <- read.csv(file, stringsAsFactors = FALSE)
  gene_lists[[method_name]] <- data$variable
  cat(method_name, ": ", length(data$variable), "个基因\n")
}

cat("\n共有", length(gene_lists), "种方法\n")

# ============================================================================
# 2. 计算所有5种方法组合的交集基因数量并保存
# ============================================================================

all_combinations <- combn(method_names, 5, simplify = FALSE)
cat("\n共有", length(all_combinations), "种5方法组合\n")

results <- data.frame(
  combination = character(length(all_combinations)),
  intersection_count = integer(length(all_combinations)),
  intersection_genes = character(length(all_combinations)),
  custom_code = character(length(all_combinations)),
  stringsAsFactors = FALSE
)

cat("\n正在计算所有组合的交集基因数量...\n")

for (i in seq_along(all_combinations)) {
  combo <- all_combinations[[i]]
  intersection_genes <- Reduce(intersect, gene_lists[combo])
  results$combination[i] <- paste(combo, collapse = " + ")
  results$intersection_count[i] <- length(intersection_genes)
  # 交集基因名称
  results$intersection_genes[i] <- paste(intersection_genes, collapse = ", ")
  # 生成用户自定义代码格式
  results$custom_code[i] <- paste0('c("', paste(combo, collapse = '", "'), '")')
}

results <- results[order(-results$intersection_count), ]

cat("\n========================================\n")
cat("交集基因数量最多的前10种组合:\n")
cat("========================================\n")
print(head(results, 10))

# 保存所有组合结果供用户参考
write.csv(results, "all_combinations_results.csv", row.names = FALSE)
cat("\n所有组合结果已保存到: all_combinations_results.csv\n")

# ============================================================================
# 3. 使用用户自定义的方法
# ============================================================================

# 验证方法名称
invalid_methods <- custom_methods[!custom_methods %in% method_names]
if (length(invalid_methods) > 0) {
  stop("以下方法名称无效: ", paste(invalid_methods, collapse = ", "),
       "\n可用的方法名称: ", paste(method_names, collapse = ", "))
}

best_combination <- custom_methods
cat("\n========================================\n")
cat("使用自定义的", length(best_combination), "种方法:\n")
cat("========================================\n")
for (i in seq_along(best_combination)) {
  cat(i, ". ", best_combination[i], "\n", sep = "")
}

# 获取交集基因
best_gene_lists <- gene_lists[best_combination]
intersection_genes <- Reduce(intersect, best_gene_lists)

cat("\n交集基因数量:", length(intersection_genes), "\n")
cat("\n交集基因列表:\n")
print(intersection_genes)

# 保存交集基因
write.csv(data.frame(Gene = intersection_genes),
          "intersection_genes_best5.csv",
          row.names = FALSE)
cat("\n交集基因已保存到: intersection_genes_best5.csv\n")

# ============================================================================
# 4. 绘制Venn图
# ============================================================================

cat("\n正在绘制Venn图...\n")

venn_data <- best_gene_lists
names(venn_data) <- best_combination

p1 <- ggvenn(venn_data,
             fill_color = brewer.pal(5, "Set2"),
             stroke_size = 0.5,
             set_name_size = 4,
             text_size = 3,
             show_percentage = FALSE) +
  ggtitle("Machine Learning Methods - Gene Intersection") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

ggsave("Venn_best5_methods_ggvenn.pdf", p1, width = 10, height = 10)
cat("Venn图已保存: Venn_best5_methods_ggvenn.pdf\n")

# ============================================================================
# 5. 输出详细统计信息
# ============================================================================

cat("\n========================================\n")
cat("详细统计信息:\n")
cat("========================================\n")

cat("\n各方法基因数量:\n")
for (method in best_combination) {
  cat("  ", method, ": ", length(gene_lists[[method]]), "个基因\n", sep = "")
}

cat("\n两两方法交集:\n")
for (i in 1:(length(best_combination)-1)) {
  for (j in (i+1):length(best_combination)) {
    m1 <- best_combination[i]
    m2 <- best_combination[j]
    inter <- length(intersect(gene_lists[[m1]], gene_lists[[m2]]))
    cat("  ", m1, " ∩ ", m2, ": ", inter, "个基因\n", sep = "")
  }
}

# ============================================================================
# 6. 创建UpSet图
# ============================================================================

all_genes <- unique(unlist(best_gene_lists))
upset_data <- data.frame(Gene = all_genes)

for (method in best_combination) {
  upset_data[[method]] <- as.integer(all_genes %in% gene_lists[[method]])
}

pdf("UpSet_best5_methods.pdf", width = 12, height = 8, onefile = FALSE)
print(upset(upset_data[, -1],
            sets = best_combination,
            order.by = "freq",
            decreasing = TRUE,
            mb.ratio = c(0.6, 0.4),
            number.angles = 0,
            text.scale = 1.2,
            point.size = 3,
            line.size = 1,
            mainbar.y.label = "Intersection Size",
            sets.x.label = "Genes per Method",
            main.bar.color = "steelblue",
            sets.bar.color = "darkgreen"))
dev.off()

cat("\nUpSet图已保存: UpSet_best5_methods.pdf\n")
