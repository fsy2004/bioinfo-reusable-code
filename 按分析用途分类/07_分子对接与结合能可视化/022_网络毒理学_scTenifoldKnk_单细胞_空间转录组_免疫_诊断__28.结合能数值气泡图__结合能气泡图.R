# =============================================================================
# 编号       : R022
# 脚本名     : 结合能气泡图.R
# 分类       : 07_分子对接与结合能可视化
# 项目来源   : 网络毒理学_scTenifoldKnk_单细胞_空间转录组_免疫_诊断
# 用途       : 整理分子对接结合能数据并绘制结合能气泡图。
# 结果图     : 气泡图/点图
# 主要 R 包  : dplyr; ggplot2; readxl; scales
# 整理时间   : 2026-05-10
# =============================================================================
#===============================================================================
# 单个化合物与基因结合能气泡图
#===============================================================================

# 加载所需R包
library(readxl)       # 读取Excel文件
library(ggplot2)      # 绘图
library(dplyr)        # 数据处理
library(scales)       # 数据缩放

# 设置工作目录
setwd("C:/Users/fsy/Desktop/网络毒理学+单细胞+scTenifoldKnk模拟基因敲除+空间转录组+细胞轨迹+机器学习+分子对接+免疫分析+诊断模型/28.结合能数值气泡图")

# 读取数据
data <- read_excel("结合能.xlsx")

# 查看数据结构
cat("Data dimensions:", dim(data), "\n")
cat("Column names:", colnames(data), "\n")
print(head(data))

# 整理数据格式
# 第一列是基因名，第二列是结合能值
gene_names <- data[[1]]
compound_name <- colnames(data)[2]
binding_values <- as.numeric(data[[2]])

# 创建数据框
bubble_data <- data.frame(
  Gene = gene_names,
  Compound = compound_name,
  BindingEnergy = binding_values,
  stringsAsFactors = FALSE
)

# 按结合能排序
bubble_data <- bubble_data %>% arrange(BindingEnergy)
bubble_data$Gene <- factor(bubble_data$Gene, levels = bubble_data$Gene)

cat("\n整理后的数据:\n")
print(bubble_data)

#===============================================================================
# 绘制高级气泡图
#===============================================================================

# 确定颜色（负值为蓝色表示强结合，正值为红色表示弱结合）
bubble_data$Color <- ifelse(bubble_data$BindingEnergy < -5, "#1f77b4",
                            ifelse(bubble_data$BindingEnergy < -3, "#4da6ff",
                                   ifelse(bubble_data$BindingEnergy < 0, "#99ccff",
                                          ifelse(bubble_data$BindingEnergy < 3, "#ffcc99", "#ff6666"))))

# 创建气泡图
p <- ggplot(bubble_data, aes(x = BindingEnergy, y = Gene)) +
  # 添加背景网格
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8, alpha = 0.5) +

  # 绘制气泡
  geom_point(aes(fill = BindingEnergy),
             shape = 21, color = "black", stroke = 1.2, alpha = 0.8, size = 10) +

  # 添加数值标签
  geom_text(aes(label = sprintf("%.2f", BindingEnergy)),
            size = 5, fontface = "bold", color = "black", hjust = 0.5, vjust = 0.5) +

  # 设置颜色渐变（蓝色→白色→红色）
  scale_fill_gradientn(
    colors = c("#1f77b4", "#4da6ff", "#99ccff", "white", "#ffcc99", "#ff6666", "#cc0000"),
    name = "Binding Energy\n(kcal/mol)",
    guide = guide_colorbar(
      barwidth = 1.2,
      barheight = 12,
      title.position = "top",
      title.hjust = 0.5
    )
  ) +

  # 设置气泡大小
  scale_size_continuous(range = c(3, 12), guide = "none") +

  # 标签和标题
  labs(
    title = sprintf("Binding Energy: %s vs Genes", compound_name),
    x = "Binding Energy (kcal/mol)",
    y = "Gene"
  ) +

  # 主题设置
  theme_minimal(base_size = 12) +
  theme(
    # 标题设置
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15, color = "#333333"),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "#666666", margin = margin(b = 10)),

    # 坐标轴设置
    axis.text.x = element_text(face = "bold", color = "black", size = 11),
    axis.text.y = element_text(face = "bold", color = "black", size = 11),
    axis.title = element_text(face = "bold", size = 12, color = "#333333"),
    axis.line = element_line(color = "black", linewidth = 0.8),

    # 图例设置
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11, color = "#333333"),
    legend.text = element_text(size = 10, color = "#555555"),
    legend.background = element_rect(fill = "white", color = "gray80", linewidth = 0.5),
    legend.margin = margin(10, 10, 10, 10),

    # 面板设置
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),

    # 边距设置
    plot.margin = margin(15, 15, 15, 15)
  )

# 保存气泡图
ggsave("Binding_Energy_Bubble_Plot.pdf", plot = p, width = 10, height = 8, dpi = 300)

cat("\n气泡图已保存为: Binding_Energy_Bubble_Plot.pdf\n")

#===============================================================================
# 数据统计分析
#===============================================================================

cat("\n=== 结合能统计分析 ===\n")

# 基本统计
cat("\n化合物:", compound_name, "\n")
cat("基因数量:", nrow(bubble_data), "\n")
cat("平均结合能:", sprintf("%.2f", mean(bubble_data$BindingEnergy)), "kcal/mol\n")
cat("中位数结合能:", sprintf("%.2f", median(bubble_data$BindingEnergy)), "kcal/mol\n")
cat("标准差:", sprintf("%.2f", sd(bubble_data$BindingEnergy)), "kcal/mol\n")

# 最强结合（最负值）
strongest_idx <- which.min(bubble_data$BindingEnergy)
cat("\n最强结合:\n")
cat("  基因:", bubble_data$Gene[strongest_idx], "\n")
cat("  结合能:", sprintf("%.2f", bubble_data$BindingEnergy[strongest_idx]), "kcal/mol\n")

# 最弱结合（最正值）
weakest_idx <- which.max(bubble_data$BindingEnergy)
cat("\n最弱结合:\n")
cat("  基因:", bubble_data$Gene[weakest_idx], "\n")
cat("  结合能:", sprintf("%.2f", bubble_data$BindingEnergy[weakest_idx]), "kcal/mol\n")

# 强结合基因（结合能 < -5）
strong_binding <- bubble_data %>% filter(BindingEnergy < -5)
cat("\n强结合基因 (结合能 < -5 kcal/mol):\n")
if(nrow(strong_binding) > 0) {
  for(i in 1:nrow(strong_binding)) {
    cat("  ", strong_binding$Gene[i], ": ", sprintf("%.2f", strong_binding$BindingEnergy[i]), " kcal/mol\n", sep="")
  }
} else {
  cat("  无\n")
}

cat("\n分析完成！\n")
cat("\n生成的文件:\n")
cat("  Binding_Energy_Bubble_Plot.pdf\n")

