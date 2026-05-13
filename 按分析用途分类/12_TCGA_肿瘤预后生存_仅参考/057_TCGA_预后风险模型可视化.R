# ==========================================================================
# 脚本名     : TCGA预后风险模型可视化.R
# 分类       : 12_TCGA_肿瘤预后生存_仅参考
# 项目来源   : 从压缩包 414.TCGA预后诊断模型的分析.rar 整理
# 原始文件   : 414.TCGA预后诊断模型的分析\模型分析-生存曲线roc基因分布生存状态图风险分分布.R
# 用途       : 基于 TCGA 风险评分文件绘制风险分布、生存状态、风险基因热图、KM生存曲线和1/3/5年时间依赖ROC。
# 结果图     : riskScore风险分布图；survStat生存状态散点图；风险基因表达热图；KM生存曲线；1/3/5年timeROC曲线；AUC表
# 非肿瘤消化适配: 肿瘤参考。非肿瘤消化系统若有队列随访和风险评分可改造，否则主要借鉴图形样式。
# 主要 R 包  : survival; survminer; timeROC; ggplot2; viridis; ggsci; MetBrewer; ComplexHeatmap; circlize
# 整理日期   : 2026-05-13
# 备注       : 保留原始代码逻辑，仅添加统一说明头；运行前请把 workDir/setwd 和输入文件名改成当前项目路径。
# ==========================================================================
# 加载必要的R包
library(survival)        # 生存分析
library(survminer)       # 生存曲线可视化
library(timeROC)         # 时间依赖ROC
library(ggplot2)         # ggplot2绘图
library(viridis)         # 渐变色
library(ggsci)           # SCI色板
library(MetBrewer)       # MetBrewer色板
library(ComplexHeatmap)  # 高级热图
library(circlize)        # 热图色彩

# 设置工作目录
setwd("H:\\常用分析生信\\414.TCGA预后诊断模型的分析")

# ==========================================
#        1. 风险评分分布图
# ==========================================
risk_data <- read.table("risk.csv", sep=",", header=TRUE, row.names=1, check.names=FALSE)
risk_data$futime <- risk_data$futime / 365
risk_data <- risk_data[order(risk_data$riskScore), ]
n_samples <- nrow(risk_data)

# 计算低风险组和高风险组的样本数
risk_class <- risk_data[, "risk"]
n_low  <- sum(risk_class == "low")
n_high <- sum(risk_class == "high")

# 限制风险评分最大值为10
risk_data$riskScore[risk_data$riskScore > 10] <- 10
risk_data$patient_rank <- 1:n_samples
risk_data$risk_group <- factor(risk_data$risk, levels = c("low", "high"))

# 计算中位数
median_risk <- median(risk_data$riskScore)

# 绘制风险评分分布图
pdf(file="riskScore.pdf", width = 10, height = 3.5)
print(
  ggplot(risk_data, aes(x=patient_rank, y=riskScore, color=risk_group)) +
    geom_point(size=2.2, alpha=0.85) +
    scale_color_manual(values = c("low" = "#0072B5", "high" = "#BC3C29")) +
    geom_vline(xintercept=n_low, linetype="dashed", color="grey50", size=1) +
    geom_hline(yintercept=median_risk, linetype="dashed", color="grey50", size=1) +
    annotate("text",
             x = n_samples/2,
             y = median_risk,
             label = paste0("Median = ", round(median_risk, 2)),
             vjust = -1,
             hjust = 0.5,
             size = 5,
             color = "grey30") +
    theme_bw(base_size=16) +
    labs(x="Patients (increasing risk score)", y="Risk score") +
    theme(
      legend.position="top",
      axis.line = element_line(color="black", size=0.8),
      axis.text = element_text(color="black", size=14),
      axis.title = element_text(color="black", size=16, face="bold"),
      panel.grid = element_blank(),
      plot.background = element_rect(fill="white", color="white")
    ) +
    guides(color=guide_legend(title="Risk Group"))
)
dev.off()

# ==========================================
#        2. 生存状态散点图
# ==========================================
risk_data$survival_status <- factor(risk_data$fustat, levels = c(0,1), labels = c("Alive", "Dead"))
risk_data$risk_label <- factor(risk_data$risk, levels = c("low", "high"), labels = c("Low risk", "High risk"))
color_risk   <- c("Low risk" = "#0072B5", "High risk" = "#BC3C29")
shape_status <- c("Alive" = 16, "Dead" = 17)

pdf(file="survStat.pdf", width = 10, height = 3.5)
print(
  ggplot(risk_data, aes(x=patient_rank, y=futime, color=risk_label, shape=survival_status)) +
    geom_point(size=3, alpha=0.85) +
    scale_color_manual(values = color_risk) +
    scale_shape_manual(values = shape_status) +
    geom_vline(xintercept=sum(risk_data$risk_label=="Low risk"),
               linetype="dashed", color="grey50", size=1) +
    theme_bw(base_size=16) +
    labs(
      x="Patients (increasing risk score)",
      y="Survival time (years)",
      color="Risk group",
      shape="Status"
    ) +
    theme(
      legend.position="top",
      legend.title=element_text(size=14, face="bold"),
      legend.text=element_text(size=12),
      axis.line = element_line(color="black", size=0.8),
      axis.text = element_text(color="black", size=14),
      axis.title = element_text(color="black", size=16, face="bold"),
      panel.grid = element_blank(),
      plot.background = element_rect(fill="white", color="white")
    ) +
    guides(
      color=guide_legend(override.aes = list(size=4)),
      shape=guide_legend(override.aes = list(size=4))
    )
)
dev.off()

# ==========================================
#        3. 基因表达热图
# ==========================================
numeric_cols <- sapply(risk_data, is.numeric)
non_expr_cols <- c("riskScore", "futime", "fustat", "patient_rank")
numeric_cols[names(numeric_cols) %in% non_expr_cols] <- FALSE
expr_matrix <- log2(risk_data[, numeric_cols] + 0.01)
expr_matrix <- t(expr_matrix)

heatmap_anno <- data.frame(type = risk_data$risk)
rownames(heatmap_anno) <- rownames(risk_data)
top_anno <- HeatmapAnnotation(df = heatmap_anno, col = list(type = c("low"="#0072B5", "high"="#BC3C29")))

pdf(file="heatmap.pdf", width = 10, height = 3)
Heatmap(expr_matrix,
        name = "log2(Expression)",
        top_annotation = top_anno,
        col = colorRamp2(
          c(min(expr_matrix, na.rm=TRUE), median(expr_matrix, na.rm=TRUE), max(expr_matrix, na.rm=TRUE)),
          c("#0072B5", "white", "#BC3C29")
        ),
        show_row_names = TRUE,
        show_column_names = FALSE,
        cluster_columns = FALSE,
        row_names_gp = gpar(fontsize = 11),
        column_names_gp = gpar(fontsize = 3))
dev.off()

# ==========================================
#        4. Kaplan-Meier 生存分析
# ==========================================
surv_data <- read.table("risk.csv", header=TRUE, sep=",", check.names=FALSE)
surv_data$futime <- surv_data$futime / 365
surv_data$risk <- factor(surv_data$risk, levels = c("low", "high"))

km_fit  <- survfit(Surv(futime, fustat) ~ risk, data=surv_data)
cox_fit <- coxph(Surv(futime, fustat) ~ risk, data=surv_data)

cox_summary  <- summary(cox_fit)
hazard_ratio <- cox_summary$coefficients[1, 2]
hr_ci        <- cox_summary$conf.int[1, c("lower .95", "upper .95")]
hr_pvalue    <- cox_summary$coefficients[1, 5]

hr_text <- paste0("HR = ", sprintf("%.2f", hazard_ratio), " (95% CI: ",
                  sprintf("%.2f", hr_ci[1]), "-", sprintf("%.2f", hr_ci[2]), ")")
p_text  <- if (hr_pvalue < 0.001) "p < 0.001" else paste0("p = ", sprintf("%.3f", hr_pvalue))
annotation_text <- paste(hr_text, p_text, sep = "\n")

km_plot <- ggsurvplot(
  km_fit,
  data = surv_data,
  conf.int = TRUE,
  pval = annotation_text,
  pval.size = 5,
  risk.table = TRUE,
  legend.labs = c("Low risk", "High risk"),
  legend.title = "Risk",
  xlab = "Time (years)",
  break.time.by = 1,
  risk.table.title = "",
  palette = c("#0072B5", "#BC3C29"),
  risk.table.height = .25,
  ggtheme = theme_bw(base_size = 18)
)
pdf(file = "survival.pdf", onefile = FALSE, width = 9, height = 7)
print(km_plot)
dev.off()

# ==========================================
#        5. 时间依赖ROC曲线
# ==========================================
roc_data <- read.table("risk.csv", header=TRUE, sep=",")
roc_data$futime <- roc_data$futime / 365

roc_result <- timeROC(T=roc_data$futime, delta=roc_data$fustat,
                      marker=roc_data$riskScore, cause=1,
                      weighting='aalen',
                      times=c(1, 3, 5), ROC=TRUE)

roc_curve_df <- data.frame(
  FPR  = c(roc_result$FP[,1], roc_result$FP[,2], roc_result$FP[,3]),
  TPR  = c(roc_result$TP[,1], roc_result$TP[,2], roc_result$TP[,3]),
  Time = factor(rep(c("1 year", "3 years", "5 years"), each=nrow(roc_result$FP)))
)
auc_labels <- paste0("AUC at ", c(1,3,5), " years: ", round(roc_result$AUC, 3))

color_roc <- c("#0072B5", "#BC3C29", "#20854E")

pdf(file="ROC.pdf", width=6, height=6)
print(
  ggplot(roc_curve_df, aes(x=FPR, y=TPR, color=Time)) +
    geom_line(size=2) +
    scale_color_manual(values=color_roc) +
    geom_abline(intercept=0, slope=1, linetype="dashed", color="grey50", size=1) +
    theme_bw(base_size=16) +
    labs(x="False Positive Rate", y="True Positive Rate", color="Time") +
    annotate("text", x=0.3, y=seq(0.2,0.05,length.out=3),
             label=auc_labels, hjust=0, size=5,
             color=color_roc) +
    theme(
      legend.position="top",
      axis.line = element_line(color="black", size=0.8),
      axis.text = element_text(color="black", size=14),
      axis.title = element_text(color="black", size=16, face="bold"),
      panel.grid = element_blank(),
      plot.background = element_rect(fill="white", color="white")
    )
)
dev.off()

# ROC曲线AUC表
auc_table <- data.frame(
  Time = c(1, 3, 5),
  AUC  = roc_result$AUC
)
write.csv(auc_table, file="roc_auc.csv", row.names=FALSE)
