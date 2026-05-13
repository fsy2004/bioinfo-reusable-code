# ==========================================================================
# 脚本名     : GEO疾病诊断模型分析.R
# 分类       : 05_诊断模型与验证
# 项目来源   : 从压缩包 493.模型分析GEO.rar 整理
# 原始文件   : 493.模型分析GEO\GEO疾病的诊断模型分析.R
# 用途       : 基于固定表达矩阵和 Genes.csv/gene.txt 构建 GEO 疾病诊断模型，输出联合模型和单基因诊断验证。
# 结果图     : 联合模型ROC；单基因ROC；表达箱线图；列线图；校准曲线；DCA决策曲线；OR森林图；预测概率表
# 非肿瘤消化适配: 很适合。非肿瘤消化系统GEO诊断模型可以直接作为验证模块。
# 主要 R 包  : rms; rmda; pROC; tools; ggplot2; ggpubr
# 整理日期   : 2026-05-13
# 备注       : 保留原始代码逻辑，仅添加统一说明头；运行前请把 workDir/setwd 和输入文件名改成当前项目路径。
# ==========================================================================
# ==== 诊断模型分析 pipeline：固定表达矩阵 + Genes.csv 基因列表 ====

suppressPackageStartupMessages({
  library(rms)
  library(rmda)
  library(pROC)
  library(tools)
  library(ggplot2)
  library(ggpubr)
})

# Step 0: 设置分析主目录
working_directory <- "H:\\常用分析生信\\493.模型分析GEO"
if (!dir.exists(working_directory)) {
  stop(paste("指定的工作目录不存在：", working_directory))
}
setwd(working_directory)

result_root <- file.path(working_directory, "Analysis_Result")
if (!dir.exists(result_root)) {
  dir.create(result_root, recursive = TRUE)
}

message("==== 疾病诊断模型分析启动 ====")

read_gene_list_txt <- function(file_path) {
  if (!file.exists(file_path)) {
    stop(sprintf("未找到基因文件: %s", file_path))
  }
  genes <- readLines(file_path, warn = FALSE, encoding = "UTF-8")
  genes <- trimws(genes)
  genes <- genes[genes != ""]
  genes <- unique(genes)
  if (length(genes) < 2) {
    stop("gene.txt 中有效基因数过少")
  }
  genes
}

read_gene_list_csv <- function(file_path) {
  if (!file.exists(file_path)) {
    stop(sprintf("未找到基因文件: %s", file_path))
  }
  gene_df <- tryCatch(
    read.csv(
      file_path,
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      fileEncoding = "UTF-8"
    ),
    error = function(e) {
      stop(sprintf("基因文件读取失败: %s", basename(file_path)))
    }
  )
  if (ncol(gene_df) == 0) {
    stop("基因文件为空")
  }
  gene_col <- if ("Gene" %in% colnames(gene_df)) "Gene" else colnames(gene_df)[1]
  genes <- trimws(as.character(gene_df[[gene_col]]))
  genes <- genes[!is.na(genes) & genes != ""]
  genes <- unique(genes)
  if (length(genes) < 2) {
    stop("Genes.csv 中有效基因数过少")
  }
  genes
}

read_expression_matrix <- function(file_path) {
  tryCatch(
    read.table(
      file_path,
      sep = ",",
      header = TRUE,
      row.names = 1,
      check.names = FALSE,
      fileEncoding = "UTF-8"
    ),
    error = function(e) {
      stop(sprintf("表达文件读取失败: %s", basename(file_path)))
    }
  )
}

extract_group_from_samples <- function(sample_names) {
  groups <- gsub("(.*)_([A-Za-z0-9]+)$", "\\2", sample_names)
  if (length(unique(groups)) < 2) {
    stop("样本名中无法识别至少两个分组，请检查样本名后缀格式")
  }
  groups
}

safe_auc <- function(labels, predictor) {
  roc_obj <- roc(labels, predictor, levels = c(0, 1), direction = "auto", quiet = TRUE)
  auc_val <- as.numeric(auc(roc_obj))
  if (auc_val < 0.5) {
    flipped_direction <- ifelse(roc_obj$direction == ">", "<", ">")
    roc_obj <- roc(labels, predictor, levels = c(0, 1), direction = flipped_direction, quiet = TRUE)
    auc_val <- as.numeric(auc(roc_obj))
  }
  list(roc = roc_obj, auc = auc_val)
}

run_single_dataset <- function(input_expr, feature_genes, result_dir) {
  if (!dir.exists(result_dir)) {
    dir.create(result_dir, recursive = TRUE)
  }

  dat_expr <- read_expression_matrix(input_expr)
  rownames(dat_expr) <- gsub("-", "_", rownames(dat_expr))

  found_genes <- feature_genes %in% rownames(dat_expr)
  if (sum(found_genes) == 0) {
    stop("Genes.csv 中基因在表达矩阵中均找不到")
  }
  if (any(!found_genes)) {
    warning(paste0("以下基因不在表达矩阵中：", paste(feature_genes[!found_genes], collapse = ",")))
  }

  dat_expr_filt <- dat_expr[feature_genes[found_genes], , drop = FALSE]
  df_expr <- as.data.frame(t(dat_expr_filt))
  sample_names <- rownames(df_expr)
  groups <- extract_group_from_samples(sample_names)
  df_expr$GroupType <- groups

  dd_name <- "ddinfo_current"
  assign(dd_name, datadist(df_expr), envir = .GlobalEnv)
  old_opt <- options(datadist = dd_name)
  on.exit({
    options(old_opt)
    if (exists(dd_name, envir = .GlobalEnv, inherits = FALSE)) {
      rm(list = dd_name, envir = .GlobalEnv)
    }
  }, add = TRUE)

  model_vars <- setdiff(colnames(df_expr), "GroupType")
  if (length(model_vars) < 2) {
    stop("可用于建模的基因数不足")
  }

  reg_formula <- as.formula(paste("GroupType ~", paste(model_vars, collapse = " + ")))
  lrm_fit <- lrm(reg_formula, data = df_expr, x = TRUE, y = TRUE)

  nomo_obj <- nomogram(
    lrm_fit,
    fun = plogis,
    fun.at = c(0.001, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99),
    lp = FALSE,
    funlabel = "Disease Risk"
  )
  pdf(file.path(result_dir, "Nomogram_Plot.pdf"), width = 11, height = 6)
  plot(nomo_obj)
  dev.off()

  calibrate_obj <- calibrate(lrm_fit, method = "boot", B = 800)
  pdf(file.path(result_dir, "Calibration_Curve.pdf"), width = 5.5, height = 5.5)
  cal_data <- calibrate_obj[, c("predy", "calibrated.orig", "calibrated.corrected")]
  plot(
    cal_data[, "predy"],
    cal_data[, "calibrated.corrected"],
    type = "l",
    lwd = 2,
    col = "black",
    xlab = "Predicted",
    ylab = "Observed",
    xlim = c(0, 1),
    ylim = c(0, 1)
  )
  abline(a = 0, b = 1, lty = 2, col = "gray50")
  lines(cal_data[, "predy"], cal_data[, "calibrated.orig"], lty = 1, lwd = 1, col = "darkblue")
  legend(
    "bottomright",
    legend = c("Ideal", "Apparent", "Bias-corrected"),
    lty = c(2, 1, 1),
    lwd = c(1, 1, 2),
    col = c("gray50", "darkblue", "black"),
    bty = "n"
  )
  dev.off()

  df_expr$GroupType <- as.numeric(factor(df_expr$GroupType)) - 1
  if (!all(df_expr$GroupType %in% c(0, 1))) {
    stop("GroupType 不是严格的 0/1")
  }

  pred_probs <- predict(lrm_fit, newdata = df_expr, type = "fitted")
  df_expr_dca <- df_expr
  df_expr_dca$Diagnostic_Model <- pred_probs

  dca_combined <- decision_curve(
    GroupType ~ Diagnostic_Model,
    data = df_expr_dca,
    thresholds = seq(0, 1, by = 0.01),
    family = binomial(link = "logit"),
    bootstraps = 100
  )

  dca_gene_list <- list()
  for (g in model_vars) {
    gene_formula <- as.formula(paste("GroupType ~", g))
    dca_gene_list[[g]] <- tryCatch(
      decision_curve(
        gene_formula,
        data = df_expr_dca,
        thresholds = seq(0, 1, by = 0.01),
        family = binomial(link = "logit"),
        bootstraps = 100
      ),
      error = function(e) NULL
    )
  }
  dca_gene_list <- dca_gene_list[!vapply(dca_gene_list, is.null, logical(1))]
  all_dca <- c(list(dca_combined), dca_gene_list)

  dca_colors <- c("red", "deepskyblue", "forestgreen", "orange", "purple", "magenta", "gold", "brown", "gray40", "cyan")
  while (length(dca_colors) < length(all_dca)) {
    dca_colors <- c(dca_colors, rainbow(length(all_dca) - length(dca_colors)))
  }
  dca_colors <- dca_colors[1:length(all_dca)]
  dca_labels <- c("Diagnostic_Model", names(dca_gene_list))

  pdf(file.path(result_dir, "DCA.pdf"), width = 9, height = 8)
  layout(matrix(c(1, 2), nrow = 2), heights = c(0.8, 5))
  par(mar = c(0, 4, 0.5, 2))
  plot.new()
  legend(
    "center",
    legend = c(dca_labels, "Treat All", "Treat None"),
    col = c(dca_colors, "black", "black"),
    lty = c(rep(1, length(dca_labels)), 1, 2),
    lwd = 2,
    ncol = 3,
    bty = "n",
    cex = 1.15
  )
  par(mar = c(5, 4, 0.5, 2))
  plot_decision_curve(
    all_dca,
    xlab = "Threshold Probability",
    col = dca_colors,
    confidence.intervals = FALSE,
    standardize = TRUE,
    cost.benefit.axis = TRUE,
    legend.position = "none"
  )
  dev.off()

  write.csv(df_expr, file = file.path(result_dir, "Filtered_Expression_Matrix.csv"), quote = FALSE)
  write.csv(as.data.frame(dca_combined$derived.data), file = file.path(result_dir, "Decision_Curve_Data.csv"))
  write.csv(
    data.frame(Sample = rownames(df_expr), GroupType = df_expr$GroupType, Predicted_Prob = pred_probs),
    file = file.path(result_dir, "Sample_Predicted_Probabilities.csv"),
    row.names = FALSE
  )

  coef_df <- as.data.frame(summary(lrm_fit))
  if ("Effect" %in% names(coef_df) && "S.E." %in% names(coef_df)) {
    coef_df$Effect <- as.numeric(as.character(coef_df$Effect))
    coef_df$`S.E.` <- as.numeric(as.character(coef_df$`S.E.`))
    coef_df$OR <- exp(coef_df$Effect)
    coef_df$OR_low <- exp(coef_df$Effect - 1.96 * coef_df$`S.E.`)
    coef_df$OR_high <- exp(coef_df$Effect + 1.96 * coef_df$`S.E.`)
    write.csv(coef_df, file = file.path(result_dir, "Model_Coefficients.csv"), row.names = FALSE)
    if ("P" %in% names(coef_df)) {
      write.csv(coef_df[!is.na(coef_df$P) & coef_df$P < 0.05, ], file = file.path(result_dir, "Model_Significant_Coefficients.csv"), row.names = FALSE)
    }
  }

  roc_y <- as.numeric(df_expr$GroupType)
  model_roc <- safe_auc(roc_y, pred_probs)
  model_roc_df <- data.frame(
    Threshold = model_roc$roc$thresholds,
    Specificity = model_roc$roc$specificities,
    Sensitivity = model_roc$roc$sensitivities,
    FPR = 1 - model_roc$roc$specificities
  )
  write.csv(model_roc_df, file = file.path(result_dir, "Diagnostic_Model_ROC_Data.csv"), row.names = FALSE)

  pdf(file.path(result_dir, "Diagnostic_Model_ROC.pdf"), width = 7, height = 7)
  par(mar = c(5, 5, 4, 2) + 0.1, cex = 1.2)
  plot(
    1 - model_roc$roc$specificities,
    model_roc$roc$sensitivities,
    type = "l",
    xlim = c(0, 1),
    ylim = c(0, 1),
    xlab = "1 - Specificity",
    ylab = "Sensitivity",
    main = "Diagnostic Model ROC Curve",
    col = "red",
    lwd = 3
  )
  abline(0, 1, lty = 2, col = "gray70", lwd = 2)
  legend(
    "bottomright",
    legend = sprintf("Diagnostic Model  AUC=%.3f", model_roc$auc),
    col = "red",
    lwd = 3,
    lty = 1,
    bty = "n"
  )
  dev.off()

  my_cols <- c("red", "deepskyblue", "forestgreen", "orange", "purple", "gray40", "black", "magenta", "gold", "brown")
  while (length(my_cols) < length(model_vars) + 1) {
    my_cols <- c(my_cols, rainbow(length(model_vars) + 1 - length(my_cols)))
  }

  roc_list <- list()
  auc_list <- c()
  leglab <- c()
  auc_dir <- c()

  pdf(file.path(result_dir, "IndividualGenes_ROC_SCI.pdf"), width = 8, height = 8)
  par(mar = c(5, 6, 4, 2) + 0.1, cex = 1.3)
  plot(
    0, 0,
    type = "n",
    xlim = c(0, 1),
    ylim = c(0, 1),
    xlab = expression("1 - Specificity"),
    ylab = "Sensitivity",
    main = "Gene ROC Curves",
    cex.lab = 1.4,
    cex.axis = 1.15,
    cex.main = 1.45
  )
  abline(0, 1, lty = 2, col = "gray70", lwd = 2)

  for (i in seq_along(model_vars)) {
    g <- model_vars[i]
    cur_roc <- safe_auc(roc_y, df_expr[[g]])
    lines(1 - cur_roc$roc$specificities, cur_roc$roc$sensitivities, col = my_cols[i + 1], lwd = 3)
    roc_list[[g]] <- cur_roc$roc
    auc_list <- c(auc_list, cur_roc$auc)
    mean0 <- mean(df_expr[[g]][roc_y == 0], na.rm = TRUE)
    mean1 <- mean(df_expr[[g]][roc_y == 1], na.rm = TRUE)
    direction_note <- ifelse(mean1 > mean0, "(组1高)", "(组0高)")
    auc_dir <- c(auc_dir, direction_note)
    leglab <- c(leglab, sprintf("%s  AUC=%.3f", g, cur_roc$auc))
  }
  legend("bottomright", legend = leglab, col = my_cols[2:(length(model_vars) + 1)], lwd = 3, lty = 1, bty = "n", cex = 1.15, y.intersp = 1.18)
  dev.off()

  write.csv(
    data.frame(
      Marker = c("Diagnostic_Model", model_vars),
      AUC = c(model_roc$auc, auc_list),
      Direction = c("Predicted_Prob", auc_dir),
      stringsAsFactors = FALSE
    ),
    file = file.path(result_dir, "ROC_AUC_Summary.csv"),
    row.names = FALSE
  )

  df_boxplot <- do.call(rbind, lapply(model_vars, function(gene) {
    data.frame(
      Gene = gene,
      Expression = df_expr[[gene]],
      Group = factor(df_expr$GroupType, labels = c("Control", "Disease"))
    )
  }))

  p_values <- sapply(model_vars, function(gene) {
    group0 <- df_expr[[gene]][df_expr$GroupType == 0]
    group1 <- df_expr[[gene]][df_expr$GroupType == 1]
    t.test(group0, group1)$p.value
  })
  p_labels <- sapply(p_values, function(p) {
    if (p < 0.001) {
      "***"
    } else if (p < 0.01) {
      "**"
    } else if (p < 0.05) {
      "*"
    } else {
      "ns"
    }
  })

  y_max <- max(df_boxplot$Expression, na.rm = TRUE)
  y_min <- min(df_boxplot$Expression, na.rm = TRUE)
  y_range <- y_max - y_min

  pdf(file.path(result_dir, "Diagnostic_Genes_Boxplot.pdf"), width = max(3, length(model_vars) * 0.8), height = 6)
  p <- ggplot(df_boxplot, aes(x = Gene, y = Expression, fill = Group)) +
    geom_boxplot(
      width = 0.6,
      outlier.shape = 21,
      outlier.size = 2,
      outlier.fill = "white",
      outlier.stroke = 1,
      alpha = 0.8,
      linewidth = 0.8
    ) +
    scale_fill_manual(values = c("Control" = "#E8E8E8", "Disease" = "#FF6B6B"), name = "Group") +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 13, face = "bold"),
      axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 11),
      legend.position = "top",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11),
      plot.title = element_text(size = 8, face = "bold", hjust = 0.5)
    ) +
    labs(y = "Expression Level", title = "Diagnostic Model Genes Expression Comparison")

  y_pos <- y_max + y_range * 0.05
  for (i in seq_along(model_vars)) {
    p <- p + annotate("text", x = i, y = y_pos, label = p_labels[i], size = 5, fontface = "bold")
  }
  print(p)
  dev.off()

  write.csv(
    data.frame(Gene = model_vars, P_Value = p_values, Significance = p_labels),
    file = file.path(result_dir, "Diagnostic_Genes_PValues.csv"),
    row.names = FALSE
  )
}

input_expr <- file.path(working_directory, "Sample Type Matrix.csv")
gene_file <- file.path(working_directory, "Genes.csv")

if (!file.exists(input_expr)) {
  stop(sprintf("未找到表达矩阵文件: %s", input_expr))
}

feature_genes <- read_gene_list_csv(gene_file)
message(sprintf("读取到 %d 个研究基因", length(feature_genes)))

analysis_summary <- list()
dataset_name <- file_path_sans_ext(basename(input_expr))
result_dir <- file.path(result_root, dataset_name)

message("========================================")
message(sprintf("开始分析: %s", basename(input_expr)))
message(sprintf("使用基因文件: %s", basename(gene_file)))
message(sprintf("输出目录: %s", result_dir))

status <- tryCatch({
  run_single_dataset(input_expr, feature_genes, result_dir)
  "OK"
}, error = function(e) {
  message(sprintf("分析失败: %s", e$message))
  paste("FAILED:", e$message)
})

analysis_summary[[length(analysis_summary) + 1]] <- data.frame(
  Dataset = basename(input_expr),
  GeneFile = basename(gene_file),
  ResultDir = result_dir,
  Status = status,
  stringsAsFactors = FALSE
)

summary_df <- do.call(rbind, analysis_summary)
write.csv(summary_df, file = file.path(result_root, "Analysis_Summary.csv"), row.names = FALSE)

message("==== 全流程执行完毕！====")
message(sprintf("汇总结果: %s", file.path(result_root, "Analysis_Summary.csv")))
