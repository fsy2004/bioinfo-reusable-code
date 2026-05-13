# ==========================================================================
# 脚本名     : SHAP机器学习解释分析.R
# 分类       : 04_机器学习筛选特征基因
# 项目来源   : 从压缩包 329.SHAP分析.rar 整理
# 原始文件   : 329.SHAP分析\shap.R
# 用途       : 训练多种机器学习分类模型，选择最优模型后使用 kernelshap/shapviz 解释特征贡献。
# 结果图     : 多模型ROC曲线；SHAP重要性条形图；SHAP蜂群图；SHAP dependence图；SHAP密度图；SHAP-表达散点图；SHAP热图；累计贡献曲线；分组SHAP条形图；permutation importance；waterfall图；force图
# 非肿瘤消化适配: 适合。可用于非肿瘤消化系统疾病诊断模型、GEO队列分类、特征基因解释。
# 主要 R 包  : caret; DALEX; ggplot2; randomForest; kernlab; kernelshap; pROC; shapviz; xgboost; klaR; pheatmap; gbm; patchwork
# 整理日期   : 2026-05-13
# 备注       : 保留原始代码逻辑，仅添加统一说明头；运行前请把 workDir/setwd 和输入文件名改成当前项目路径。
# ==========================================================================
# ===================== 包加载 =====================
library(caret)
library(DALEX)
library(ggplot2)
library(randomForest)
library(kernlab)
library(kernelshap)
library(pROC)
library(shapviz)
library(xgboost)
library(klaR)
library(RColorBrewer)
library(pheatmap)
library(gbm)  # 用于GBM模型
library(patchwork)  # 用于组合图表

# ===================== 参数设置 =====================
work_dir <- "H:\\常用分析生信\\329.SHAP分析"
result_dir <- "results"  # 结果输出文件夹

# 创建子文件夹
subdir_model <- "01_Model_Evaluation"
subdir_importance <- "02_Feature_Importance"
subdir_dependence <- "03_Dependence_Plots"
subdir_density <- "04_Density_Plots"
subdir_scatter <- "05_Scatter_Plots"
# subdir_interaction <- "06_Interaction_Plots"  # 已删除交互图生成
subdir_heatmap <- "07_Heatmap"
subdir_summary <- "08_Summary_Statistics"

random_seed <- 12345
num_colors <- 20  # pastel色系最多9色，超出自动渐变扩展
train_ratio <- 0.7
cv_folds <- 5

# 基于子文件夹的文件路径
roc_file <- file.path(result_dir, subdir_model, "01_model_roc_curve.pdf")
barplot_file <- file.path(result_dir, subdir_importance, "01_shap_importance_barplot.pdf")
bee_file <- file.path(result_dir, subdir_importance, "02_shap_importance_beeswarm.pdf")
dependence_file <- file.path(result_dir, subdir_dependence, "00_shap_dependence_all.pdf")
waterfall_file <- file.path(result_dir, subdir_summary, "01_shap_waterfall_sample1.pdf")
force_file <- file.path(result_dir, subdir_summary, "02_shap_force_sample1.pdf")
heatmap_file <- file.path(result_dir, subdir_heatmap, "01_shap_value_heatmap.pdf")
cumulative_contribution_file <- file.path(result_dir, subdir_summary, "03_cumulative_contribution_curve.pdf")
group_importance_file <- file.path(result_dir, subdir_importance, "03_shap_group_importance_barplot.pdf")
permutation_importance_file <- file.path(result_dir, subdir_importance, "04_feature_permutation_importance.pdf")
# interaction_scatter_file <- file.path(result_dir, subdir_interaction, "01_shap_interaction_scatter.pdf")  # 已删除
top20_file <- file.path(result_dir, subdir_summary, "shap_top_genes.csv")
feature_importance_csv <- file.path(result_dir, subdir_summary, "feature_importance.csv")
groupwise_importance_csv <- file.path(result_dir, subdir_summary, "groupwise_importance.csv")
auc_summary_csv <- file.path(result_dir, subdir_summary, "model_auc_summary.csv")
influential_samples_csv <- file.path(result_dir, subdir_summary, "influential_samples.csv")
permutation_importance_csv <- file.path(result_dir, subdir_summary, "permutation_importance.csv")

data_file <- "geneexp.csv"
pdf_width <- 6
pdf_height <- 6
progress_style <- 3

ml_methods <- data.frame(
  ModelName = c("RF", "SVM", "XGB", "GBM", "KNN"),
  MethodID = c("rf", "svmRadial", "xgbTree", "gbm", "knn")
)

# ===================== 步骤1：环境准备 =====================
cat("步骤1/8: 设置工作目录和随机种子...\n")
if (!dir.exists(work_dir)) stop("工作目录不存在，请检查路径！")
setwd(work_dir)
set.seed(random_seed)

# 创建结果输出文件夹及其子文件夹
if (!dir.exists(result_dir)) dir.create(result_dir)
if (!dir.exists(file.path(result_dir, subdir_model))) dir.create(file.path(result_dir, subdir_model), recursive=TRUE)
if (!dir.exists(file.path(result_dir, subdir_importance))) dir.create(file.path(result_dir, subdir_importance), recursive=TRUE)
if (!dir.exists(file.path(result_dir, subdir_dependence))) dir.create(file.path(result_dir, subdir_dependence), recursive=TRUE)
if (!dir.exists(file.path(result_dir, subdir_density))) dir.create(file.path(result_dir, subdir_density), recursive=TRUE)
if (!dir.exists(file.path(result_dir, subdir_scatter))) dir.create(file.path(result_dir, subdir_scatter), recursive=TRUE)
# if (!dir.exists(file.path(result_dir, subdir_interaction))) dir.create(file.path(result_dir, subdir_interaction), recursive=TRUE)  # 已删除交互图文件夹
if (!dir.exists(file.path(result_dir, subdir_heatmap))) dir.create(file.path(result_dir, subdir_heatmap), recursive=TRUE)
if (!dir.exists(file.path(result_dir, subdir_summary))) dir.create(file.path(result_dir, subdir_summary), recursive=TRUE)

# ===================== 步骤2：配色 =====================
cat("步骤2/8: 生成可爱风格配色方案...\n")
if (num_colors > 9) {
  color_palette <- colorRampPalette(brewer.pal(9, "Pastel1"))(num_colors)
} else {
  color_palette <- brewer.pal(num_colors, "Pastel1")
}

# ===================== 步骤3：数据读取与预处理 =====================
cat("步骤3/8: 读取并预处理表达数据...\n")
if (!file.exists(data_file)) stop("数据文件不存在！")
raw_data <- read.csv(data_file, header=TRUE, check.names=FALSE, row.names=1)
if (nrow(raw_data) == 0) stop("数据文件为空！")
transposed_data <- t(raw_data)
sample_ids <- rownames(transposed_data)
if (!all(grepl("(_con$|_tra$)", sample_ids))) stop("样本名需以_con或_tra结尾！")
group_labels <- ifelse(grepl("_con$", sample_ids), "Control",
                       ifelse(grepl("_tra$", sample_ids), "Treatment", NA))
if (any(is.na(group_labels))) stop("分组标签生成失败！")
expr_data <- as.data.frame(transposed_data)
expr_data$Group <- as.factor(group_labels)

# 数据验证 - 确保所有列都是数值型
for (col in colnames(expr_data)[-ncol(expr_data)]) {
  expr_data[[col]] <- as.numeric(expr_data[[col]])
}

# ===================== 步骤4：训练集/测试集划分 =====================
cat("步骤4/8: 划分训练集和测试集...\n")
if (length(unique(expr_data$Group)) != 2) stop("分组数不是2，无法二分类！")
split_idx <- createDataPartition(y=expr_data$Group, p=train_ratio, list=FALSE)
train_set <- expr_data[split_idx, ]
test_set <- expr_data[-split_idx, ]
test_labels <- test_set$Group
test_features <- test_set[, -ncol(test_set)]

# ===================== 步骤5：模型训练与评估 =====================
cat("步骤5/8: 多模型训练与ROC评估...\n")

# 诊断：检查caret和train_set
cat("DEBUG: 检查caret包...\n")
if (!("train" %in% ls("package:caret"))) {
  cat("ERROR: train函数未找到，重新加载caret...\n")
  detach("package:caret", unload=TRUE)
  library(caret)
}

cat("DEBUG: train_set数据结构\n")
cat("  行数:", nrow(train_set), "\n")
cat("  列数:", ncol(train_set), "\n")
cat("  Group列类:\n")
print(class(train_set$Group))
cat("  前5个基因列类:\n")
print(sapply(train_set[, 1:min(5, ncol(train_set)-1)], class))

auc_results <- c()
model_auc_map <- list()
roc_colors <- color_palette[1:nrow(ml_methods)]
roc_ci_list <- list()
pb <- txtProgressBar(min=0, max=nrow(ml_methods), style=progress_style)

pdf(file=roc_file, width=pdf_width, height=pdf_height)
par(mar=c(5, 5, 4, 2) + 0.1)
for (i in 1:nrow(ml_methods)) {
  setTxtProgressBar(pb, i)
  mdl_name <- ml_methods$ModelName[i]
  mdl_id <- ml_methods$MethodID[i]
  if (nrow(train_set) < 10) stop("训练集样本太少！")
  if (mdl_id == "svmRadial") {
    mdl <- caret::train(Group ~ ., data=train_set, method=mdl_id, prob.model=TRUE,
                 trControl=caret::trainControl(method="repeatedcv", number=cv_folds, savePredictions=TRUE))
  } else {
    mdl <- caret::train(Group ~ ., data=train_set, method=mdl_id,
                 trControl=caret::trainControl(method="repeatedcv", number=cv_folds, savePredictions=TRUE))
  }
  pred_prob <- predict(mdl, newdata=test_features, type="prob")
  if (!"Treatment" %in% colnames(pred_prob)) stop("预测概率列名不含Treatment！")
  roc_obj <- roc(as.numeric(test_labels)-1, as.numeric(pred_prob[,"Treatment"]))
  auc_val <- as.numeric(roc_obj$auc)
  auc_ci <- ci.auc(roc_obj, conf.level=0.95)
  roc_ci_list[[mdl_id]] <- auc_ci
  auc_results <- c(auc_results,
                   sprintf("%s: %.03f [%.03f, %.03f]", mdl_name, auc_val, auc_ci[1], auc_ci[3]))
  model_auc_map[[mdl_id]] <- auc_val
  if (i == 1) {
    plot(roc_obj, print.auc=FALSE, legacy.axes=TRUE, main="", col=roc_colors[i], lwd=3)
  } else {
    plot(roc_obj, print.auc=FALSE, legacy.axes=TRUE, main="", col=roc_colors[i], lwd=3, add=TRUE)
  }
}
close(pb)
legend('bottomright', auc_results, col=roc_colors, lwd=3, bty='n', cex=0.9)
dev.off()

# ===================== 步骤6：最优模型选择与SHAP计算 =====================
cat("步骤6/8: 选择最优模型并计算SHAP值...\n")
auc_vec <- unlist(model_auc_map)
if (length(auc_vec) == 0) stop("AUC结果为空！")
best_method <- names(which.max(auc_vec))
cat("AUC最高模型为：", best_method, "\n")
final_model <- caret::train(Group ~ ., data=train_set, method=best_method,
                     trControl=caret::trainControl(method="repeatedcv", number=cv_folds, savePredictions=TRUE))
shap_fit <- kernelshap(
  object = final_model,
  X = train_set[, -ncol(train_set)],
  pred_fun = function(obj, newdata) {
    pred <- predict(obj, newdata, type="prob")
    if (!"Treatment" %in% colnames(pred)) stop("SHAP预测概率列名不含Treatment！")
    pred[,"Treatment"]
  }
)

# ===================== 步骤7：SHAP可视化 =====================
cat("步骤7/8: SHAP可视化...\n")
shap_obj <- shapviz(shap_fit, X_pred = train_set[, -ncol(train_set)], X = train_set[, -ncol(train_set)], interactions=TRUE)
shap_importance <- sort(colMeans(abs(shap_obj$S)), decreasing=TRUE)
top_features <- names(shap_importance)

# 统一主题
custom_theme <- theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = "right")

options(shapviz.colors = color_palette)

# SHAP重要性条形图
pdf(file=barplot_file, width=pdf_width, height=pdf_height)
print(sv_importance(shap_obj, kind="bar", show_numbers=TRUE) + custom_theme)
dev.off()

# SHAP重要性蜂群图
pdf(file=bee_file, width=pdf_width, height=pdf_height)
print(sv_importance(shap_obj, kind="bee", show_numbers=TRUE) + custom_theme)
dev.off()

# 所有Top特征依赖图 - 删除了all.pdf的生成，改为单个图，最后会组合

# 针对每个Top基因，分别绘制依赖图、密度图、散点图
cat("为每个Top基因绘制单独的依赖图、密度图、散点图...\n")

# 创建dependence图列表，用于最后组合
dependence_plots_list <- list()

for (i in seq_along(top_features)) {
  gene <- top_features[i]
  safe_gene <- tolower(gsub("[^A-Za-z0-9_]", "_", gene))

  # 添加序号 (格式：01, 02, 03...)
  idx_str <- sprintf("%02d", i)

  # 依赖图 - 保存为ggplot对象
  dependence_plot <- sv_dependence(shap_obj, v = gene) +
    custom_theme +
    ggtitle(gene) +
    theme(plot.title = element_text(size = 11, face = "bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.title = element_blank())

  dependence_plots_list[[i]] <- dependence_plot

  # 同时保存单个pdf - 添加序号
  pdf(file = file.path(result_dir, subdir_dependence, paste0(idx_str, "_", safe_gene, ".pdf")), width = pdf_width, height = pdf_height)
  print(dependence_plot)
  dev.off()

  # 密度图
  shap_df <- data.frame(
    SHAP = shap_obj$S[, gene],
    Group = train_set$Group
  )
  pdf(file = file.path(result_dir, subdir_density, paste0(idx_str, "_", safe_gene, ".pdf")), width=pdf_width, height=pdf_height)
  print(ggplot(shap_df, aes(x=SHAP, fill=Group)) +
          geom_density(alpha=0.5) +
          custom_theme +
          ggtitle(paste("SHAP Density for", gene)))
  dev.off()

  # SHAP-表达量散点图
  shap_df <- data.frame(
    SHAP = shap_obj$S[, gene],
    Expr = train_set[, gene],
    Group = train_set$Group
  )
  pdf(file = file.path(result_dir, subdir_scatter, paste0(idx_str, "_", safe_gene, ".pdf")), width=pdf_width, height=pdf_height)
  print(ggplot(shap_df, aes(x=Expr, y=SHAP, color=Group)) +
          geom_point(alpha=0.7) +
          custom_theme +
          ggtitle(paste("SHAP vs Expression for", gene)))
  dev.off()
}

# ===================== 组合依赖图 - 每行4个 =====================
cat("组合所有依赖图（每行4个）...\n")
n_cols <- 4
n_rows <- ceiling(length(dependence_plots_list) / n_cols)
pdf(file=dependence_file, width=16, height=4*n_rows)
print(wrap_plots(dependence_plots_list, ncol = n_cols))
dev.off()
cat("依赖图已保存至：", dependence_file, "\n")

# ===================== Top2基因交互图已删除 =====================

# SHAP聚类热图
pdf(heatmap_file, width=pdf_width, height=pdf_height)
pheatmap(shap_obj$S, 
         cluster_rows=TRUE, cluster_cols=TRUE, 
         show_rownames=FALSE, show_colnames=TRUE, 
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100))
dev.off()

# ===================== 步骤8：输出基因 =====================
cat("步骤8/8: 输出基因...\n")
top20_genes <- head(shap_importance, 20)
write.csv(data.frame(Gene=names(top20_genes), Importance=as.numeric(top20_genes)), 
          top20_file, row.names=FALSE)

# SHAP特征累计贡献曲线
cum_contrib <- cumsum(abs(shap_importance)) / sum(abs(shap_importance))
pdf(file.path(result_dir, subdir_summary, "03_cumulative_contribution_curve.pdf"), width = pdf_width, height = pdf_height)
plot(seq_along(cum_contrib), cum_contrib, type = "s", lwd = 2,
     xlab = "Number of top features", ylab = "Cumulative SHAP Contribution",
     main = "Cumulative SHAP Contribution Curve")
abline(h = 0.8, col = "red", lty = 2)  # 80%参考线
dev.off()

#分组SHAP重要性条形图
group_levels <- levels(train_set$Group)
shap_group <- sapply(group_levels, function(grp) 
    colMeans(abs(shap_obj$S[train_set$Group == grp, , drop = FALSE]))
)

#分组SHAP重要性条形图
shap_group <- as.data.frame(shap_group)
shap_group$Gene <- rownames(shap_group)
library(reshape2)
shap_group_long <- reshape2::melt(shap_group, id.vars = "Gene", variable.name = "Group", value.name = "SHAP_Importance")

pdf(file.path(result_dir, subdir_importance, "03_shap_group_importance_barplot.pdf"), width = max(8, pdf_width*1.5), height = pdf_height)
print(
  ggplot(shap_group_long, aes(x = reorder(Gene, SHAP_Importance), y = SHAP_Importance, fill = Group)) +
    geom_bar(stat='identity', position = "dodge") +
    coord_flip() +
    labs(title="Group-wise SHAP Feature Importance", x="Gene", y="Mean(|SHAP|)") +
    custom_theme
)
dev.off()
feature_importance_df <- data.frame(
  Gene = names(shap_importance),
  MeanAbsSHAP = as.numeric(shap_importance)
)
write.csv(feature_importance_df, file.path(result_dir, subdir_summary, "feature_importance.csv"), row.names = FALSE)
group_levels <- levels(train_set$Group)
shap_group <- sapply(group_levels, function(grp)
  colMeans(abs(shap_obj$S[train_set$Group == grp, , drop = FALSE]))
)
shap_group_df <- data.frame(Gene = rownames(shap_group), shap_group)
write.csv(shap_group_df, file.path(result_dir, subdir_summary, "groupwise_importance.csv"), row.names = FALSE)

auc_table <- data.frame(
  Model = ml_methods$ModelName,
  AUC = sapply(ml_methods$MethodID, function(id) as.numeric(model_auc_map[[id]])),
  CI_lower = sapply(ml_methods$MethodID, function(id) roc_ci_list[[id]][1]),
  CI_upper = sapply(ml_methods$MethodID, function(id) roc_ci_list[[id]][3])
)
write.csv(auc_table, file.path(result_dir, subdir_summary, "model_auc_summary.csv"), row.names = FALSE)

sample_influence <- rowSums(abs(shap_obj$S))
influential_samples <- order(sample_influence, decreasing=TRUE)[1:10]
influence_df <- data.frame(
  Sample = rownames(train_set)[influential_samples],
  InfluenceScore = sample_influence[influential_samples]
)
write.csv(influence_df, file.path(result_dir, subdir_summary, "influential_samples.csv"), row.names=FALSE)

# 将y转为0/1
y_numeric <- ifelse(train_set$Group == "Treatment", 1, 0)

# DALEX explainer - 为caret模型提供自定义预测函数
explainer <- DALEX::explain(
  model = final_model,
  data = train_set[, -ncol(train_set)],
  y = y_numeric,
  predict_function = function(model, newdata) {
    pred <- predict(model, newdata, type="prob")
    pred[, "Treatment"]
  },
  label = "BestModel"
)
vi <- DALEX::model_parts(explainer, loss_function = DALEX::loss_one_minus_auc)
vi_df <- as.data.frame(vi)
pdf(file.path(result_dir, subdir_importance, "04_feature_permutation_importance.pdf"), width=pdf_width*1.2, height=pdf_height)
print(
  ggplot(vi_df[vi_df$variable != "_full_model_", ], aes(x=reorder(variable, -dropout_loss), y=dropout_loss)) +
    geom_bar(stat="identity", fill="tomato") +
    coord_flip() +
    labs(title="Permutation Feature Importance (AUC Drop)", x="Feature", y="AUC Drop") +
    custom_theme
)
dev.off()
write.csv(vi_df, file.path(result_dir, subdir_summary, "permutation_importance.csv"), row.names=FALSE)

# 瀑布图
pdf(file=waterfall_file, width=pdf_width, height=pdf_height)
print(sv_waterfall(shap_obj, row_id=30) + custom_theme)
dev.off()

# 力图
pdf(file=force_file, width=pdf_width, height=pdf_height)
print(sv_force(shap_obj, row_id=30) + custom_theme)
dev.off()

cat("全部流程已完成！\n")

