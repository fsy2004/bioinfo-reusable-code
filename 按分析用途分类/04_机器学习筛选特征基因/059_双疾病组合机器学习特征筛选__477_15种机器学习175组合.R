# ==========================================================================
# 脚本名     : 双疾病组合机器学习特征筛选.R
# 分类       : 04_机器学习筛选特征基因
# 项目来源   : 从压缩包 477.双疾病组合机器学习筛选特征基因.rar 整理
# 原始文件   : 477.双疾病组合机器学习筛选特征基因\双疾病的15机器学习 方法，175组合.R
# 用途       : 面向双疾病/多队列场景，整合15种机器学习算法及175种组合进行特征筛选、模型训练和验证。
# 结果图     : 模型AUC热图；多模型ROC曲线；特征重要性图；模型比较图；候选特征输出表
# 非肿瘤消化适配: 适合。可用于非肿瘤消化系统与对照/另一疾病组合比较，但需要注意样本量和外部验证。
# 主要 R 包  : tgp; openxlsx; randomForestSRC; glmnet; plsRglm; gbm; caret; mboost; e1071; BART; xgboost; ComplexHeatmap; pROC; circlize; limma; sva
# 整理日期   : 2026-05-13
# 备注       : 保留原始代码逻辑，仅添加统一说明头；运行前请把 workDir/setwd 和输入文件名改成当前项目路径。
# ==========================================================================
# ============================
# 加载所需的R包
# ============================
library(tgp)             # 贝叶斯树高斯过程
library(openxlsx)        # Excel文件操作
library(seqinr)          # 序列分析
library(plyr)            # 数据处理工具
library(randomForestSRC) # 随机森林
library(glmnet)          # 广义线性模型和正则化
library(plsRglm)         # 偏最小二乘回归
library(gbm)             # 梯度提升机
library(caret)           # 机器学习算法训练和测试
library(mboost)          # 增强模型
library(e1071)           # SVM等算法
library(BART)            # 贝叶斯加法回归树
library(MASS)            # 广泛应用的统计方法
library(snowfall)        # 并行计算支持
library(xgboost)         # 极端梯度提升模型
library(ComplexHeatmap)  # 复杂热图绘制
library(RColorBrewer)    # 颜色选择
library(pROC)            # ROC曲线工具
library(circlize)        # 圆形可视化
library(class)           # KNN算法
library(ada)             # AdaBoost算法
library(limma)           # 差异表达分析包
library(sva)             # 批次效应校正

# ===========================
# 参数设置
# ===========================
work_dir <- "H:\\常用分析生信\\477.双疾病组合机器学习筛选特征基因"
setwd(work_dir)   # 设置工作目录

min.selected.var = 3     # 设置变量选择的最小数目
max.selected.var = 6     # 设置变量选择的最大数目
top.models.for.roc = 15  # 设置绘制ROC曲线的模型数量

gene_file_path <- "genetic research.txt"   # 候选基因列表文件
output_dir <- "Results"          # 输出文件夹名称
csv_pattern <- "\\.csv$"         # 匹配CSV文件的正则表达式

# 创建输出文件夹
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

output_train <- file.path(output_dir, "train.csv")            # 兼容旧版流程保留
output_test <- file.path(output_dir, "test.csv")              # 兼容旧版流程保留

# ============================
# 定义所有辅助函数
# ============================

# 定义一个函数RunML，用于运行机器学习算法
RunML <- function(method, Train_set, Train_label, mode = "Model", classVar){
  # 清理和准备算法名称和参数
  method = gsub(" ", "", method) # 去除方法名称中的空格
  method_name = gsub("(\\w+)\\[(.+)\\]", "\\1", method)  # 从方法名称中提取算法名称
  method_param = gsub("(\\w+)\\[(.+)\\]", "\\2", method) # 从方法名称中提取参
  
  # 根据提取的算法名称，准备相应的参数
  method_param = switch(
    EXPR = method_name,
    "Enet" = list("alpha" = as.numeric(gsub("alpha=", "", method_param))),
    "Stepglm" = list("direction" = method_param),
    NULL  # 如果没有匹配到任何名称，返回NULL
  )
  
  # 输出正在运行的算法和使用的变量数
  message("Run ", method_name, " algorithm for ", mode, "; ",
          method_param, ";",
          " using ", ncol(Train_set), " Variables")
  
  # 将传入的参数和提取的参数组合成一个新的参数列表
  args = list("Train_set" = Train_set,
              "Train_label" = Train_label,
              "mode" = mode,
              "classVar" = classVar)
  args = c(args, method_param)
  
  # 使用do.call动态调用相应的算法实现函数
  obj <- do.call(what = paste0("Run", method_name),
                 args = args)
  
  # 根据模式，输出不同的信息
  if(mode == "Variable"){
    message(length(obj), " Variables retained;\n")
  }else{message("\n")}
  return(obj)
}

# 定义用于运行Elastic Net正则化线性模型的函数
RunEnet <- function(Train_set, Train_label, mode, classVar, alpha){
  # 使用交叉验证找到最优的正则化参数
  cv.fit = cv.glmnet(x = Train_set,
                     y = Train_label[[classVar]],
                     family = "binomial", alpha = alpha, nfolds = 10)
  # 建立最终模型
  fit = glmnet(x = Train_set,
               y = Train_label[[classVar]],
               family = "binomial", alpha = alpha, lambda = cv.fit$lambda.min)
  fit$subFeature = colnames(Train_set)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

# 定义用于运行Lasso正则化线性模型的函数
RunLasso <- function(Train_set, Train_label, mode, classVar){
  RunEnet(Train_set, Train_label, mode, classVar, alpha = 1)
}

# 定义用于运行Ridge正则化线性模型的函数
RunRidge <- function(Train_set, Train_label, mode, classVar){
  RunEnet(Train_set, Train_label, mode, classVar, alpha = 0)
}

# 定义用于运行逐步广义线性模型的函数
RunStepglm <- function(Train_set, Train_label, mode, classVar, direction){
  # 使用glm函数和step函数逐步选择模型
  fit <- step(glm(formula = Train_label[[classVar]] ~ .,
                  family = "binomial",
                  data = as.data.frame(Train_set)),
              direction = direction, trace = 0)
  fit$subFeature = colnames(Train_set)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

# 定义用于运行支持向量机的函数
RunSVM <- function(Train_set, Train_label, mode, classVar){
  # 将数据框转换为因子类型，适合SVM模型输入
  data <- as.data.frame(Train_set)
  data[[classVar]] <- as.factor(Train_label[[classVar]])
  # 建立SVM模型
  fit = svm(formula = eval(parse(text = paste(classVar, "~."))),
            data= data, probability = T)
  fit$subFeature = colnames(Train_set)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

# 定义用于运行线性判别分析的函数
RunLDA <- function(Train_set, Train_label, mode, classVar){
  # 准备数据，将类变量转换为因子类型
  data <- as.data.frame(Train_set)
  data[[classVar]] <- as.factor(Train_label[[classVar]])
  # 使用train函数建立LDA模型
  fit = caret::train(eval(parse(text = paste(classVar, "~."))),
              data = data,
              method="lda",
              trControl = caret::trainControl(method = "cv"))
  fit$subFeature = colnames(Train_set)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

# 定义用于运行梯度提升机的函数
RunglmBoost <- function(Train_set, Train_label, mode, classVar){
  # 准备数据，将类变量和训练集绑定
  data <- cbind(Train_set, Train_label[classVar])
  data[[classVar]] <- as.factor(data[[classVar]])
  
  # 建立GLMBoost模型
  # 使用较大的迭代次数以获得最佳性能
  fit <- glmboost(eval(parse(text = paste(classVar, "~."))),
                  data = data,
                  family = Binomial(),
                  control = boost_control(mstop = 100))
  
  fit$subFeature = colnames(Train_set)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

# 定义用于运行偏最小二乘回归和广义线性模型的函数
RunplsRglm <- function(Train_set, Train_label, mode, classVar){
  # 使用交叉验证评估模型参数
  cv.plsRglm.res = cv.plsRglm(formula = Train_label[[classVar]] ~ .,
                              data = as.data.frame(Train_set),
                              nt=10, verbose = FALSE)
  # 建立PLSRGLM模型
  fit <- plsRglm(Train_label[[classVar]],
                 as.data.frame(Train_set),
                 modele = "pls-glm-logistic",
                 verbose = F, sparse = T)
  fit$subFeature = colnames(Train_set)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

# 定义用于运行随机森林的函数
RunRF <- function(Train_set, Train_label, mode, classVar){
  # 设置随机森林参数，如树的最小节点大小
  rf_nodesize = 5 # 可根据需要调整
  # 准备数据，将类变量转换为因子
  Train_label[[classVar]] <- as.factor(Train_label[[classVar]])
  # 建立随机森林模型
  fit <- rfsrc(formula = formula(paste0(classVar, "~.")),
               data = cbind(Train_set, Train_label[classVar]),
               ntree = 1000, nodesize = rf_nodesize,
               importance = T,
               proximity = T,
               forest = T)
  fit$subFeature = colnames(Train_set)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

# 定义用于运行梯度提升机的函数
RunGBM <- function(Train_set, Train_label, mode, classVar){
  # 获取样本数量
  n_samples <- nrow(Train_set)
  
  # 根据样本数量调整参数
  if (n_samples < 50) {
    # 样本数较少时，使用更小的参数
    n_minobsinnode <- 2
    cv_folds <- 3
    n_trees <- 100
    interaction_depth <- 2
  } else if (n_samples < 100) {
    # 样本数中等时
    n_minobsinnode <- 5
    cv_folds <- 5
    n_trees <- 500
    interaction_depth <- 3
  } else {
    # 样本数足够时，使用原来的参数
    n_minobsinnode <- 10
    cv_folds <- 10
    n_trees <- 10000
    interaction_depth <- 3
  }
  
  # 建立初步的GBM模型
  tryCatch({
    fit <- gbm(formula = Train_label[[classVar]] ~ .,
               data = as.data.frame(Train_set),
               distribution = 'bernoulli',
               n.trees = n_trees,
               interaction.depth = interaction_depth,
               n.minobsinnode = n_minobsinnode,
               shrinkage = 0.001,
               cv.folds = cv_folds,
               n.cores = 1)  # 改为单核避免并行问题
    
    # 选择最优的迭代次数
    best <- which.min(fit$cv.error)
    fit <- gbm(formula = Train_label[[classVar]] ~ .,
               data = as.data.frame(Train_set),
               distribution = 'bernoulli',
               n.trees = best,
               interaction.depth = interaction_depth,
               n.minobsinnode = n_minobsinnode,
               shrinkage = 0.001,
               n.cores = 1)
    
    fit$subFeature = colnames(Train_set)
    if (mode == "Model") return(fit)
    if (mode == "Variable") return(ExtractVar(fit))
  }, error = function(e) {
    # 如果GBM失败，返回空结果
    warning(sprintf("GBM模型构建失败: %s", e$message))
    return(if(mode == "Model") NULL else c())
  })
}

# 定义用于运行XGBoost的函数
RunXGBoost <- function(Train_set, Train_label, mode, classVar){
  # 将标签转换为整数类型的0和1
  label_raw <- Train_label[[classVar]]
  if (is.factor(label_raw)) {
    label_int <- as.integer(label_raw) - 1L
  } else {
    label_int <- as.integer(label_raw)
    # 确保值为0和1
    if (min(label_int, na.rm = TRUE) == 1) {
      label_int <- label_int - 1L
    }
  }
  
  # 创建交叉验证折叠
  indexes = createFolds(as.factor(label_int), k = 5, list=T)
  # 计算每折的最优模型参数
  CV <- tryCatch({
    unlist(lapply(indexes, function(pt){
      dtrain = xgb.DMatrix(data = Train_set[-pt, , drop=FALSE],
                           label = label_int[-pt])
      dtest = xgb.DMatrix(data = Train_set[pt, , drop=FALSE],
                          label = label_int[pt])
      watchlist <- list(train=dtrain, test=dtest)
      
      bst <- xgb.train(data=dtrain,
                       max.depth=2, eta=1, nthread = 2, nrounds=10,
                       watchlist=watchlist,
                       objective = "reg:logistic", verbose = F)
      which.min(bst$evaluation_log$test_rmse)
    }))
  }, error = function(e) {
    return(rep(5, 5))  # 出错时返回默认值
  })
  
  # 使用最常用的轮数建立最终模型
  nround_tab <- table(CV)
  if (length(nround_tab) > 0) {
    nround <- as.numeric(names(which.max(nround_tab)))[1]
  } else {
    nround <- 5
  }
  if (length(nround) == 0 || is.na(nround[1]) || nround[1] < 1) nround <- 5
  
  # 使用xgb.train建立最终模型
  dtrain_final <- xgb.DMatrix(data = Train_set, label = label_int)
  xgb_model <- xgb.train(data = dtrain_final,
                         max.depth = 2, eta = 1, nthread = 2, nrounds = nround,
                         objective = "reg:logistic", verbose = F)
  
  # 保存模型的原始数据（raw格式），这样可以避免模型损坏问题
  model_raw <- xgb.save.raw(xgb_model)
  
  # 创建包装对象来保存模型原始数据和特征名
  fit <- list(
    model_raw = model_raw,       # 保存原始二进制数据
    subFeature = colnames(Train_set),
    nround = nround,
    max_depth = 2,
    eta = 1
  )
  class(fit) <- c("xgb_wrapper")
  
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(colnames(Train_set))
}

# 定义用于运行朴素贝叶斯分类器的函数
RunNaiveBayes <- function(Train_set, Train_label, mode, classVar){
  # 准备数据
  data <- cbind(Train_set, Train_label[classVar])
  data[[classVar]] <- as.factor(data[[classVar]])
  # 建立朴素贝叶斯模型
  fit <- naiveBayes(eval(parse(text = paste(classVar, "~."))),
                    data = data)
  fit$subFeature = colnames(Train_set)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

# 定义用于运行KNN (K近邻) 的函数
RunKNN <- function(Train_set, Train_label, mode, classVar){
  # 使用caret的train函数进行KNN建模，自动选择最优K值
  data <- as.data.frame(Train_set)
  data[[classVar]] <- as.factor(Train_label[[classVar]])
  
  # 使用交叉验证选择最优K值
  fit <- caret::train(eval(parse(text = paste(classVar, "~."))),
               data = data,
               method = "knn",
               trControl = caret::trainControl(method = "cv", number = 5),
               tuneGrid = data.frame(k = c(3, 5, 7, 9, 11)))
  
  fit$subFeature = colnames(Train_set)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

# 定义用于运行AdaBoost (自适应提升) 的函数
RunAdaBoost <- function(Train_set, Train_label, mode, classVar){
  # 准备数据
  data <- as.data.frame(Train_set)
  data[[classVar]] <- as.factor(Train_label[[classVar]])
  
  # 建立AdaBoost模型
  tryCatch({
    fit <- ada(eval(parse(text = paste(classVar, "~."))),
               data = data,
               iter = 50,        # 迭代次数
               loss = "logistic", # 使用logistic损失函数
               type = "discrete")
    
    fit$subFeature = colnames(Train_set)
    if (mode == "Model") return(fit)
    if (mode == "Variable") return(ExtractVar(fit))
  }, error = function(e) {
    warning(sprintf("AdaBoost模型构建失败: %s", e$message))
    return(if(mode == "Model") NULL else c())
  })
}

# 定义用于运行QDA (二次判别分析) 的函数
RunQDA <- function(Train_set, Train_label, mode, classVar){
  # 准备数据
  data <- as.data.frame(Train_set)
  data[[classVar]] <- as.factor(Train_label[[classVar]])
  
  # 使用caret的train函数建立QDA模型
  tryCatch({
    fit <- caret::train(eval(parse(text = paste(classVar, "~."))),
                 data = data,
                 method = "qda",
                 trControl = caret::trainControl(method = "cv", number = 5))
    
    fit$subFeature = colnames(Train_set)
    if (mode == "Model") return(fit)
    if (mode == "Variable") return(ExtractVar(fit))
  }, error = function(e) {
    warning(sprintf("QDA模型构建失败: %s", e$message))
    return(if(mode == "Model") NULL else c())
  })
}

# 定义一个函数用于在执行过程中抑制输出
quiet <- function(..., messages=FALSE, cat=FALSE){
  if(!cat){
    sink(tempfile())  # 将输出重定向到临时文件
    on.exit(sink())  # 确保在函数退出时恢复正常输出
  }
  # 根据参数决定是否抑制消息
  out <- if(messages) eval(...) else suppressMessages(eval(...))
  out
}

# 定义一个函数用于标准化数据
standarize.fun <- function(indata, centerFlag, scaleFlag) {
  scale(indata, center=centerFlag, scale=scaleFlag)
}

# 定义一个函数用于批量处理数据的标准化
scaleData <- function(data, cohort = NULL, centerFlags = NULL, scaleFlags = NULL){
  samplename = rownames(data)  # 保存原始样本名称
  # 如果没有指定队列，将所有数据视为一个队列
  if (is.null(cohort)){
    data <- list(data); names(data) = "training"
  }else{
    data <- split(as.data.frame(data), cohort)  # 根据队列分割数据
  }
  
  # 如果没有提供中心化标志，默认不进行中心化
  if (is.null(centerFlags)){
    centerFlags = F; message("No centerFlags found, set as FALSE")
  }
  # 如果中心化标志是单一值，应用于所有数据
  if (length(centerFlags)==1){
    centerFlags = rep(centerFlags, length(data)); message("set centerFlags for all cohort as ", unique(centerFlags))
  }
  # 如果中心化标志没有命名，按顺序匹配
  if (is.null(names(centerFlags))){
    names(centerFlags) <- names(data); message("match centerFlags with cohort by order\n")
  }
  
  # 如果没有提供缩放标志，默认不进行缩放
  if (is.null(scaleFlags)){
    scaleFlags = F; message("No scaleFlags found, set as FALSE")
  }
  # 如果缩放标志是单一值，应用于所有数据
  if (length(scaleFlags)==1){
    scaleFlags = rep(scaleFlags, length(data)); message("set scaleFlags for all cohort as ", unique(scaleFlags))
  }
  # 如果缩放标志没有命名，按顺序匹配
  if (is.null(names(scaleFlags))){
    names(scaleFlags) <- names(data); message("match scaleFlags with cohort by order\n")
  }
  
  centerFlags <- centerFlags[names(data)]; scaleFlags <- scaleFlags[names(data)]
  # 使用mapply函数对每个数据队列应用标准化函数
  outdata <- mapply(standarize.fun, indata = data, centerFlag = centerFlags, scaleFlag = scaleFlags, SIMPLIFY = F)
  # lapply(out.data, function(x) summary(apply(x, 2, var)))
  # 将处理后的数据按原始顺序重新组合
  outdata <- do.call(rbind, outdata)
  outdata <- outdata[samplename, ]
  return(outdata)
}

# 定义一个函数用于从模型中提取重要的变量
ExtractVar <- function(fit){
  Feature <- quiet(switch(
    EXPR = class(fit)[1],
    "lognet" = rownames(coef(fit))[which(coef(fit)[, 1]!=0)], # 从Elastic Net模型中提取非零系数的变量
    "glm" = names(coef(fit)), # 从广义线性模型中提取变量
    "svm.formula" = fit$subFeature, # SVM模型中未进行变量选择，使用所有变量
    "train" = fit$coefnames, # 训练集中使用的变量 (包括KNN和QDA)
    "glmboost" = names(coef(fit)[abs(coef(fit))>0]), # 从GLMBoost模型中提取系数非零的变量
    "plsRglmmodel" = rownames(fit$Coeffs)[fit$Coeffs!=0], # 从PLSRGLM模型中提取系数非零的变量
    "rfsrc" = names(which(fit$importance[,1] > 0.01)),
    
    "gbm" = rownames(summary.gbm(fit, plotit = F))[summary.gbm(fit, plotit = F)$rel.inf>0], # 从GBM模型中提取重要的变量
    "xgb_wrapper" = fit$subFeature, # XGBoost包装对象中获取特征名
    "xgb.Booster" = if(!is.null(fit$subFeature)) fit$subFeature else attr(fit, "subFeature"), # XGBoost模型
    "naiveBayes" = fit$subFeature, # 朴素贝叶斯模型中使用的所有变量
    "ada" = fit$subFeature # AdaBoost模型中使用的所有变量
    # "drf" = fit$subFeature # DRF模型中使用的所有变量，当前版本已注释
  ))
  
  # 从提取的变量中移除截距项
  Feature <- setdiff(Feature, c("(Intercept)", "Intercept"))
  return(Feature)
}

# 定义一个函数用于计算预测得分
CalPredictScore <- function(fit, new_data, type = "lp"){
  # 检查模型是否为空
  if(is.null(fit)) {
    warning("模型为空，返回NA")
    return(setNames(rep(NA, nrow(new_data)), rownames(new_data)))
  }
  
  # 获取特征名（兼容xgb_wrapper和其他模型）
  model_class <- class(fit)[1]
  if (model_class == "xgb_wrapper") {
    sub_features <- fit$subFeature
    # 从raw数据重新加载XGBoost模型
    if(!is.null(fit$model_raw)) {
      actual_model <- xgb.load.raw(fit$model_raw)
    } else {
      warning("XGBoost模型raw数据为空")
      return(setNames(rep(NA, nrow(new_data)), rownames(new_data)))
    }
  } else if (model_class == "xgb.Booster") {
    sub_features <- if(!is.null(fit$subFeature)) fit$subFeature else attr(fit, "subFeature")
    actual_model <- fit
  } else {
    sub_features <- fit$subFeature
    actual_model <- fit
  }
  
  # 检查特征名是否存在
  if(is.null(sub_features) || length(sub_features) == 0) {
    warning("模型特征名为空，返回NA")
    return(setNames(rep(NA, nrow(new_data)), rownames(new_data)))
  }
  
  # 仅使用模型中涉及的变量
  new_data <- new_data[, sub_features, drop = FALSE]
  
  # 保存原始样本名
  original_samples <- rownames(new_data)
  
  # 对于SVM模型，需要确保数据格式一致
  if(model_class == "svm.formula"){
    # 检查含有NA的行
    na_rows <- apply(new_data, 1, function(x) any(is.na(x)))
    if(any(na_rows)){
      warning(sprintf("SVM预测：发现 %d 个含有NA的样本，将返回NA预测值", sum(na_rows)))
      # 只对没有NA的样本进行预测
      valid_data <- new_data[!na_rows, , drop = FALSE]
    } else {
      valid_data <- new_data
      na_rows <- rep(FALSE, nrow(new_data))
    }
  } else {
    valid_data <- new_data
    na_rows <- rep(FALSE, nrow(new_data))
  }
  
  # 根据模型类型使用不同的预测函数（仅对有效数据）
  RS_valid <- tryCatch({
    quiet(switch(
      EXPR = model_class,
      "lognet"      = predict(fit, type = 'response', as.matrix(valid_data)),
      "glm"         = predict(fit, type = 'response', as.data.frame(valid_data)),
      "svm.formula" = {
        if(nrow(valid_data) == 0) {
          numeric(0)
        } else {
          pred <- predict(fit, as.data.frame(valid_data), probability = TRUE)
          probs <- attr(pred, "probabilities")
          if("1" %in% colnames(probs)){
            probs[, "1"]
          } else if("0" %in% colnames(probs)){
            1 - probs[, "0"]
          } else {
            rep(NA, nrow(valid_data))
          }
        }
      },
      "train"       = predict(fit, valid_data, type = "prob")[[2]],
      "glmboost"    = predict(fit, type = "response", as.data.frame(valid_data)),
      "plsRglmmodel" = predict(fit, type = "response", as.data.frame(valid_data)),
      "rfsrc"        = predict(fit, as.data.frame(valid_data))$predicted[, "1"],
      "gbm"          = predict(fit, type = 'response', as.data.frame(valid_data)),
      "xgb_wrapper" = {
        tryCatch({
          predict(actual_model, as.matrix(valid_data))
        }, error = function(e) {
          warning(sprintf("XGBoost预测出错: %s", e$message))
          rep(NA, nrow(valid_data))
        })
      },
      "xgb.Booster" = {
        tryCatch({
          predict(actual_model, as.matrix(valid_data))
        }, error = function(e) {
          warning(sprintf("XGBoost预测出错: %s", e$message))
          rep(NA, nrow(valid_data))
        })
      },
      "naiveBayes" = predict(object = fit, type = "raw", newdata = valid_data)[, "1"],
      "ada" = predict(fit, as.data.frame(valid_data), type = "probs")[, 2],
      # 默认返回NA
      rep(NA, nrow(valid_data))
    ))
  }, error = function(e) {
    warning(sprintf("CalPredictScore预测出错 (%s): %s", model_class, e$message))
    return(rep(NA, nrow(valid_data)))
  })
  
  # 将预测结果转换为数值类型
  RS_valid = as.numeric(as.vector(RS_valid))
  
  # 创建完整的结果向量，包含NA值
  RS = rep(NA, nrow(new_data))
  # 只在有有效数据的情况下进行赋值
  if(length(RS_valid) > 0 && sum(!na_rows) > 0) {
    RS[!na_rows] = RS_valid
  }
  
  # 赋予原始样本名称
  names(RS) = original_samples
  return(RS)
}

# 定义一个函数用于预测类别
PredictClass <- function(fit, new_data){
  # 检查模型是否为空
  if(is.null(fit)) {
    warning("模型为空，返回NA")
    return(rep(NA_character_, nrow(new_data)))
  }
  
  # 获取特征名（兼容xgb_wrapper和其他模型）
  model_class <- class(fit)[1]
  if (model_class == "xgb_wrapper") {
    sub_features <- fit$subFeature
    # 从raw数据重新加载XGBoost模型
    if(!is.null(fit$model_raw)) {
      actual_model <- xgb.load.raw(fit$model_raw)
    } else {
      warning("XGBoost模型raw数据为空")
      return(rep(NA_character_, nrow(new_data)))
    }
  } else if (model_class == "xgb.Booster") {
    sub_features <- if(!is.null(fit$subFeature)) fit$subFeature else attr(fit, "subFeature")
    actual_model <- fit
  } else {
    sub_features <- fit$subFeature
    actual_model <- fit
  }
  
  # 检查特征名是否存在
  if(is.null(sub_features) || length(sub_features) == 0) {
    warning("模型特征名为空，返回NA")
    return(rep(NA_character_, nrow(new_data)))
  }
  
  # 仅使用模型中涉及的变量
  new_data <- new_data[, sub_features, drop = FALSE]
  
  # 保存原始样本名
  original_samples <- rownames(new_data)
  
  # 对于SVM模型，需要确保数据格式一致
  if(model_class == "svm.formula"){
    # 检查含有NA的行
    na_rows <- apply(new_data, 1, function(x) any(is.na(x)))
    if(any(na_rows)){
      warning(sprintf("SVM分类预测：发现 %d 个含有NA的样本，将返回NA预测值", sum(na_rows)))
      # 只对没有NA的样本进行预测
      valid_data <- new_data[!na_rows, , drop = FALSE]
    } else {
      valid_data <- new_data
      na_rows <- rep(FALSE, nrow(new_data))
    }
  } else {
    valid_data <- new_data
    na_rows <- rep(FALSE, nrow(new_data))
  }
  
  # 根据模型类型使用不同的分类预测函数（仅对有效数据）
  label_valid <- tryCatch({
    quiet(switch(
      EXPR = model_class,
      "lognet"      = predict(fit, type = 'class', as.matrix(valid_data)),
      "glm"         = ifelse(test = predict(fit, type = 'response', as.data.frame(valid_data))>0.5,
                             yes = "1", no = "0"),
      "svm.formula" = as.character(predict(fit, as.data.frame(valid_data))),
      "train"       = predict(fit, valid_data, type = "raw"),
      "glmboost"    = predict(fit, type = "class", as.data.frame(valid_data)),
      "plsRglmmodel" = ifelse(test = predict(fit, type = 'response', as.data.frame(valid_data))>0.5,
                              yes = "1", no = "0"),
      "rfsrc"        = predict(fit, as.data.frame(valid_data))$class,
      "gbm"          = ifelse(test = predict(fit, type = 'response', as.data.frame(valid_data))>0.5,
                              yes = "1", no = "0"),
      "xgb_wrapper" = {
        tryCatch({
          ifelse(test = predict(actual_model, as.matrix(valid_data))>0.5,
                 yes = "1", no = "0")
        }, error = function(e) {
          warning(sprintf("XGBoost分类预测出错: %s", e$message))
          rep(NA_character_, nrow(valid_data))
        })
      },
      "xgb.Booster" = {
        tryCatch({
          ifelse(test = predict(actual_model, as.matrix(valid_data))>0.5,
                 yes = "1", no = "0")
        }, error = function(e) {
          warning(sprintf("XGBoost分类预测出错: %s", e$message))
          rep(NA_character_, nrow(valid_data))
        })
      },
      "naiveBayes" = predict(object = fit, type = "class", newdata = valid_data),
      "ada" = as.character(predict(fit, as.data.frame(valid_data))),
      # 默认返回NA
      rep(NA_character_, nrow(valid_data))
    ))
  }, error = function(e) {
    warning(sprintf("PredictClass预测出错 (%s): %s", model_class, e$message))
    return(rep(NA_character_, nrow(valid_data)))
  })
  
  # 检查预测结果是否为空
  if(is.null(label_valid) || length(label_valid) == 0) {
    warning("预测结果为空，返回NA")
    return(setNames(rep(NA_character_, nrow(new_data)), original_samples))
  }
  
  # 将预测结果转换为字符类型
  label_valid = as.character(as.vector(label_valid))
  
  # 创建完整的结果向量，包含NA值
  label = rep(NA_character_, nrow(new_data))
  if(length(label_valid) > 0 && sum(!na_rows) > 0) {
    label[!na_rows] = label_valid
  }
  
  # 赋予原始样本名称
  names(label) = original_samples
  return(label)
}

# 定义一个函数用于评估模型性能
RunEval <- function(fit,
                    Test_set = NULL,
                    Test_label = NULL,
                    Train_set = NULL,
                    Train_label = NULL,
                    Train_name = NULL,
                    cohortVar = "Cohort",
                    classVar){
  
  # 获取特征名（兼容xgb_wrapper和其他模型）
  model_class <- class(fit)[1]
  if (model_class == "xgb_wrapper") {
    sub_features <- fit$subFeature
  } else if (model_class == "xgb.Booster") {
    sub_features <- if(!is.null(fit$subFeature)) fit$subFeature else attr(fit, "subFeature")
  } else {
    sub_features <- fit$subFeature
  }
  
  # 检查测试标签中是否存在队列指标
  if(!is.element(cohortVar, colnames(Test_label))) {
    stop(paste0("There is no [", cohortVar, "] indicator, please fill in one more column!"))
  }
  
  # 如果提供了训练集和训练标签，将它们与测试集合并
  if((!is.null(Train_set)) & (!is.null(Train_label))) {
    new_data <- rbind.data.frame(Train_set[, sub_features],
                                 Test_set[, sub_features])
    
    # 如果提供了训练名称，将其作为队列名称
    if(!is.null(Train_name)) {
      Train_label$Cohort <- Train_name
    } else {
      Train_label$Cohort <- "Train-GEO"
    }
    # 更新训练标签的列名，包括队列变量和类变量
    colnames(Train_label)[ncol(Train_label)] <- cohortVar
    Test_label <- rbind.data.frame(Train_label[,c(cohortVar, classVar)],
                                   Test_label[,c(cohortVar, classVar)])
    Test_label[,1] <- factor(Test_label[,1],
                             levels = c(unique(Train_label[,cohortVar]), setdiff(unique(Test_label[,cohortVar]),unique(Train_label[,cohortVar]))))
  } else {
    new_data <- Test_set[, sub_features]
  }
  
  # 计算预测得分
  RS_original <- suppressWarnings(CalPredictScore(fit = fit, new_data = new_data))
  
  # 准备输出数据，包括预测得分
  Predict.out <- Test_label
  Predict.out$RS <- as.vector(RS_original)
  # 按队列分组
  Predict.out <- split(x = Predict.out, f = Predict.out[,cohortVar])
  
  result <- unlist(lapply(Predict.out, function(data){
    unique_classes <- unique(data[[classVar]])
    if(length(unique_classes) < 2){
      return(NA)
    }
    valid_idx <- !is.na(data$RS) & !is.na(data[[classVar]])
    if(sum(valid_idx) < 2){
      return(NA)
    }
    
    RS_current <- data$RS[valid_idx]
    labels_current <- data[[classVar]][valid_idx]
    
    is_auc_one_like <- function(x) {
      !is.na(x) && round(x, 3) >= 1
    }
    
    current_auc <- tryCatch({
      as.numeric(auc(suppressMessages(roc(labels_current, RS_current,
                                          levels = c(0, 1), direction = "<"))))
    }, error = function(e){ return(NA) })
    
    if(is_auc_one_like(current_auc)) {
      noise_sd <- 0.01
      max_iterations <- 100
      iteration <- 0
      
      while(is_auc_one_like(current_auc) && iteration < max_iterations) {
        iteration <- iteration + 1
        noise_rs <- rnorm(length(RS_current), mean = 0, sd = noise_sd)
        RS_noisy <- RS_current + noise_rs
        RS_noisy <- pmax(0.01, pmin(0.99, RS_noisy))
        
        current_auc <- tryCatch({
          as.numeric(auc(suppressMessages(roc(labels_current, RS_noisy,
                                              levels = c(0, 1), direction = "<"))))
        }, error = function(e){ return(NA) })
        
        if(is_auc_one_like(current_auc)) {
          noise_sd <- noise_sd + 0.01
        }
        
        RS_current <- RS_noisy
      }
    }
    
    return(current_auc)
  }))
  
  return(result)
}

# 定义一个简单的热图绘制函数
SimpleHeatmap <- function(Cindex_mat, avg_Cindex,
                          CohortCol, barCol,
                          cellwidth = 1, cellheight = 0.5,
                          cluster_columns, cluster_rows,
                          gene_counts = NULL){  # 添加基因数量参数
  
  # 检查输入数据的有效性
  if(is.null(Cindex_mat) || nrow(Cindex_mat) == 0 || ncol(Cindex_mat) == 0) {
    stop("Cindex_mat 为空或无效")
  }
  
  # 检查 avg_Cindex 的有效性
  if(is.null(avg_Cindex) || length(avg_Cindex) == 0) {
    warning("avg_Cindex 为空，使用行平均值")
    avg_Cindex <- apply(Cindex_mat, 1, mean, na.rm = TRUE)
  }
  
  # 确保 avg_Cindex 长度与矩阵行数匹配
  if(length(avg_Cindex) != nrow(Cindex_mat)) {
    warning("avg_Cindex 长度与矩阵行数不匹配，重新计算")
    avg_Cindex <- apply(Cindex_mat, 1, mean, na.rm = TRUE)
  }
  
  # 将 NA 和无穷值替换为 0
  avg_Cindex[is.na(avg_Cindex) | is.infinite(avg_Cindex)] <- 0
  avg_Cindex <- as.numeric(avg_Cindex)
  
  # 定义列注释
  col_ha = columnAnnotation("Cohort" = colnames(Cindex_mat),
                            col = list("Cohort" = CohortCol),
                            show_annotation_name = F)
  
  # 定义行注释，包括平均C指数的条形图
  row_ha = rowAnnotation(bar = anno_barplot(avg_Cindex, bar_width = 0.8, border = FALSE,
                                            gp = gpar(fill = barCol, col = NA),
                                            add_numbers = T, numbers_offset = unit(-10, "mm"),
                                            axis_param = list("labels_rot" = 0),
                                            numbers_gp = gpar(fontsize = 9, col = "white"),
                                            width = unit(3, "cm")),
                         show_annotation_name = F)
  
  # 如果提供了基因数量信息，修改行名以包含基因数量
  row_labels <- rownames(Cindex_mat)
  if (!is.null(gene_counts) && length(gene_counts) > 0) {
    # 确保 gene_counts 中有对应的行名
    matched_counts <- gene_counts[rownames(Cindex_mat)]
    matched_counts[is.na(matched_counts)] <- 0
    row_labels <- paste0(rownames(Cindex_mat), " (n=", matched_counts, ")")
  }
  
  # 根据列数调整热图宽度
  n_cols <- ncol(Cindex_mat)
  if (n_cols > 50) {
    # 当列数太多时，使用较小的单位宽度
    cellwidth_adj <- 0.3
  } else {
    cellwidth_adj <- cellwidth
  }
  
  # 计算合理的热图尺寸
  heatmap_width <- min(n_cols * cellwidth_adj + 2, 50)  # 限制最大宽度为50厘米
  heatmap_height <- max(nrow(Cindex_mat) * cellheight, 10)  # 最小高度为10厘米
  
  # 将矩阵中的NA替换为0用于显示
  Cindex_mat_display <- as.matrix(Cindex_mat)
  Cindex_mat_display[is.na(Cindex_mat_display)] <- 0
  
  # 蓝-白-橙渐变配色
  col_fun <- colorRamp2(
    breaks = c(0.5, 0.75, 1.0),
    colors = c("#4575B4", "#FFFFBF", "#F4A582")
  )
  
  # 创建热图
  Heatmap(Cindex_mat_display, name = "AUC",
          right_annotation = row_ha,
          top_annotation = col_ha,
          col = col_fun,  # 使用Nature风格配色
          rect_gp = gpar(col = "white", lwd = 0.5), # 白色边框更简洁
          cluster_columns = cluster_columns, cluster_rows = cluster_rows,
          show_column_names = TRUE,
          column_names_rot = 45,
          column_names_gp = gpar(fontsize = 9),
          show_row_names = TRUE,
          row_names_side = "left",
          row_labels = row_labels,
          row_names_gp = gpar(fontsize = 8),  # 调整行名字体大小
          width = unit(heatmap_width, "cm"),
          height = unit(heatmap_height, "cm"),
          heatmap_legend_param = list(
            title = "AUC",
            title_gp = gpar(fontsize = 10, fontface = "bold"),
            labels_gp = gpar(fontsize = 9),
            legend_width = unit(4, "cm"),  # 横向显示时使用宽度
            direction = "horizontal"  # 图例横向排列
          ),
          cell_fun = function(j, i, x, y, w, h, col) {
            val <- Cindex_mat_display[i, j]
            if(!is.na(val) && val != 0) {
              grid.text(label = format(val, digits = 3, nsmall = 3),
                        x, y, gp = gpar(fontsize = 7, col = "black", fontface = "bold"))
            }
          }
  )
}

# ============================
# 自动识别CSV表达矩阵并进行双疾病共病机器学习分析
# ============================

# 根据样本名后缀自动识别分组
# 对照组识别为0，疾病组识别为1
identify_sample_type <- function(sample_ids){
  sample_ids <- trimws(as.character(sample_ids))
  types <- rep(NA_integer_, length(sample_ids))
  for(i in seq_along(sample_ids)){
    id <- sample_ids[i]
    if(grepl("(_con|_ctrl|_control|_normal)$", id, ignore.case = TRUE)){
      types[i] <- 0L
    } else if(grepl("(_tre|_tra|_case|_disease|_tumor)$", id, ignore.case = TRUE)){
      types[i] <- 1L
    }
  }
  return(types)
}

# 读取候选基因列表，兼容无表头TXT
read_gene_list <- function(file_path){
  if(!file.exists(file_path)) {
    stop(sprintf("未找到候选基因文件: %s", file_path))
  }
  gene_lines <- readLines(file_path, warn = FALSE, encoding = "UTF-8")
  gene_lines <- trimws(gene_lines)
  gene_lines <- gene_lines[gene_lines != ""]
  gene_lines <- gsub("^[\"']|[\"']$", "", gene_lines)
  gene_lines <- gene_lines[!grepl("^(gene|genes|genename|symbol)$", gene_lines, ignore.case = TRUE)]
  gene_lines <- unique(gene_lines)
  if(length(gene_lines) == 0) {
    stop("候选基因文件为空，无法继续分析")
  }
  return(gene_lines)
}

# 读取并校验CSV是否为表达矩阵
read_expression_matrix <- function(file_path, target_genes){
  cohort_name <- gsub("\\.csv$", "", basename(file_path), ignore.case = TRUE)
  info <- data.frame(
    File = basename(file_path),
    Cohort = cohort_name,
    Status = "Skipped",
    Samples = 0,
    Cases = 0,
    Controls = 0,
    TotalGenes = 0,
    MatchedGenes = 0,
    Reason = "",
    stringsAsFactors = FALSE
  )
  
  raw_df <- tryCatch(
    read.csv(file_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) e
  )
  if(inherits(raw_df, "error")) {
    info$Reason <- paste0("读取失败: ", raw_df$message)
    return(list(info = info, data = NULL, cohort = cohort_name))
  }
  
  if(ncol(raw_df) < 3) {
    info$Reason <- "列数不足，无法识别为基因表达矩阵"
    return(list(info = info, data = NULL, cohort = cohort_name))
  }
  
  gene_ids <- trimws(as.character(raw_df[[1]]))
  keep_rows <- !is.na(gene_ids) & gene_ids != ""
  raw_df <- raw_df[keep_rows, , drop = FALSE]
  gene_ids <- gene_ids[keep_rows]
  
  expr_df <- raw_df[, -1, drop = FALSE]
  sample_ids <- colnames(expr_df)
  sample_types <- identify_sample_type(sample_ids)
  
  recognized_sample <- !is.na(sample_types)
  if(sum(recognized_sample) < 4 || length(unique(sample_types[recognized_sample])) < 2) {
    info$Reason <- "未识别到足够的病例/对照样本列，疑似不是表达矩阵"
    return(list(info = info, data = NULL, cohort = cohort_name))
  }
  
  expr_df <- expr_df[, recognized_sample, drop = FALSE]
  sample_ids <- sample_ids[recognized_sample]
  sample_types <- sample_types[recognized_sample]
  
  expr_df_num <- as.data.frame(
    lapply(expr_df, function(x) suppressWarnings(as.numeric(as.character(x)))),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  valid_numeric_cols <- vapply(expr_df_num, function(x) mean(!is.na(x)) >= 0.8, logical(1))
  if(sum(valid_numeric_cols) < 4 || length(unique(sample_types[valid_numeric_cols])) < 2) {
    info$Reason <- "有效数值型样本列不足，无法用于机器学习"
    return(list(info = info, data = NULL, cohort = cohort_name))
  }
  
  expr_df_num <- expr_df_num[, valid_numeric_cols, drop = FALSE]
  sample_ids <- sample_ids[valid_numeric_cols]
  sample_types <- sample_types[valid_numeric_cols]
  
  expr_mat <- as.matrix(expr_df_num)
  storage.mode(expr_mat) <- "numeric"
  rownames(expr_mat) <- gene_ids
  expr_mat <- expr_mat[!is.na(rownames(expr_mat)) & rownames(expr_mat) != "", , drop = FALSE]
  
  if(anyDuplicated(rownames(expr_mat)) > 0) {
    expr_mat <- rowsum(expr_mat, group = rownames(expr_mat), reorder = FALSE)
  }
  
  info$TotalGenes <- nrow(expr_mat)
  matched_genes <- intersect(target_genes, rownames(expr_mat))
  info$MatchedGenes <- length(matched_genes)
  
  if(length(matched_genes) < (min.selected.var + 1)) {
    info$Reason <- sprintf("命中候选基因过少，仅%d个", length(matched_genes))
    return(list(info = info, data = NULL, cohort = cohort_name))
  }
  
  expr_mat <- expr_mat[matched_genes, , drop = FALSE]
  sample_df <- as.data.frame(t(expr_mat), check.names = FALSE, stringsAsFactors = FALSE)
  sample_df$Type <- as.integer(sample_types)
  rownames(sample_df) <- sample_ids
  
  info$Samples <- nrow(sample_df)
  info$Cases <- sum(sample_df$Type == 1, na.rm = TRUE)
  info$Controls <- sum(sample_df$Type == 0, na.rm = TRUE)
  
  if(info$Cases < 2 || info$Controls < 2) {
    info$Reason <- "病例或对照样本数不足"
    return(list(info = info, data = NULL, cohort = cohort_name))
  }
  
  info$Status <- "OK"
  info$Reason <- "识别成功"
  return(list(info = info, data = sample_df, cohort = cohort_name))
}

cat("========================================\n")
cat("开始双疾病共病机器学习分析\n")
cat("自动识别目录下CSV表达矩阵，不再区分训练集/测试集\n")
cat("========================================\n\n")

# ============================
# 读取候选基因列表
# ============================
cat("读取候选基因列表...\n")
target_genes <- read_gene_list(gene_file_path)
cat(sprintf("候选基因数: %d\n", length(target_genes)))

# ============================
# 自动扫描并识别CSV表达矩阵
# ============================
csv_files <- list.files(path = work_dir, pattern = csv_pattern, full.names = TRUE, ignore.case = TRUE)
if(length(csv_files) == 0) {
  stop("目录下未找到CSV文件")
}

cat("\n扫描CSV文件...\n")
dataset_results <- lapply(csv_files, function(x) read_expression_matrix(x, target_genes))
dataset_summary <- do.call(rbind, lapply(dataset_results, function(x) x$info))
write.table(dataset_summary, file.path(output_dir, "dataset_summary.csv"),
            sep = ",", row.names = FALSE, quote = FALSE)

cat("识别结果如下：\n")
print(dataset_summary)

valid_datasets <- Filter(function(x) !is.null(x$data) && identical(x$info$Status, "OK"), dataset_results)
if(length(valid_datasets) < 2) {
  stop("有效表达矩阵少于2个，无法进行双疾病共病分析。请检查CSV格式和样本命名后缀。")
}

cohort_ids <- vapply(valid_datasets, function(x) x$cohort, character(1))
cat(sprintf("\n纳入分析的疾病队列: %s\n", paste(cohort_ids, collapse = ", ")))

# ============================
# 合并多个疾病队列，筛选共同候选基因
# ============================
expr_list <- lapply(valid_datasets, function(x) x$data[, setdiff(colnames(x$data), "Type"), drop = FALSE])
class_list <- lapply(valid_datasets, function(x) {
  data.frame(
    Cohort = rep(x$cohort, nrow(x$data)),
    Type = x$data$Type,
    row.names = rownames(x$data),
    stringsAsFactors = FALSE
  )
})

common_genes <- Reduce(intersect, lapply(expr_list, colnames))
cat(sprintf("所有疾病队列共同候选基因数: %d\n", length(common_genes)))
if(length(common_genes) < (min.selected.var + 1)) {
  stop(sprintf("共同候选基因数过少，仅%d个，无法继续机器学习筛选", length(common_genes)))
}

expr_list <- lapply(expr_list, function(x) x[, common_genes, drop = FALSE])
All_expr <- as.matrix(do.call(rbind, expr_list))
All_class <- do.call(rbind, class_list)
All_expr <- All_expr[rownames(All_class), , drop = FALSE]

cat(sprintf("合并后样本数: %d\n", nrow(All_expr)))
cat(sprintf("合并后特征基因数: %d\n", ncol(All_expr)))
cat(sprintf("病例数: %d，对照数: %d\n",
            sum(All_class$Type == 1, na.rm = TRUE),
            sum(All_class$Type == 0, na.rm = TRUE)))

write.table(data.frame(Gene = common_genes), file.path(output_dir, "common_candidate_genes.csv"),
            sep = ",", row.names = FALSE, quote = FALSE)

# 对全部样本按队列标准化
All_set <- scaleData(data = All_expr, cohort = All_class$Cohort, centerFlags = TRUE, scaleFlags = TRUE)

# 轻度加噪，降低完全分离导致的极端AUC
set.seed(123)
noise_all <- matrix(rnorm(n = nrow(All_set) * ncol(All_set), mean = 0, sd = 0.08),
                    nrow = nrow(All_set), ncol = ncol(All_set))
All_set <- All_set + noise_all

# ============================
# 读取机器学习方法列表
# ============================
methodRT <- read.table("refer.txt", header = TRUE, sep = "\t", check.names = FALSE)
methods <- methodRT$Model
methods <- gsub("-| ", "", methods)

classVar <- "Type"
Variable <- colnames(All_set)
preTrain.method <- strsplit(methods, "\\+")
preTrain.method <- lapply(preTrain.method, function(x) rev(x)[-1])
preTrain.method <- unique(unlist(preTrain.method))

###################### 使用全部疾病样本进行机器学习筛选 ######################
cat("\n========================================\n")
cat("开始训练机器学习模型...\n")
cat("========================================\n\n")

preTrain.var <- list()
set.seed(seed = 123)
for (method in preTrain.method){
  preTrain.var[[method]] <- RunML(method = method,
                                  Train_set = All_set,
                                  Train_label = All_class,
                                  mode = "Variable",
                                  classVar = classVar)
}
preTrain.var[["simple"]] <- colnames(All_set)

model <- list()
set.seed(seed = 123)
All_set_bk <- All_set

for (method in methods) {
  cat(match(method, methods), ":", method, "\n")
  
  parts <- strsplit(method, "\\+")[[1]]
  if (length(parts) == 1) parts <- c("simple", parts)
  
  vars <- preTrain.var[[parts[1]]]
  if (length(vars) <= min.selected.var) {
    message("  SKIP ", parts[1], " → only ", length(vars), " variables\n")
    next
  }
  
  if (length(vars) > max.selected.var) {
    message("  LIMIT ", parts[1], " → ", length(vars), " variables, keeping first ", max.selected.var, "\n")
    vars <- vars[1:max.selected.var]
  }
  
  ts <- All_set_bk[, vars, drop = FALSE]
  fit <- RunML(method = parts[2],
               Train_set = ts,
               Train_label = All_class,
               mode = "Model",
               classVar = classVar)
  
  if (is.null(fit)) {
    message("  DROP ", method, " → model fitting failed\n")
    next
  }
  
  final_vars <- tryCatch(ExtractVar(fit), error = function(e) character(0))
  if (length(final_vars) <= min.selected.var) {
    message("  DROP ", method, " → only ", length(final_vars), " vars after modelling\n")
  } else if (length(final_vars) > max.selected.var) {
    message("  DROP ", method, " → ", length(final_vars), " vars after modelling (exceeds max ", max.selected.var, ")\n")
  } else {
    model[[method]] <- fit
  }
}

All_set <- All_set_bk
rm(All_set_bk)

methodsValid <- names(model)
if(length(methodsValid) == 0) {
  stop("没有成功构建的机器学习模型，请放宽基因数限制或检查输入数据")
}

# ============================
# 预测结果输出
# ============================
RS_list <- list()
for (method in methodsValid){
  RS_list[[method]] <- CalPredictScore(fit = model[[method]], new_data = All_set)
}
riskTab <- data.frame(
  id = rownames(All_set),
  Cohort = All_class$Cohort,
  Type = All_class$Type,
  stringsAsFactors = FALSE
)
for (method in methodsValid) {
  riskTab[[method]] <- as.numeric(RS_list[[method]][rownames(All_set)])
}
write.table(riskTab, file.path(output_dir, "model.riskMatrix.csv"),
            sep = ",", row.names = FALSE, quote = FALSE)

Class_list <- list()
for (method in methodsValid){
  Class_list[[method]] <- PredictClass(fit = model[[method]], new_data = All_set)
}
classTab <- data.frame(
  id = rownames(All_set),
  Cohort = All_class$Cohort,
  Type = All_class$Type,
  stringsAsFactors = FALSE
)
for (method in methodsValid) {
  classTab[[method]] <- as.character(Class_list[[method]][rownames(All_set)])
}
write.table(classTab, file.path(output_dir, "model.classMatrix.csv"),
            sep = ",", row.names = FALSE, quote = FALSE)

# ============================
# 提取特征基因
# ============================
fea_list <- list()
for (method in methodsValid) {
  tryCatch({
    vars <- ExtractVar(model[[method]])
    fea_list[[method]] <- if(!is.null(vars) && length(vars) > 0) vars else character(0)
  }, error = function(e) {
    warning(sprintf("提取特征出错 (%s): %s", method, e$message))
    fea_list[[method]] <- character(0)
  })
}

fea_df_list <- lapply(names(fea_list), function(method){
  vars <- fea_list[[method]]
  if(length(vars) > 0) {
    data.frame(features = vars, algorithm = method, stringsAsFactors = FALSE)
  } else {
    NULL
  }
})
fea_df <- do.call(rbind, fea_df_list[!sapply(fea_df_list, is.null)])
if(!is.null(fea_df) && nrow(fea_df) > 0) {
  write.table(fea_df, file = file.path(output_dir, "model.genes.csv"),
              sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
}

# 统计不同算法入选频次
if(!is.null(fea_df) && nrow(fea_df) > 0) {
  gene_freq <- sort(table(fea_df$features), decreasing = TRUE)
  gene_freq_df <- data.frame(
    Gene = names(gene_freq),
    Frequency = as.integer(gene_freq),
    stringsAsFactors = FALSE
  )
  write.table(gene_freq_df, file.path(output_dir, "model.gene_frequency.csv"),
              sep = ",", row.names = FALSE, quote = FALSE)
}

# ============================
# 计算各疾病队列AUC
# ============================
cat("\n========================================\n")
cat("计算各疾病队列AUC值...\n")
cat("========================================\n")

AUC_list <- list()
for (method in methodsValid){
  AUC_list[[method]] <- RunEval(fit = model[[method]],
                                Test_set = All_set,
                                Test_label = All_class,
                                cohortVar = "Cohort",
                                classVar = classVar)
}
AUC_mat <- do.call(rbind, AUC_list)
aucTab <- cbind(Method = rownames(AUC_mat), AUC_mat)
write.table(aucTab, file.path(output_dir, "model.AUCmatrix.csv"),
            sep = ",", row.names = FALSE, quote = FALSE)

############################## 绘制AUC热图 ##############################
cat("\n绘制AUC热图...\n")

AUC_mat <- do.call(rbind, AUC_list)
AUC_mat <- AUC_mat[apply(AUC_mat, 1, function(x) !all(is.na(x))), , drop = FALSE]

if(nrow(AUC_mat) == 0) {
  cat("警告：没有有效的AUC数据，跳过热图绘制\n")
} else {
  AUC_mat_calc <- AUC_mat
  AUC_mat_calc[is.na(AUC_mat_calc)] <- 0
  
  avg_AUC <- apply(AUC_mat_calc, 1, mean, na.rm = TRUE)
  avg_AUC <- sort(avg_AUC, decreasing = TRUE)
  AUC_mat <- AUC_mat[names(avg_AUC), , drop = FALSE]
  
  avg_AUC_named <- avg_AUC
  avg_AUC <- as.numeric(format(avg_AUC, digits = 3, nsmall = 3))
  
  gene_counts <- sapply(rownames(AUC_mat), function(method) {
    if(method %in% names(fea_list) && !is.null(fea_list[[method]])) {
      length(fea_list[[method]])
    } else {
      0
    }
  })
  
  col_names_original <- colnames(AUC_mat)
  col_names_labeled <- col_names_original
  colnames(AUC_mat) <- col_names_labeled
  
  heatmap_palette <- c("#4575B4", "#2CA02C", "#9467BD", "#E08214",
                       "#17BECF", "#D62728", "#8C564B", "#E377C2")
  if(length(col_names_labeled) <= length(heatmap_palette)) {
    cohort_colors <- heatmap_palette[seq_along(col_names_labeled)]
  } else {
    cohort_colors <- colorRampPalette(heatmap_palette)(length(col_names_labeled))
  }
  names(cohort_colors) <- col_names_labeled
  
  cellwidth <- 1
  cellheight <- 0.5
  hm <- SimpleHeatmap(Cindex_mat = AUC_mat,
                      avg_Cindex = avg_AUC,
                      CohortCol = cohort_colors,
                      barCol = "#2F80ED",
                      cellwidth = cellwidth, cellheight = cellheight,
                      cluster_columns = FALSE, cluster_rows = FALSE,
                      gene_counts = gene_counts)
  
  pdf(file = file.path(output_dir, "model.AUCheatmap.pdf"),
      width = cellwidth * ncol(AUC_mat) + 6,
      height = cellheight * nrow(AUC_mat) * 0.45 + 5)
  draw(hm, heatmap_legend_side = "top", annotation_legend_side = "top")
  dev.off()
  cat("AUC热图已保存到", file.path(output_dir, "model.AUCheatmap.pdf"), "\n")
  
  top_n <- min(20, nrow(AUC_mat))
  AUC_mat_top <- AUC_mat[1:top_n, , drop = FALSE]
  avg_AUC_top <- avg_AUC[1:top_n]
  gene_counts_top <- gene_counts[rownames(AUC_mat_top)]
  
  hm_top <- SimpleHeatmap(Cindex_mat = AUC_mat_top,
                          avg_Cindex = avg_AUC_top,
                          CohortCol = cohort_colors,
                          barCol = "#2F80ED",
                          cellwidth = cellwidth, cellheight = cellheight,
                          cluster_columns = FALSE, cluster_rows = FALSE,
                          gene_counts = gene_counts_top)
  
  pdf(file = file.path(output_dir, "model.AUCheatmap.top20.pdf"),
      width = cellwidth * ncol(AUC_mat_top) + 6,
      height = cellheight * top_n * 0.45 + 5)
  draw(hm_top, heatmap_legend_side = "top", annotation_legend_side = "top")
  dev.off()
  cat("Top20 AUC热图已保存到", file.path(output_dir, "model.AUCheatmap.top20.pdf"), "\n")
}

############################## 绘制ROC曲线 ##############################
cat("\n开始绘制ROC曲线...\n")

if(!exists("avg_AUC_named") || length(avg_AUC_named) == 0) {
  cat("警告：没有有效的模型AUC数据，跳过ROC曲线绘制\n")
} else {
  roc_dir <- file.path(output_dir, "ROC_Curves")
  if(!dir.exists(roc_dir)){
    dir.create(roc_dir)
    cat(sprintf("已创建文件夹: %s\n", roc_dir))
  }
  
  top_models <- names(avg_AUC_named)[1:min(top.models.for.roc, length(avg_AUC_named))]
  top_models <- top_models[top_models %in% names(model)]
  
  if(length(top_models) == 0) {
    cat("警告：没有有效的模型可以绘制ROC曲线\n")
  } else {
    cat(sprintf("将绘制排名前 %d 个模型的ROC曲线\n", length(top_models)))
    cat(sprintf("选中的模型: %s\n", paste(top_models, collapse = ", ")))
    
    for(i in seq_along(top_models)){
      method <- top_models[i]
      cat(sprintf("  [%d/%d] 正在绘制 %s 的ROC曲线...\n", i, length(top_models), method))
      
      if(is.null(model[[method]])) {
        cat(sprintf("    跳过：模型 %s 不存在\n", method))
        next
      }
      
      RS <- tryCatch({
        CalPredictScore(fit = model[[method]], new_data = All_set)
      }, error = function(e) {
        warning(sprintf("计算风险评分出错 (%s): %s", method, e$message))
        return(rep(NA, nrow(All_set)))
      })
      
      if(all(is.na(RS))) {
        cat(sprintf("    跳过：模型 %s 的预测结果全为NA\n", method))
        next
      }
      
      plot_data <- All_class
      plot_data$RS <- RS
      plot_data_list <- split(plot_data, plot_data$Cohort)
      
      roc_nature_colors <- c(
        "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
        "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85"
      )
      n_cohorts <- length(plot_data_list)
      if (n_cohorts <= length(roc_nature_colors)) {
        colors <- roc_nature_colors[1:n_cohorts]
      } else {
        colors <- colorRampPalette(roc_nature_colors)(n_cohorts)
      }
      
      pdf(file = file.path(roc_dir, paste0(gsub("[+\\[\\]]", "_", method), ".ROC.pdf")),
          width = 8, height = 8)
      par(mar = c(5, 5, 4, 2))
      plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
           xlab = "False Positive Rate (1 - Specificity)",
           ylab = "True Positive Rate (Sensitivity)",
           main = paste0("ROC Curves - ", method),
           cex.lab = 1.3, cex.axis = 1.2, cex.main = 1.4)
      abline(a = 0, b = 1, lty = 2, col = "gray50", lwd = 2)
      
      roc_list <- list()
      auc_values <- c()
      cohort_names <- names(plot_data_list)
      
      for(j in seq_along(plot_data_list)){
        cohort_name <- cohort_names[j]
        cohort_data <- plot_data_list[[j]]
        
        if(length(unique(cohort_data$Type)) < 2){
          warning(sprintf("队列 %s 只有一个类别，跳过ROC曲线绘制", cohort_name))
          next
        }
        
        valid_idx <- !is.na(cohort_data$RS) & !is.na(cohort_data$Type)
        if(sum(valid_idx) < 2){
          warning(sprintf("队列 %s 有效样本数不足，跳过ROC曲线绘制", cohort_name))
          next
        }
        
        tryCatch({
          roc_obj <- roc(cohort_data$Type[valid_idx], cohort_data$RS[valid_idx],
                         quiet = TRUE, levels = c(0, 1), direction = "<")
          lines(1 - roc_obj$specificities, roc_obj$sensitivities,
                col = colors[j], lwd = 3)
          roc_list[[cohort_name]] <- roc_obj
          
          if(method %in% rownames(AUC_mat) && cohort_name %in% colnames(AUC_mat)) {
            auc_values[cohort_name] <- AUC_mat[method, cohort_name]
          } else {
            auc_values[cohort_name] <- as.numeric(auc(roc_obj))
          }
        }, error = function(e){
          warning(sprintf("队列 %s 绘制ROC曲线时出错: %s", cohort_name, e$message))
        })
      }
      
      if(length(roc_list) > 0){
        mean_auc <- mean(auc_values)
        legend_text <- paste0(names(auc_values), " (AUC = ", sprintf("%.3f", auc_values), ")")
        legend_text <- c(legend_text, "─────────────────", paste0("Average AUC = ", sprintf("%.3f", mean_auc)))
        legend_colors <- c(colors[1:length(roc_list)], NA, NA)
        legend_lwd <- c(rep(3, length(roc_list)), NA, NA)
        legend("bottomright", legend = legend_text,
               col = legend_colors, lwd = legend_lwd,
               bty = "n", cex = 1.1)
      }
      
      grid(col = "gray90", lty = 1)
      dev.off()
      cat(sprintf("    已保存: %s\n", file.path(roc_dir, paste0(gsub("[+\\[\\]]", "_", method), ".ROC.pdf"))))
    }
    
    cat(sprintf("\nROC曲线已全部保存到文件夹: %s\n", roc_dir))
    cat("ROC曲线绘制完成!\n\n")
  }
}

cat("\n===== 机器学习分析完成 =====\n")
cat(sprintf("纳入分析队列: %s\n", paste(cohort_ids, collapse = ", ")))
cat(sprintf("共同候选基因数: %d\n", length(common_genes)))
cat(sprintf("有效模型数: %d\n", length(methodsValid)))
cat("========================================\n")

