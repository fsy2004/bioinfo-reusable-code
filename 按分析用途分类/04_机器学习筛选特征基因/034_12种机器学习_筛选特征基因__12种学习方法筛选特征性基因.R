# =============================================================================
# 编号       : R034
# 脚本名     : 12种学习方法筛选特征性基因.R
# 分类       : 04_机器学习筛选特征基因
# 项目来源   : 12种机器学习_筛选特征基因
# 用途       : 整合 12 种机器学习方法筛选特征基因，输出各模型重要性、ROC 和交集图。
# 结果图     : ROC曲线；Venn图；UpSet图；条形图/柱状图；森林图；特征重要性图
# 主要 R 包  : C50; caret; DALEX; gbm; ggplot2; grid; gridExtra; kernelshap; kernlab; klaR; pdftools; png; pROC; randomForest; UpSetR; xgboost
# 整理时间   : 2026-05-10
# =============================================================================

library(UpSetR)


#导入包
library(caret)
library(DALEX)
library(ggplot2)
library(randomForest)
library(kernlab)
library(kernelshap)
library(pROC)
library(xgboost)
library(klaR)
library(gbm)
library(C50)

set.seed(123)      #设置随机数
expFile="Sample Type Matrix.csv"      #表达式数据文件
setwd("H:\\常用分析生信\\343.12种学习方法筛选特征性基因")      #设置工作目录

#读取表达式数据文件(CSV格式)
data=read.csv(expFile, header=T, check.names=F, row.names=1)
row.names(data)=gsub("-", "_", row.names(data))

#读取样品信息(对照组和实验组)
data=t(data)
group=gsub("(.*)_(con|tre)", "\\2", row.names(data))
group=gsub("con", "Control", group)
group=gsub("tre", "Treatment", group)
data=as.data.frame(data)


data$Type=as.factor(group)

# 检查数据维度
cat("Total samples:", nrow(data), "\n")
cat("Total features (genes + Type):", ncol(data), "\n")

#数据分割分析(训练集和测试集)
inTrain<-createDataPartition(y=data$Type, p=0.7, list=F)
train<-data[inTrain,]
test<-data[-inTrain,]

cat("Training samples:", nrow(train), "\n")
cat("Testing samples:", nrow(test), "\n")

#获取测试预测类别
yTestClass=test$Type
yTest=ifelse(yTestClass=="Control", 0, 1)

# 移除test的Type列用于预测特征提取
test_features<-test[,-ncol(test)]

#设置训练参数
control=trainControl(method="repeatedcv", number=5, savePredictions=TRUE)

#定义12种机器学习方法
methodRT <- data.frame(
  Name = c("Lasso", "RF", "SVM", "LDA", "GBM", "ElasticNet", "NeuralNet", "PLS", "AdaBoost", "Logistic", "NaiveBayes", "C5.0"),
  Method = c("glmnet", "rf", "svmRadial", "lda", "gbm", "glmnet", "nnet", "pls", "AdaBoost.M1", "LogitBoost", "nb", "C5.0"),
  stringsAsFactors = FALSE
)

#使用各种机器学习方法循环,建立机器学习模型
modelList=list()
AUC=c()
ROCcolor=rainbow(nrow(methodRT))
for(i in 1:nrow(methodRT)){
	name=methodRT[i,"Name"]
	method=methodRT[i,"Method"]

	#执行机器学习模型(训练集)
	tryCatch({
		cat("Processing:", name, "...\n")
		cat("  Train data shape:", nrow(train), "x", ncol(train), "\n")
		cat("  Train Type column class:", class(train$Type), "\n")
		if(name=="SVM"){
			model=train(Type ~ ., data = train, method=method, prob.model=TRUE, trControl = control, verbose=FALSE)
		} else if(name=="NeuralNet"){
			model=train(Type ~ ., data = train, method=method, trControl = control, verbose=FALSE, linout=FALSE, maxit=1000)
		} else if(name=="AdaBoost"){
			model=train(Type ~ ., data = train, method=method, trControl = control, verbose=FALSE)
		} else if(name=="C5.0"){
			model=train(Type ~ ., data = train, method=method, trControl = control, verbose=FALSE)
		} else if(name=="Lasso"){
			model=train(Type ~ ., data = train, method=method, trControl = control, verbose=FALSE, tuneGrid=expand.grid(alpha=1, lambda=seq(0.001, 0.1, length=10)))
		} else {
			model=train(Type ~ ., data = train, method=method, trControl = control, verbose=FALSE)
		}

		#得到每个样品预测的结果(概率值)
		pred=predict(model, newdata=test_features, type="prob")
		#绘制ROC曲线
		roc=roc(yTest, as.numeric(pred[,2]))
		AUC=c(AUC, paste0(name, ': ', sprintf("%.03f",roc$auc)))
		modelList[[method]]=as.numeric(roc$auc)
		if(i==1){
			pdf(file="ROC.pdf", width=5.5, height=5)
			plot(roc, print.auc=F, legacy.axes=T, main="", col=ROCcolor[i], lwd=3)
		}else{
			plot(roc, print.auc=F, legacy.axes=T, main="", col=ROCcolor[i], lwd=3, add=T)
		}

		#定义预测的函数
		p_fun=function(object, newdata){
			predict(object, newdata=newdata, type="prob")[,2]
		}

		#得到模型预测的结果
		explainer=explain(model, label = method,
							data = test, y = yTest,
							predict_function = p_fun,
							verbose = FALSE)

		#重要性分析,得到每个基因的重要性贡献度
		importance<-variable_importance(
		  explainer,
		  loss_function = loss_root_mean_square
		)

		#绘制重要性基因的图表(Nature风格)
		geneNum=10     #自定义展示基因数目
		pdf(file=paste0("importance.", name, ".pdf"), width=6, height=6)

		# 获取重要性数据
		# 过滤掉 baseline 和 Type
	imp_filtered <- importance[!grepl("baseline|full.model|^Type$", importance$variable), ]

	# 按基因聚合，计算平均dropout_loss
	imp_agg <- aggregate(dropout_loss ~ variable, data=imp_filtered, FUN=mean)
	imp_agg$label <- imp_filtered$label[1]

	# 获取top基因数据（按dropout_loss降序排序）
	if(nrow(imp_agg) > geneNum){
	  imp_sorted <- imp_agg[order(imp_agg$dropout_loss, decreasing=TRUE),]
	  imp_data <- imp_sorted[1:geneNum,]
	} else {
	  imp_data <- imp_agg
	}

		
		nature_colors <- c("#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
		                   "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#1a9850")
		p <- ggplot2::ggplot(imp_data, ggplot2::aes(x = reorder(variable, dropout_loss), y = dropout_loss)) +
		  ggplot2::geom_bar(stat = "identity", fill = "#17becf", color = "black", size = 0.5) +
		  ggplot2::coord_flip() +
		  ggplot2::theme_minimal() +
		  ggplot2::theme(
		    axis.title = ggplot2::element_text(size = 12, face = "bold", color = "black"),
		    axis.text = ggplot2::element_text(size = 10, color = "black"),
		    axis.line = ggplot2::element_line(color = "black", size = 0.8),
		    axis.ticks = ggplot2::element_line(color = "black", size = 0.8),
		    panel.border = ggplot2::element_rect(fill = NA, color = "black", size = 0.8),
		    panel.grid.major.x = ggplot2::element_line(color = "gray90"),
		    panel.grid.minor = ggplot2::element_blank(),
		    plot.title = ggplot2::element_text(size = 13, face = "bold", hjust = 0.5, color = "black")
		  ) +
		  ggplot2::labs(x = "Gene", y = "Feature Importance", title = paste0("Feature Importance - ", name))

		print(p)
		dev.off()

		#输出最重要的特定数量的基因
		geneNum=10
		importance$variable=gsub("_", "-", importance$variable)
		# 过滤掉 baseline 和 Type
	imp_filtered <- importance[!grepl("baseline|full.model|^Type$", importance$variable), ]

	# 按基因聚合，计算平均dropout_loss
	imp_agg <- aggregate(dropout_loss ~ variable, data=imp_filtered, FUN=mean)
	imp_agg$label <- imp_filtered$label[1]

	# 按dropout_loss降序排序，取最重要的前geneNum个基因
	if(nrow(imp_agg) > geneNum){
	  imp_sorted <- imp_agg[order(imp_agg$dropout_loss, decreasing=TRUE),]
	  imp_output <- imp_sorted[1:geneNum,]
	} else {
	  imp_output <- imp_agg
	}
	write.csv(imp_output, file=paste0("importanceGene.", name, ".csv"), row.names=F)
	}, error=function(e){
		cat("ERROR: ", name, " failed with message: ", e$message, "\n")
		cat("Trying alternative parameters for ", name, "...\n")
	})
}
legend('bottomright', AUC, col=ROCcolor, lwd=3, bty = 'n', cex=0.9)
dev.off()


#====================================================================
# 提取12种方法的关键基因，取交集，绘制Venn图
#====================================================================

# 定义12种方法的名称
methodNames <- c("Lasso", "RF", "SVM", "LDA", "GBM", "ElasticNet", "NeuralNet", "PLS", "AdaBoost", "Logistic", "NaiveBayes", "C5.0")

# 读取每种方法的基因列表
geneList <- list()
for(method in methodNames){
	fileName <- paste0("importanceGene.", method, ".csv")
	if(file.exists(fileName)){
		data <- read.csv(fileName, header=TRUE, stringsAsFactors=FALSE)
		geneList[[method]] <- data$variable
		cat("Loaded", method, ":", length(data$variable), "genes\n")
	} else {
		cat("Warning: File not found -", fileName, "\n")
	}
}

# 计算所有方法的交集基因
if(length(geneList) > 0){
	intersectGenes <- Reduce(intersect, geneList)
	cat("\n交集基因数量：", length(intersectGenes), "\n")
	cat("交集基因列表：\n")
	print(intersectGenes)

	# 保存交集基因到文件
	write.table(data.frame(Gene=intersectGenes), file="intersect_genes.txt", sep="\t", quote=FALSE, row.names=FALSE)
	cat("\n交集基因已保存到: intersect_genes.txt\n")
}

#====================================================================
# 绘制12种方法的UpSet图
#====================================================================
if(!require(UpSetR)){
	install.packages("UpSetR")
	library(UpSetR)
}

# 重新读取基因列表用于UpSet图
geneList_upset <- list()
for(method in methodNames){
	fileName <- paste0("importanceGene.", method, ".csv")
	if(file.exists(fileName)){
		data <- read.csv(fileName, header=TRUE, stringsAsFactors=FALSE)
		geneList_upset[[method]] <- data$variable
	}
}

if(length(geneList_upset) >= 2){
	# 创建所有基因的并集
	allGenes <- unique(unlist(geneList_upset))

	# 创建矩阵，每行是一个基因，每列是一个方法
	geneMatrix <- matrix(0, nrow=length(allGenes), ncol=length(geneList_upset))
	rownames(geneMatrix) <- allGenes
	colnames(geneMatrix) <- names(geneList_upset)

	for(method in names(geneList_upset)){
		geneMatrix[geneList_upset[[method]], method] <- 1
	}

	# 转换为UpSetR需要的格式
	upset_data <- as.data.frame(geneMatrix)

	# 绘制UpSet图
	pdf(file="UpSet_12methods.pdf", width=14, height=8)
	print(upset(upset_data,
		order.by = "freq",
		number.angles = 0,
		nsets = 12,
		nintersects = 40,
		sets = colnames(upset_data),
		mainbar.y.label = "Gene Intersection Size",
		sets.x.label = "Genes Per Method",
		point.size = 3,
		line.size = 1.2,
		mb.ratio = c(0.6, 0.4),
		text.scale = c(1.5, 1.3, 1.2, 1.2, 1.5, 1.2)
	))
	dev.off()

	# 提取第二页并覆盖原文件
	if(!require(pdftools)){
		install.packages("pdftools")
		library(pdftools)
	}
	if(!require(png)){
		install.packages("png")
		library(png)
	}

	upset_pdf <- "UpSet_12methods.pdf"
	pdf_info <- pdf_info(upset_pdf)
	if(pdf_info$pages >= 2){
		img <- pdf_render_page(upset_pdf, page = 2, dpi = 300)
		png_temp <- "temp_upset_page2.png"
		writePNG(img, png_temp)

		img_read <- readPNG(png_temp)
		pdf(file = upset_pdf, width = 14, height = 8)
		par(mar = c(0,0,0,0))
		plot(NA, xlim = c(0,1), ylim = c(0,1), xaxs = "i", yaxs = "i", axes = FALSE, xlab = "", ylab = "")
		rasterImage(img_read, 0, 0, 1, 1)
		dev.off()

		file.remove(png_temp)
		cat("UpSet图第二页已提取并覆盖原文件\n")
	}
	cat("\nUpSet图已保存: UpSet_12methods.pdf\n")
}

#====================================================================
# 汇总12种方法的重要性图为一个大图 (3行×4列)
#====================================================================
if(!require(pdftools)){
	install.packages("pdftools")
	library(pdftools)
}
if(!require(png)){
	install.packages("png")
	library(png)
}
if(!require(grid)){
	library(grid)
}
if(!require(gridExtra)){
	install.packages("gridExtra")
	library(gridExtra)
}

# 定义12种方法的名称（按顺序排列）
methodNames <- c("Lasso", "RF", "SVM", "LDA", "GBM", "ElasticNet",
                 "NeuralNet", "PLS", "AdaBoost", "Logistic", "NaiveBayes", "C5.0")

# 读取所有PDF并转换为图像
plot_list <- list()
for(i in 1:length(methodNames)){
	pdf_file <- paste0("importance.", methodNames[i], ".pdf")
	if(file.exists(pdf_file)){
		# 将PDF转换为PNG临时文件
		png_temp <- paste0("temp_", methodNames[i], ".png")
		img_bitmap <- pdf_render_page(pdf_file, page = 1, dpi = 150)
		writePNG(img_bitmap, png_temp)

		# 读取PNG并转换为raster对象
		img_png <- readPNG(png_temp)
		plot_list[[i]] <- rasterGrob(img_png, interpolate = TRUE)

		# 删除临时文件
		file.remove(png_temp)
		cat("已读取:", pdf_file, "\n")
	} else {
		# 如果文件不存在，创建空白占位符
		plot_list[[i]] <- textGrob(paste("Missing:", methodNames[i]))
		cat("文件不存在:", pdf_file, "\n")
	}
}

# 创建3×4的大图
if(length(plot_list) == 12){
	pdf(file = "importance_combined_3x4.pdf", width = 18, height = 24)
	grid.arrange(
		plot_list[[1]], plot_list[[2]], plot_list[[3]],
		plot_list[[4]], plot_list[[5]], plot_list[[6]],
		plot_list[[7]], plot_list[[8]], plot_list[[9]],
		plot_list[[10]], plot_list[[11]], plot_list[[12]],
		ncol = 3, nrow = 4
	)
	dev.off()
	cat("\n汇总图已保存: importance_combined_3x4.pdf\n")
}

#####
