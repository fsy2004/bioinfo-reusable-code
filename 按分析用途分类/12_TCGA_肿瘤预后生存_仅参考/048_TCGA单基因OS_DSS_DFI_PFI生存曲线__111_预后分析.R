# ==========================================================================
# 脚本名     : TCGA单基因OS_DSS_DFI_PFI生存曲线.R
# 分类       : 12_TCGA_肿瘤预后生存_仅参考
# 项目来源   : 从压缩包 111一个代码绘制单基因预后指标OS、DSS、DFI、PFI 生存曲线.rar 整理
# 原始文件   : 111一个代码绘制单基因预后指标OS、DSS、DFI、PFI 生存曲线\111预后分析.R
# 用途       : 对 TCGA 单基因表达进行 OS、DSS、DFI、PFI 四类生存终点的 Cox 分析和高低表达组 Kaplan-Meier 曲线绘制。
# 结果图     : OS生存曲线；DSS生存曲线；DFI生存曲线；PFI生存曲线；HR/95%CI/p值标注
# 非肿瘤消化适配: 肿瘤参考。非肿瘤消化系统通常没有 TCGA 生存终点，除非有随访/复发/事件时间数据，否则只借鉴图形和代码结构。
# 主要 R 包  : limma; survival; survminer; dplyr
# 整理日期   : 2026-05-13
# 备注       : 保留原始代码逻辑，仅添加统一说明头；运行前请把 workDir/setwd 和输入文件名改成当前项目路径。
# ==========================================================================
#################### 参数配置与包加载 ####################

# 参数设定区
gene_of_interest <- "JPH3"          # 指定探究基因名
control_sample_n <- 41              # 对照组样本数
tumor_sample_n <- 473               # 肿瘤组样本数
raw_exp_file <- "COAD.txt"          # 标准化表达数据文件名
survival_file <- "TCGASurvival.csv"     # 临床生存信息文件
data_id_tag <- "TCGA"               # 数据集标签(用于输出)
progress_total_steps <- 8           # 预定义进度条总数

# 设置工作目录（自行更改）
setwd("H:\\常用分析生信\\111一个代码绘制单基因预后指标OS、DSS、DFI、PFI 生存曲线")

# 包加载
library(limma)
library(survival)
library(survminer)
library(dplyr)       # 用于distinct去重
library(utils)       # 用于write.table等
library(grDevices)   # PDF输出

#################### 进度显示工具 ####################
show_progress <- function(step, total, msg){
  cat(sprintf("【%d/%d】%s\n", step, total, msg))
}

#################### 步骤一：表达矩阵预处理 ####################
show_progress(1, progress_total_steps, "正在读取并处理表达矩阵...")
# 读取基因表达原始数据
raw_matrix <- read.table(raw_exp_file, sep="\t", header=TRUE, check.names=FALSE)
# 转换为矩阵
mat_exprs <- as.matrix(raw_matrix)
# 行名设置为首列
rownames(mat_exprs) <- mat_exprs[, 1]
# 提取表达区，去掉第一列
mat_exprs <- mat_exprs[, -1]
# 重新赋行列名防止丢失
dimnames(mat_exprs) <- list(rownames(mat_exprs), colnames(mat_exprs))
# 全部转为数值
exp_num <- apply(mat_exprs, 2, as.numeric)
rownames(exp_num) <- rownames(mat_exprs)
colnames(exp_num) <- colnames(mat_exprs)

# 添加冗余判空校验
if(nrow(exp_num)==0 | ncol(exp_num)==0){ stop("表达矩阵有误！") }

# 样本分组信息
sample_categories <- c(rep("CTRL", control_sample_n), rep("TUMOR", tumor_sample_n))
# 选定的基因行合并上ID和分组信息, 防止ID偏移
gene_profile_data <- data.frame(SID=colnames(exp_num), GVAL=exp_num[gene_of_interest, ], GROUP=sample_categories)

# 列名重命名（更直观）
colnames(gene_profile_data) <- c("SampleID", gene_of_interest, "Group")
# 保存该基因表达
write.table(gene_profile_data, file="singleGeneProfile.txt", sep="\t", quote=FALSE, row.names=FALSE)

#################### 步骤二：样本名标准化处理 ####################
show_progress(2, progress_total_steps, "样本ID标准化处理...")
# 读取刚写好的文件
gene_profile_data <- read.table("singleGeneProfile.txt", sep="\t", header=TRUE, check.names=FALSE)
# 只取ID前15字符（保证样本一致性）
gene_profile_data$SampleID <- substr(gene_profile_data$SampleID, 1, 15)
# 替换写回
write.table(gene_profile_data, file="singleGeneProfile_trimmed.txt", sep="\t", quote=FALSE, row.names=FALSE)

#################### 步骤三：重复测序样本去重保存 ####################
show_progress(3, progress_total_steps, "剔除重复样本，仅保留独一...")
gene_no_dup <- read.table("singleGeneProfile_trimmed.txt", header=TRUE, sep="\t", check.names=FALSE)
gene_no_dup <- dplyr::distinct(gene_no_dup, SampleID, .keep_all=TRUE)
write.table(gene_no_dup, file="singleGeneProfile_unique.txt", sep="\t", quote=FALSE, row.names=FALSE)

#################### 步骤四：统一生存分析函数封装 ####################
show_progress(4, progress_total_steps, "准备生存分析函数定义...")

# ---- 主体封装 ----
do_single_survival <- function(time_name, status_name, y_label, pdffile, coxoutfile){
  # time_name/status_name: 生存天/生存状态列名字符串
  # y_label: PDF y轴说明
  # pdffile/coxoutfile: 输出PDF和cox表路径
  
  # 读取临床数据
  clin_df <- tryCatch({
    read.table(survival_file, header=TRUE, sep=",", check.names=FALSE, row.names=1)
  }, error=function(e){
    stop("无法读取生存临床数据文件！")
  })
  # 检查列名存在性
  if(!(time_name %in% colnames(clin_df)) || !(status_name %in% colnames(clin_df))){
    cat("临床文件缺少指定列！跳过...\n"); return(NULL)
  }
  clin_df_red <- clin_df[, c(time_name, status_name)]
  clin_df_red <- na.omit(clin_df_red)
  colnames(clin_df_red) <- c("surv_days", "surv_status")
  
  # 读取表达+分组文件
  gene_df <- read.table("singleGeneProfile_unique.txt", header=TRUE, sep="\t", check.names=FALSE, row.names=1)
  # 仅保留肿瘤
  if(!("Group" %in% colnames(gene_df))){ cat("表达文件缺Group列！"); return(NULL) }
  gene_df <- gene_df[gene_df$Group == "TUMOR", ]
  if(nrow(gene_df) == 0){ cat("肿瘤样本为零，跳过...\n"); return(NULL) }
  
  # 合并表达与临床公共样本
  comm_samples <- intersect(rownames(clin_df_red), rownames(gene_df))
  clin_df_red <- clin_df_red[comm_samples, , drop=FALSE]
  gene_df <- gene_df[comm_samples, , drop=FALSE]
  merge_df <- cbind(clin_df_red, gene_df)
  merge_df$years <- merge_df$surv_days/365
  gene_col <- gene_of_interest     # 这里gene_of_interest可以动态设定
  
  # 基因高低分组(可调换更复杂逻辑)
  merge_df$EXPGROUP <- ifelse(merge_df[, gene_col] > median(merge_df[, gene_col]), "High", "Low")
  
  # 冗余判断
  if(length(unique(merge_df$EXPGROUP))<2){ cat("分组缺失，跳过...\n"); return(NULL) }
  
  # 生存分析差异
  surv_diff_test <- survdiff(Surv(years, surv_status) ~ EXPGROUP, data=merge_df)
  pval_raw <- 1 - pchisq(surv_diff_test$chisq, df=1)
  # 格式化输出P值
  form_pval <- ifelse(pval_raw<0.001, "p < 0.001", paste0("p = ", sprintf("%.3f", pval_raw)))
  # COX回归
  cox_obj <- coxph(Surv(years, surv_status) ~ merge_df[, gene_col], data=merge_df)
  coxtable <- summary(cox_obj)
  HR_coef <- coxtable$conf.int[, "exp(coef)"]
  HR_LW <- coxtable$conf.int[, "lower .95"]
  HR_UP <- coxtable$conf.int[, "upper .95"]
  PCOX <- coxtable$coefficients[, "Pr(>|z|)"]
  
  # Cox结果数据表输出
  cox_out <- data.frame(code=data_id_tag, HR=HR_coef, HR95L=HR_LW, HR95H=HR_UP, p=PCOX)
  write.table(cox_out, file=coxoutfile, sep="\t", quote=FALSE, row.names=FALSE)
  
  # 生存曲线
  fitobj <- survfit(Surv(years, surv_status) ~ EXPGROUP, data=merge_df)
  HR_string <- sprintf("HR=%.2f (%.2f-%.2f)", HR_coef, HR_LW, HR_UP)
  
  # 进度显示
  cat(sprintf("当前处理: %s，样本数%d, 高组%d, 低组%d\n",
              y_label,
              nrow(merge_df),
              sum(merge_df$EXPGROUP=="High"),
              sum(merge_df$EXPGROUP=="Low")
  ))
  
  # 绘图
  surplt <- ggsurvplot(
    fitobj, 
    data=merge_df,
    conf.int=TRUE,
    pval=form_pval,
    pval.size=6,
    legend.labs=c("High", "Low"),
    legend.title=gene_col,
    xlab="Time (years)",
    ylab=y_label,
    break.time.by=1,
    risk.table.title="",
    palette=c("#F17C34", "#377EB8"),
    risk.table=TRUE,
    risk.table.height=0.28,
    ggtheme=theme_classic()
  )
  surplt$plot <- surplt$plot + annotate("text",
                                        x=max(merge_df$years)*0.88, y=max(fitobj$surv)*0.83,
                                        label=HR_string, size=4.8, color="black", hjust=1, vjust=1)
  pdf(file=pdffile, width=6.6, height=5.2)
  print(surplt)
  dev.off()
}

#################### 步骤五：OS分析 ####################
show_progress(5, progress_total_steps, "进行总生存(OS)分析...")
do_single_survival(time_name="OS.time", status_name="OS", 
                   y_label="Overall Survival", pdffile="OS_survplot.pdf", 
                   coxoutfile="TCGA_OS_cox.txt")

#################### 步骤六：DSS分析 ####################
show_progress(6, progress_total_steps, "进行DSS生存曲线分析...")
do_single_survival(time_name="DSS.time", status_name="DSS", 
                   y_label="Disease Specific Survival", pdffile="DSS_survplot.pdf", 
                   coxoutfile="TCGA_DSS_cox.txt")

#################### 步骤七：DFI分析 ####################
show_progress(7, progress_total_steps, "进行DFI生存曲线分析...")
do_single_survival(time_name="DFI.time", status_name="DFI", 
                   y_label="Disease Free Interval", pdffile="DFI_survplot.pdf", 
                   coxoutfile="TCGA_DFI_cox.txt")

#################### 步骤八：PFI生存试验 ####################
show_progress(8, progress_total_steps, "进行PFI生存曲线分析...")
do_single_survival(time_name="PFI.time", status_name="PFI", 
                   y_label="Progression Free Interval", pdffile="PFI_survplot.pdf", 
                   coxoutfile="TCGA_PFI_cox.txt")
# 函数：提取第二页并覆盖原PDF
extract_and_replace_second_page <- function(pdf_file) {
  # 检查文件是否存在
  if (!file.exists(pdf_file)) {
    cat("文件不存在：", pdf_file, "\n")
    return(FALSE)
  }
  # 临时输出文件
  temp_pdf <- tempfile(fileext = ".pdf")
  # 使用qpdf提取第2页的命令
  cmd <- sprintf('qpdf "%s" --pages . 2 -- "%s"', pdf_file, temp_pdf)
  # 执行系统命令
  result <- system(cmd)
  if (result != 0) {
    cat("qpdf命令运行失败：", pdf_file, "\n")
    return(FALSE)
  }
  # 覆盖原文件
  file.copy(temp_pdf, pdf_file, overwrite = TRUE)
  file.remove(temp_pdf)
  cat("[成功] 已提取", pdf_file, "的第2页\n")
  TRUE
}

# 步骤1：遍历当前目录下所有PDF文件（不区分名称，可定制pattern）
pdf_files <- list.files(pattern = "\\.pdf$", ignore.case = TRUE)

# 步骤2：循环批量处理
cat("发现", length(pdf_files), "个PDF文件。\n")
success_count <- 0
for (i in seq_along(pdf_files)) {
  cat(sprintf("【%d/%d】正在处理: %s\n", i, length(pdf_files), pdf_files[i]))
  res <- extract_and_replace_second_page(pdf_files[i])
  if (res) success_count <- success_count + 1
}
cat("完成。共提取并覆盖", success_count, "个PDF文件的第二页。\n")


#################### 结束提示 ####################
cat("全部生存分析已完成，结果已保存在当前目录下。\n")
