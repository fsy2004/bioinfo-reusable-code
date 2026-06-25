# =============================================================================
# 503 · Cross-cohort generalization & honest robustness
# -----------------------------------------------------------------------------
# 泛化诚实三件套(多数生信不做,可信度大加分;思路来源:NETs 多组学范文):
#   (1) REML 随机效应 meta(metafor):signature 效应跨队列合并 + 森林图 + I^2/tau^2
#   (2) LODO 留一队列外推:轮流留一队列做独立测试,LASSO 训练其余 -> 各留一 AUC(诚实泛化)
#   (3) meta-weighted score:用 meta 合并效应作权重的跨平台稳健评分 -> 各队列分组分离
# Turnkey: Rscript 503_generalization_robustness.R  (默认读 example_data/, 写 results/+assets/)
# 复用 _framework/theme_pub.R;代码注释中文,图中文字英文。
# =============================================================================
suppressWarnings(suppressMessages({
  library(metafor); library(glmnet); library(pROC); library(ggplot2)
}))

## ---- 定位脚本目录 + 载入框架 ----------------------------------------------
.this <- tryCatch({ a <- commandArgs(FALSE); m <- grep("^--file=", a)
  if(length(m)) dirname(normalizePath(sub("^--file=","",a[m[1]]))) else getwd() },
  error=function(e) getwd())
.fw <- NULL; .p <- .this
for(i in 1:6){ cf <- file.path(.p,"_framework","theme_pub.R"); if(file.exists(cf)){ .fw <- cf; break }; .p <- dirname(.p) }
if(!is.null(.fw)){ source(.fw) } else {
  theme_pub <- function(base_size=11, ...) theme_bw(base_size=base_size)
  pal_pub   <- function(n=NULL, name="npg") scales::hue_pal()(ifelse(is.null(n),6,n))
  save_fig  <- function(plot, file, width=7, height=6, dpi=300){
    ggsave(paste0(file,".pdf"), plot, width=width, height=height); ggsave(paste0(file,".png"), plot, width=width, height=height, dpi=dpi) }
}

set.seed(42)
DIR <- .this; DDAT <- file.path(DIR,"example_data"); DRES <- file.path(DIR,"results"); DAST <- file.path(DIR,"assets")
for(d in c(DDAT,DRES,DAST)) dir.create(d, showWarnings=FALSE, recursive=TRUE)

## ---- 1. 读 / 生成多队列示例(模拟跨平台异质) -------------------------------
# 3 个队列,各 con/tre 两组,共享 6 个信号基因(效应量队列间略异 + 队列基线 shift)。
fM <- file.path(DDAT,"cohorts.rds")
if(!file.exists(fM)){
  mk <- function(ns=60, shift=0, eff=1.2){
    ng <- 40; genes <- sprintf("GENE%02d",1:ng); samp <- sprintf("S%02d",1:ns)
    grp <- rep(c("con","tre"), each=ns/2)
    M <- matrix(rnorm(ng*ns, 8+shift, 1.3), ng, ns, dimnames=list(genes,samp))
    sig <- genes[1:6]; M[sig, grp=="tre"] <- M[sig, grp=="tre"] + eff   # 信号:队列间效应异质
    list(expr=M, group=grp)
  }
  cohorts <- list(C1=mk(60,0,1.4), C2=mk(56,0.6,1.0), C3=mk(64,-0.4,1.2))  # 不同样本量/基线/效应
  saveRDS(cohorts, fM); cat("[gen] 3 synthetic cohorts (shared 6-gene signature, heterogeneous effects)\n")
}
cohorts <- readRDS(fM)
SIG <- sprintf("GENE%02d",1:6)   # signature 基因(实际用时换成你的 signature)

## ---- 2. REML 随机效应 meta:signature 效应跨队列合并 -----------------------
# 每队列:signature score = 6 基因 z-score 均值;算 con vs tre 的标准化均差 SMD。
smd <- lapply(names(cohorts), function(cn){
  ch <- cohorts[[cn]]; sc <- as.numeric(rowMeans(scale(t(ch$expr[SIG,]))))  # 每样本 signature score
  g  <- factor(ch$group)
  es <- escalc(measure="SMD",
               m1i=mean(sc[g=="tre"]), m2i=mean(sc[g=="con"]),
               sd1i=sd(sc[g=="tre"]),  sd2i=sd(sc[g=="con"]),
               n1i=sum(g=="tre"),      n2i=sum(g=="con"))
  data.frame(cohort=cn, yi=as.numeric(es$yi), vi=as.numeric(es$vi))
})
smd <- do.call(rbind, smd)
fit <- tryCatch(rma(yi, vi, data=smd, method="REML"),
                error=function(e) rma(yi, vi, data=smd, method="DL"))
cat(sprintf("[meta] pooled SMD=%.2f (95%%CI %.2f-%.2f), I^2=%.0f%%, tau^2=%.3f, p=%.1e\n",
            fit$b, fit$ci.lb, fit$ci.ub, fit$I2, fit$tau2, fit$pval))
# 森林图(用框架配色,自绘 ggplot 版)
fp <- data.frame(cohort=c(smd$cohort,"Pooled (REML)"),
                 est=c(smd$yi, as.numeric(fit$b)),
                 lo =c(smd$yi-1.96*sqrt(smd$vi), fit$ci.lb),
                 hi =c(smd$yi+1.96*sqrt(smd$vi), fit$ci.ub),
                 pooled=c(rep(FALSE,nrow(smd)), TRUE))
fp$cohort <- factor(fp$cohort, levels=rev(fp$cohort))
pf <- ggplot(fp, aes(est, cohort, color=pooled)) +
  geom_vline(xintercept=0, linetype=2, color="grey60") +
  geom_errorbarh(aes(xmin=lo, xmax=hi), height=0.2) + geom_point(aes(size=pooled)) +
  scale_color_manual(values=c("FALSE"=pal_pub(2)[1],"TRUE"=pal_pub(2)[2]), guide="none") +
  scale_size_manual(values=c("FALSE"=2.5,"TRUE"=4), guide="none") +
  labs(x="Standardized mean difference (signature)", y=NULL,
       title=sprintf("REML meta-analysis (I^2=%.0f%%)", fit$I2)) + theme_pub(base_size=11)
save_fig(pf, file.path(DAST,"meta_forest"), width=6.5, height=3.6)

## ---- 3. LODO 留一队列外推(LASSO) ------------------------------------------
auc_lodo <- sapply(names(cohorts), function(test){
  tr <- setdiff(names(cohorts), test)
  Xtr <- do.call(cbind, lapply(tr, function(c) cohorts[[c]]$expr)); Xtr <- t(Xtr)
  ytr <- factor(unlist(lapply(tr, function(c) cohorts[[c]]$group)))
  Xte <- t(cohorts[[test]]$expr); yte <- factor(cohorts[[test]]$group)
  cv  <- cv.glmnet(Xtr, ytr, family="binomial", alpha=1)
  pr  <- as.numeric(predict(cv, Xte, s="lambda.min", type="response"))
  as.numeric(auc(roc(yte, pr, quiet=TRUE, levels=c("con","tre"), direction="<")))
})
lodo <- data.frame(held_out=c(names(auc_lodo),"mean"), AUC=c(auc_lodo, mean(auc_lodo)))
write.csv(lodo, file.path(DRES,"LODO_AUC.csv"), row.names=FALSE)
cat(sprintf("[LODO] leave-one-dataset-out AUC: %s | mean=%.3f\n",
            paste(sprintf("%s=%.3f",names(auc_lodo),auc_lodo),collapse=", "), mean(auc_lodo)))
pl <- ggplot(lodo, aes(reorder(held_out,AUC), AUC, fill=held_out=="mean")) + geom_col(width=0.65) +
  geom_text(aes(label=sprintf("%.3f",AUC)), hjust=-0.1, size=3.3) +
  geom_hline(yintercept=0.5, linetype=2, color="grey60") + coord_flip(ylim=c(0,1)) +
  scale_fill_manual(values=c("FALSE"=pal_pub(2)[1],"TRUE"=pal_pub(2)[2]), guide="none") +
  labs(x="Held-out cohort", y="Test AUC (LASSO trained on the rest)", title="LODO generalization") + theme_pub(base_size=11)
save_fig(pl, file.path(DAST,"LODO_auc"), width=6, height=3.4)

## ---- 4. meta-weighted score(跨平台稳健评分) -------------------------------
# 每基因跨队列 meta 合并效应作权重,对样本 z-score 加权求和。
wgt <- sapply(rownames(cohorts[[1]]$expr), function(g){
  d <- lapply(names(cohorts), function(cn){ ch<-cohorts[[cn]]; x<-ch$expr[g,]; gp<-factor(ch$group)
    es<-escalc(measure="SMD", m1i=mean(x[gp=="tre"]), m2i=mean(x[gp=="con"]), sd1i=sd(x), sd2i=sd(x),
               n1i=sum(gp=="tre"), n2i=sum(gp=="con")); c(es$yi, es$vi) })
  d <- do.call(rbind, d); as.numeric(rma(d[,1], d[,2], method="REML")$b)
})
ms <- do.call(rbind, lapply(names(cohorts), function(cn){
  ch<-cohorts[[cn]]; z<-t(scale(t(ch$expr))); score<-colSums(z*wgt[rownames(z)])
  data.frame(cohort=cn, group=ch$group, score=as.numeric(score)) }))
write.csv(data.frame(gene=names(wgt), meta_weight=as.numeric(wgt)), file.path(DRES,"meta_weights.csv"), row.names=FALSE)
pm <- ggplot(ms, aes(cohort, score, fill=group)) +
  geom_boxplot(outlier.size=0.6, width=0.6) +
  scale_fill_manual(values=pal_pub(2), name="group") +
  labs(x="Cohort", y="Meta-weighted score", title="Meta-weighted score separates groups across cohorts") + theme_pub(base_size=11)
save_fig(pm, file.path(DAST,"meta_weighted_score"), width=6, height=3.8)

cat("[fig] assets/meta_forest, LODO_auc, meta_weighted_score (.pdf/.png)\n")
