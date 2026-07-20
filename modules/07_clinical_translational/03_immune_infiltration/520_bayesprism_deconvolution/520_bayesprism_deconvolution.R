# =============================================================================
# 520 · BayesPrism bulk 反卷积 (Bayesian cell-type deconvolution, scRNA reference)
# -----------------------------------------------------------------------------
# 用 BayesPrism(真包,本机已装)以 scRNA 参考对 bulk RNA-seq 做贝叶斯反卷积,估计
# 各样本的细胞类型比例;合成 bulk 为已知 Dirichlet 混合 → 可【验证】估计 vs 真值。
# BayesPrism 优势:对参考-bulk 批次差异稳健、贝叶斯给后验,优于 NNLS/线性反卷积。
#
# ★诚实评估:反卷积务必与真值(或已知混合)对照,报 Pearson r + RMSE,不可只给比例图。
# Turnkey: Rscript 520_bayesprism_deconvolution.R   (合成参考+bulk→results/+assets/)
#          换数据: --reference ref_counts.csv --labels ref_labels.csv --bulk bulk_counts.csv
# 复用 _framework/theme_pub.R;无条形图(散点 + 热图)。需要 BayesPrism(见 README 服务器说明)。
# =============================================================================
suppressWarnings(suppressMessages({ library(ggplot2) }))

.this <- tryCatch({ a <- commandArgs(FALSE); m <- grep("^--file=", a)
  if(length(m)) dirname(normalizePath(sub("^--file=","",a[m[1]]))) else getwd() },
  error=function(e) getwd())
.fw <- NULL; .p <- .this
for(i in 1:6){ cand <- file.path(.p,"_framework","theme_pub.R"); if(file.exists(cand)){ .fw <- cand; break }; .p <- dirname(.p) }
if(!is.null(.fw)) source(.fw) else stop("需要 _framework/theme_pub.R")

if(!requireNamespace("BayesPrism", quietly=TRUE))
  stop("需要 BayesPrism:remotes::install_github('Danko-Lab/BayesPrism/BayesPrism')(见 README 服务器说明)")
suppressWarnings(suppressMessages(library(BayesPrism)))

set.seed(42)
DIR <- .this
DDAT <- file.path(DIR,"example_data"); DRES <- file.path(DIR,"results"); DAST <- file.path(DIR,"assets")
for(d in c(DDAT,DRES,DAST)) dir.create(d, showWarnings=FALSE, recursive=TRUE)

## ---- 1. 合成:scRNA 参考(4 类)+ bulk(已知 Dirichlet 混合)----------------
types <- c("Tcell","Bcell","Myeloid","Epithelial"); K<-length(types)
ng<-200; genes<-sprintf("G%03d",1:ng)
# 每类型一段 marker 程序的"真实表达谱"(基因均值)
prof <- matrix(rgamma(K*ng, shape=1.2, scale=1.0), K, ng, dimnames=list(types,genes))
for(i in 1:K) prof[i, ((i-1)*40+1):((i-1)*40+40)] <- prof[i, ((i-1)*40+1):((i-1)*40+40)] + rgamma(40,8,1)
# 参考:每类型 100 细胞,泊松抽样
ncell<-100; ref<-NULL; lab<-NULL
for(i in 1:K){ M<-matrix(rpois(ncell*ng, rep(prof[i,]*4, each=ncell)), ncell, ng, byrow=FALSE)
  ref<-rbind(ref,M); lab<-c(lab, rep(types[i], ncell)) }
rownames(ref)<-sprintf("cell%04d",1:nrow(ref)); colnames(ref)<-genes
# bulk:12 样本,每样本一个已知比例向量(Dirichlet),按比例混合 profile + 文库&噪声
nbulk<-12; frac_true<-matrix(rgamma(nbulk*K,2,1), nbulk, K); frac_true<-frac_true/rowSums(frac_true)
colnames(frac_true)<-types; rownames(frac_true)<-sprintf("bulk%02d",1:nbulk)
bulk<-matrix(0,nbulk,ng,dimnames=list(rownames(frac_true),genes))
for(s in 1:nbulk){ mu<-colSums(frac_true[s,]*prof)*800; bulk[s,]<-rpois(ng, mu) }
write.csv(ref, file.path(DDAT,"reference_counts.csv")); write.csv(data.frame(cell=rownames(ref),label=lab), file.path(DDAT,"reference_labels.csv"), row.names=FALSE)
write.csv(bulk, file.path(DDAT,"bulk_counts.csv")); write.csv(round(frac_true,3), file.path(DDAT,"true_fractions.csv"))
cat(sprintf("[gen] synthetic: ref %dx%d (4 types), bulk %dx%d with known fractions (demo only)\n", nrow(ref),ng,nbulk,ng))

## ---- 2. BayesPrism 反卷积 ---------------------------------------------------
myPrism <- new.prism(reference=ref, mixture=bulk, input.type="count.matrix",
                     cell.type.labels=lab, cell.state.labels=lab, key=NULL)
bp <- run.prism(prism=myPrism, n.cores=1)
theta <- get.fraction(bp=bp, which.theta="final", state.or.type="type")  # 样本 x 类型 估计比例
theta <- theta[rownames(frac_true), colnames(frac_true)]
write.csv(round(theta,4), file.path(DRES,"estimated_fractions.csv"))

## ---- 3. 验证:估计 vs 真值 -------------------------------------------------
val <- data.frame(sample=rep(rownames(frac_true),K),
                  celltype=factor(rep(colnames(frac_true), each=nbulk), levels=types),
                  true=as.vector(frac_true), est=as.vector(theta))
r  <- cor(val$true, val$est); rmse <- sqrt(mean((val$true-val$est)^2))
write.csv(data.frame(metric=c("Pearson_r","RMSE"), value=round(c(r,rmse),4)),
          file.path(DRES,"deconvolution_accuracy.csv"), row.names=FALSE)
cat(sprintf("[bayesprism] estimated vs true fractions: Pearson r=%.3f, RMSE=%.3f\n", r, rmse))

## ---- 4. 出图(散点验证 + 比例热图;无条形图)------------------------------
col <- pal_pub(K,"okabe_ito")
# Fig1: 估计 vs 真值 散点(y=x;按细胞类型着色)
p1 <- ggplot(val, aes(true, est, color=celltype)) +
  geom_abline(slope=1, intercept=0, linetype=2, color="grey55") +
  geom_point(size=2.4, alpha=0.85) +
  scale_color_manual(values=col, name="Cell type") +
  coord_equal(xlim=c(0,max(val$true,val$est)), ylim=c(0,max(val$true,val$est))) +
  labs(x="True fraction", y="BayesPrism estimate",
       title=sprintf("Deconvolution accuracy (r=%.2f, RMSE=%.3f)", r, rmse)) +
  theme_pub(base_size=11)
save_fig(p1, file.path(DAST,"accuracy_scatter"), width=5.2, height=4.0)

# Fig2: 估计比例热图(样本 x 类型)
hd <- data.frame(sample=rep(rownames(theta),K),
                 celltype=factor(rep(colnames(theta),each=nbulk),levels=types),
                 frac=as.vector(theta))
p2 <- ggplot(hd, aes(celltype, sample, fill=frac)) +
  geom_tile(color="white", linewidth=0.4) +
  scale_fill_viridis_c(option="D", name="Fraction") +
  labs(x=NULL, y=NULL, title="Estimated cell-type composition") +
  theme_pub(base_size=11) + theme(axis.text.x=element_text(angle=30, hjust=1))
save_fig(p2, file.path(DAST,"fraction_heatmap"), width=4.6, height=4.4)

cat("[fig] assets/: accuracy_scatter, fraction_heatmap (.pdf+.png)\n")
sink(file.path(DRES,"sessionInfo.txt")); print(sessionInfo()); sink()
