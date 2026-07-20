# =============================================================================
# 511 · 转录因子多证据收口 (regulon activity × JASPAR motif × DepMap essentiality)
# -----------------------------------------------------------------------------
# 把判断"某 TF 是否为真核心调控子"的三条正交证据收敛为一个可排序的收口分:
#   ① regulon 活性 = SCENIC/pySCENIC 的调控子 AUCell 平均活性(目标群);
#   ② motif 支持   = 该 TF 的 JASPAR 结合 motif 在其 regulon 靶基因启动子的命中比例;
#   ③ DepMap 必需性 = CRISPR gene-effect(越负越必需,转为正向"必需度")。
# 三者各自 rank 归一后取均值 = 收口分;三线齐高者=高可信核心 TF(抗单证据假阳性)。
#
# ★为 turnkey,三类证据为合成;接真数据时:regulon←pySCENIC(模块081),
#   motif←JASPAR2024/匹配,essentiality←DepMap CRISPR(Chronos gene effect)。见 README。
# Turnkey: Rscript 511_tf_convergence_depmap_jaspar.R   (合成→results/+assets/)
# 复用 _framework/theme_pub.R;无条形图(scatter+heatmap+lollipop)。
# =============================================================================
suppressWarnings(suppressMessages({ library(ggplot2) }))

.this <- tryCatch({ a <- commandArgs(FALSE); m <- grep("^--file=", a)
  if(length(m)) dirname(normalizePath(sub("^--file=","",a[m[1]]))) else getwd() },
  error=function(e) getwd())
.fw <- NULL; .p <- .this
for(i in 1:6){ cand <- file.path(.p,"_framework","theme_pub.R"); if(file.exists(cand)){ .fw <- cand; break }; .p <- dirname(.p) }
if(!is.null(.fw)) source(.fw) else stop("需要 _framework/theme_pub.R")

set.seed(42)
DIR <- .this
DDAT <- file.path(DIR,"example_data"); DRES <- file.path(DIR,"results"); DAST <- file.path(DIR,"assets")
for(d in c(DDAT,DRES,DAST)) dir.create(d, showWarnings=FALSE, recursive=TRUE)

rank01 <- function(x) (rank(x)-1)/(length(x)-1)   # rank 归一到 0-1

## ---- 1. 合成/读取三类证据 -------------------------------------------------
fT <- file.path(DDAT,"tf_evidence.csv")
if(!file.exists(fT)){
  tfs <- sprintf("TF%02d",1:12)
  # 设计:TF01-03 三线齐高(真核心);TF04-05 仅 1-2 线高(单证据陷阱);其余背景
  regulon  <- runif(12,0.2,0.5); motif <- runif(12,0.1,0.4); essential <- runif(12,-0.2,0.1)
  regulon[1:3]  <- runif(3,0.75,0.95); motif[1:3] <- runif(3,0.65,0.9); essential[1:3] <- runif(3,-1.1,-0.7)
  regulon[4]    <- 0.9                      # 仅 regulon 高(motif/essential 平)
  motif[5]      <- 0.85                     # 仅 motif 高
  essential[6]  <- -1.0                     # 仅 essential 高(泛必需,非特异调控)
  d <- data.frame(TF=tfs, regulon_activity=round(regulon,3),
                  motif_support=round(motif,3), depmap_gene_effect=round(essential,3))
  write.csv(d, fT, row.names=FALSE)
  cat("[gen] synthetic TF evidence (12 TFs; TF01-03 = true convergent) (demo only)\n")
}
d <- read.csv(fT)

## ---- 2. 收口分(三证据 rank 归一取均值;essential 取负→必需度越大越好)-----
d$essentiality <- -d$depmap_gene_effect                      # 越大越必需
d$r_regulon <- rank01(d$regulon_activity)
d$r_motif   <- rank01(d$motif_support)
d$r_essen   <- rank01(d$essentiality)
d$convergence <- rowMeans(d[,c("r_regulon","r_motif","r_essen")])
d <- d[order(-d$convergence),]
d$converged <- d$r_regulon>0.6 & d$r_motif>0.6 & d$r_essen>0.6  # 三线齐高
write.csv(d, file.path(DRES,"tf_convergence.csv"), row.names=FALSE)
cat(sprintf("[converge] 三线齐高的核心 TF: %s\n", paste(d$TF[d$converged], collapse=", ")))

## ---- 3. 出图(scatter + heatmap + lollipop;无条形图)---------------------
col <- pal_pub(3,"npg")
# Fig1: 收敛散点(motif x regulon,颜色=必需度,三线齐高者描边 + 标注)
p1 <- ggplot(d, aes(motif_support, regulon_activity)) +
  geom_point(aes(color=essentiality, size=convergence)) +
  geom_point(data=subset(d, converged), shape=21, size=6, stroke=1.1, color="black", fill=NA) +
  ggrepel::geom_text_repel(data=subset(d, converged), aes(label=TF), size=3.4, fontface="bold") +
  scale_color_viridis_c(option="C", name="Essentiality\n(-DepMap)") +
  scale_size(range=c(2,7), name="Convergence") +
  labs(x="JASPAR motif support (frac targets)", y="Regulon activity (AUCell)",
       title="TF convergence: motif × regulon × essentiality") +
  theme_pub(base_size=11)
save_fig(p1, file.path(DAST,"tf_convergence_scatter"), width=5.8, height=4.2)

# Fig2: 证据热图(TF x 三证据,各 rank 归一,viridis)
hm <- d[,c("TF","r_regulon","r_motif","r_essen")]
hd <- reshape(hm, direction="long", varying=list(2:4), v.names="score",
              times=c("Regulon","Motif","Essentiality"), timevar="evidence", idvar="TF")
hd$TF <- factor(hd$TF, levels=rev(d$TF))
hd$evidence <- factor(hd$evidence, levels=c("Regulon","Motif","Essentiality"))
p2 <- ggplot(hd, aes(evidence, TF, fill=score)) +
  geom_tile(color="white", linewidth=0.6) +
  scale_fill_viridis_c(option="D", name="Rank\nscore") +
  labs(x=NULL, y=NULL, title="Convergent evidence per TF") +
  theme_pub(base_size=11)
save_fig(p2, file.path(DAST,"evidence_heatmap"), width=4.4, height=4.4)

# Fig3: 收口分 lollipop(核心 TF 高亮)
d$TF <- factor(d$TF, levels=rev(d$TF))
p3 <- ggplot(d, aes(convergence, TF)) +
  geom_segment(aes(x=0, xend=convergence, yend=TF, color=converged), linewidth=1.2, alpha=0.7) +
  geom_point(aes(color=converged), size=4) +
  scale_color_manual(values=c(`FALSE`="grey65",`TRUE`=col[1]),
                     labels=c("partial","convergent (3/3)"), name=NULL) +
  labs(x="Convergence score (mean of 3 rank-normalised evidences)", y=NULL,
       title="TF convergence ranking") + theme_pub(base_size=11) +
  theme(legend.position=c(0.75,0.2))
save_fig(p3, file.path(DAST,"convergence_lollipop"), width=5.6, height=4.0)

cat("[fig] assets/: tf_convergence_scatter, evidence_heatmap, convergence_lollipop (.pdf+.png)\n")
sink(file.path(DRES,"sessionInfo.txt")); print(sessionInfo()); sink()
