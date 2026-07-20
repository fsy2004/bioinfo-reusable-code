# =============================================================================
# 516 · 复合多 panel 图 (composite "Figure 1") — UMAP + 火山 + 热图 + 森林
# -----------------------------------------------------------------------------
# 演示如何把多种图型拼成一张带 A/B/C/D 角标的发表级复合图(顶刊 Figure 1 的做法):
#   A UMAP 聚类  B 火山图(DEG)  C 顶 DEG 热图(行 z-score)  D 关键基因效应森林。
# 每个 panel 用框架统一主题,compose_panels() 原生拼版(非贴位图,避免拉伸失真)。
# 既给"复合多panel"模板,也含"火山+热图组合"。
#
# Turnkey: Rscript 516_composite_multipanel.R   (全合成→assets/)
# 复用 _framework/theme_pub.R(theme_pub/compose_panels/save_fig);无条形图。
# =============================================================================
suppressWarnings(suppressMessages({ library(ggplot2); library(ggrepel) }))

.this <- tryCatch({ a <- commandArgs(FALSE); m <- grep("^--file=", a)
  if(length(m)) dirname(normalizePath(sub("^--file=","",a[m[1]]))) else getwd() },
  error=function(e) getwd())
.fw <- NULL; .p <- .this
for(i in 1:6){ cand <- file.path(.p,"_framework","theme_pub.R"); if(file.exists(cand)){ .fw <- cand; break }; .p <- dirname(.p) }
if(!is.null(.fw)) source(.fw) else stop("需要 _framework/theme_pub.R")

set.seed(42)
DIR <- .this
DRES <- file.path(DIR,"results"); DAST <- file.path(DIR,"assets")
for(d in c(DRES,DAST)) dir.create(d, showWarnings=FALSE, recursive=TRUE)

## ---- Panel A: UMAP 聚类散点 -----------------------------------------------
nc<-600; k<-4; ctr<-matrix(c(-4,-4, 4,-3, -3,4, 4,4), k, 2, byrow=TRUE)
cl <- sample(1:k, nc, replace=TRUE)
um <- data.frame(UMAP1=ctr[cl,1]+rnorm(nc,0,1.1), UMAP2=ctr[cl,2]+rnorm(nc,0,1.1),
                 cluster=factor(paste0("C",cl)))
pA <- ggplot(um, aes(UMAP1, UMAP2, color=cluster)) +
  geom_point(size=0.9, alpha=0.8) +
  scale_color_manual(values=pal_pub(k,"okabe_ito"), name="Cluster") +
  labs(title="Single-cell clusters") + theme_pub(base_size=10) +
  theme(axis.text=element_blank(), axis.ticks=element_blank()) +
  guides(color=guide_legend(override.aes=list(size=2.5)))

## ---- Panel B: 火山图 -------------------------------------------------------
ng<-2000; lfc<-rnorm(ng,0,0.8); p<-runif(ng)
de<-sample(ng,120); lfc[de]<-lfc[de]+sample(c(-1,1),120,TRUE)*runif(120,1.5,3.5); p[de]<-p[de]*1e-4
vol<-data.frame(gene=sprintf("G%04d",1:ng), lfc=lfc, padj=p.adjust(p,"BH"))
vol$sig<-ifelse(vol$padj<0.05 & abs(vol$lfc)>1, ifelse(vol$lfc>0,"Up","Down"),"ns")
top<-vol[vol$sig!="ns",]; top<-top[order(top$padj),][1:8,]
pB <- ggplot(vol, aes(lfc, -log10(padj), color=sig)) +
  geom_point(size=0.7, alpha=0.6) +
  geom_vline(xintercept=c(-1,1), linetype=2, color="grey70") +
  geom_hline(yintercept=-log10(0.05), linetype=2, color="grey70") +
  ggrepel::geom_text_repel(data=top, aes(label=gene), size=2.4, max.overlaps=20) +
  scale_color_manual(values=c(Up="#B2182B",Down="#2166AC",ns="grey80"), name=NULL) +
  labs(x=expression(log[2]~FC), y=expression(-log[10]~FDR), title="Differential expression") +
  theme_pub(base_size=10)

## ---- Panel C: 顶 DEG 热图(行 z-score,RdBu)------------------------------
topg<-vol[order(vol$padj),"gene"][1:20]; ns<-12; grp<-rep(c("Ctrl","Case"),each=6)
H<-matrix(rnorm(20*ns), 20, ns, dimnames=list(topg, paste0(grp,1:ns)))
H[,grp=="Case"]<-H[,grp=="Case"]+rep(sample(c(-1.5,1.5),20,TRUE), times=1)  # 组间差异
Hz<-t(scale(t(H)))
hd<-expand.grid(gene=factor(topg,levels=rev(topg)), sample=factor(colnames(H),levels=colnames(H)))
hd$z<-as.vector(Hz)
pC <- ggplot(hd, aes(sample, gene, fill=z)) + geom_tile() +
  scale_fill_gradient2(low="#2166AC",mid="#F7F7F7",high="#B2182B",midpoint=0,name="z") +
  labs(x=NULL,y=NULL,title="Top DEG heatmap") + theme_pub(base_size=9) +
  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=6), axis.text.y=element_text(size=6))

## ---- Panel D: 关键基因效应森林(HR + 95%CI)------------------------------
fg<-top$gene[1:6]; hr<-exp(rnorm(6,0.2,0.4)); lo<-hr*exp(-0.3); hi<-hr*exp(0.3)
fd<-data.frame(gene=factor(fg,levels=rev(fg)), hr=hr, lo=lo, hi=hi)
pD <- ggplot(fd, aes(hr, gene)) +
  geom_vline(xintercept=1, linetype=2, color="grey60") +
  geom_errorbar(aes(xmin=lo,xmax=hi), orientation="y", width=0.2, linewidth=0.7, color=pal_pub(1)) +
  geom_point(size=2.6, color=pal_pub(1)) +
  labs(x="Hazard ratio (95% CI)", y=NULL, title="Key-gene association") +
  theme_pub(base_size=10)

## ---- 拼版(A/B/C/D 角标)+ 导出 ------------------------------------------
fig <- compose_panels(list(pA,pB,pC,pD), ncol=2, tag="A")
save_fig(fig, file.path(DAST,"composite_figure1"), width=NATURE_W2, height=6.6)
write.csv(vol[vol$sig!="ns",], file.path(DRES,"deg_significant.csv"), row.names=FALSE)
cat("[fig] assets/composite_figure1.{pdf,png} (4-panel A/B/C/D)\n")
sink(file.path(DRES,"sessionInfo.txt")); print(sessionInfo()); sink()
