# =============================================================================
# 510 · 单细胞代谢通路活性 (scMetabolism / AUCell-style metabolic scoring)
# -----------------------------------------------------------------------------
# 用秩和(AUCell/UCell 式)算法对每个细胞在多条代谢通路(KEGG/Reactome 代谢子集)
# 上打活性分,再按细胞类型聚合,刻画"哪类细胞代谢偏好哪条通路"。复现 scMetabolism
# 的核心产出:通路×细胞类型 活性热图 + 标志性 dotplot(点大小=活跃细胞比例,
# 颜色=平均活性)+ 关键差异通路的分布。
#
# Turnkey: Rscript 510_scmetabolism_pathway_activity.R   (合成 scRNA→results/+assets/)
#          换数据: --expr expr.csv --meta meta.csv --genesets pathways.gmt
# 复用 _framework/theme_pub.R;无条形图(heatmap+dotplot+violin)。
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

## ---- AUCell/UCell 式秩和打分(每细胞每通路 0-1 活性)-----------------------
score_sets <- function(expr, gsets, maxRank=1000){
  out <- sapply(gsets, function(sig){
    sig <- intersect(sig, rownames(expr)); ns <- length(sig)
    if(ns < 3) return(rep(NA_real_, ncol(expr)))
    apply(expr, 2, function(x){
      r <- rank(-x, ties.method="average"); r[r>maxRank] <- maxRank+1
      1 - (sum(r[match(sig, rownames(expr))]) - ns*(ns+1)/2)/(ns*maxRank)
    })
  })
  t(out)   # 通路 x 细胞
}

## ---- 1. 合成 scRNA(3 细胞类型,各有代谢偏好)+ 代谢通路基因集 -------------
ng<-200; nc<-360; genes<-sprintf("G%03d",1:ng); cells<-sprintf("C%03d",1:nc)
ctype <- rep(c("Tumor","CAF","Immune"), each=nc/3)
pathways <- list(
  Glycolysis        = genes[1:18],
  OXPHOS            = genes[19:38],
  Fatty_acid_oxid   = genes[39:56],
  Glutaminolysis    = genes[57:72],
  Nucleotide_synth  = genes[73:90],
  Pentose_phosphate = genes[91:106],
  One_carbon        = genes[107:122],
  Lipid_synth       = genes[123:140])
E <- matrix(rpois(ng*nc, 4), ng, nc, dimnames=list(genes,cells))
boost <- function(g, cells_idx, lam) E[g, cells_idx] <<- E[g, cells_idx] + rpois(length(g)*length(cells_idx), lam)
boost(pathways$Glycolysis, which(ctype=="Tumor"), 7)        # 肿瘤:糖酵解(Warburg)
boost(pathways$Nucleotide_synth, which(ctype=="Tumor"), 5)
boost(pathways$OXPHOS, which(ctype=="Immune"), 7)           # 免疫:OXPHOS
boost(pathways$Fatty_acid_oxid, which(ctype=="CAF"), 7)     # CAF:脂肪酸氧化
boost(pathways$Lipid_synth, which(ctype=="CAF"), 5)
write.csv(round(E,2), file.path(DDAT,"expr.csv"))
write.csv(data.frame(cell=cells, celltype=ctype), file.path(DDAT,"meta.csv"), row.names=FALSE)
cat(sprintf("[gen] synthetic scRNA %dx%d, %d metabolic pathways, 3 cell types (demo only)\n", ng,nc,length(pathways)))

## ---- 2. 打分 + 按细胞类型聚合 ---------------------------------------------
S <- score_sets(E, pathways)                                # 通路 x 细胞
types <- unique(ctype)
mean_act <- sapply(types, function(t) rowMeans(S[, ctype==t, drop=FALSE]))   # 通路 x 类型
frac_act <- sapply(types, function(t){                      # 活跃比例(>该通路中位)
  thr <- apply(S, 1, median); rowMeans(S[, ctype==t, drop=FALSE] > thr)
})
write.csv(data.frame(pathway=rownames(mean_act), mean_act, check.names=FALSE),
          file.path(DRES,"pathway_activity_by_celltype.csv"), row.names=FALSE)
# 差异通路(Kruskal-Wallis across cell types)
kw <- apply(S, 1, function(v) kruskal.test(v, factor(ctype))$p.value)
write.csv(data.frame(pathway=names(kw), KW_p=kw), file.path(DRES,"pathway_kruskal.csv"), row.names=FALSE)
top_path <- names(sort(kw))[1]
cat(sprintf("[score] most cell-type-variable pathway = %s (KW p=%.2g)\n", top_path, min(kw)))

## ---- 3. 出图(heatmap + dotplot + violin;无条形图)-----------------------
# Fig1: 通路 x 细胞类型 平均活性热图(行 z-score,RdBu 发散)
zm <- t(scale(t(mean_act)))                                 # 行内 z-score
hd <- expand.grid(pathway=rownames(zm), celltype=colnames(zm))
hd$z <- as.vector(zm)
hd$pathway <- factor(hd$pathway, levels=rev(rownames(zm)))
p1 <- ggplot(hd, aes(celltype, pathway, fill=z)) +
  geom_tile(color="white", linewidth=0.6) +
  scale_fill_gradient2(low="#2166AC", mid="#F7F7F7", high="#B2182B", midpoint=0, name="Activity\n(row z)") +
  labs(x=NULL, y=NULL, title="Metabolic pathway activity (scMetabolism-style)") +
  theme_pub(base_size=11) + theme(axis.text.x=element_text(angle=0))
save_fig(p1, file.path(DAST,"pathway_heatmap"), width=5.2, height=4.2)

# Fig2: scMetabolism 标志性 dotplot(点大小=活跃比例,颜色=平均活性)
dd <- expand.grid(pathway=rownames(mean_act), celltype=colnames(mean_act))
dd$mean <- as.vector(mean_act); dd$frac <- as.vector(frac_act)
dd$pathway <- factor(dd$pathway, levels=rev(rownames(mean_act)))
p2 <- ggplot(dd, aes(celltype, pathway)) +
  geom_point(aes(size=frac, color=mean)) +
  scale_size(range=c(1,8), name="Frac active") +
  scale_color_viridis_c(option="D", name="Mean activity") +
  labs(x=NULL, y=NULL, title="Metabolic activity dotplot") +
  theme_pub(base_size=11)
save_fig(p2, file.path(DAST,"pathway_dotplot"), width=5.4, height=4.2)

# Fig3: 最差异通路的活性分布(violin+box,按细胞类型)
vd <- data.frame(celltype=factor(ctype, levels=types), score=S[top_path,])
p3 <- ggplot(vd, aes(celltype, score, fill=celltype)) +
  geom_violin(alpha=0.55, color=NA, trim=FALSE) +
  geom_boxplot(width=0.16, outlier.shape=NA, alpha=0.9) +
  scale_fill_manual(values=pal_pub(3,"npg"), guide="none") +
  labs(x=NULL, y=sprintf("%s activity", top_path),
       title=sprintf("%s across cell types (KW p=%.2g)", top_path, min(kw))) +
  theme_pub(base_size=11)
save_fig(p3, file.path(DAST,"top_pathway_violin"), width=4.4, height=3.6)

cat("[fig] assets/: pathway_heatmap, pathway_dotplot, top_pathway_violin (.pdf+.png)\n")
sink(file.path(DRES,"sessionInfo.txt")); print(sessionInfo()); sink()
