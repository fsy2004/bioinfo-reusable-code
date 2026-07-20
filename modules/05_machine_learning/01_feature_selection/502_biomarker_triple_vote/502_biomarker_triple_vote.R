# =============================================================================
# 502 · Biomarker triple-vote feature selection
# -----------------------------------------------------------------------------
# 三法投票:网络拓扑(igraph 度中心性 top-k) × 相关中心性(与表型相关 top-k) ×
# 机器学习(Boruta Confirmed),只保留 >=2 法共同选中的基因 -> 高可信 biomarker 候选,
# 抗 single-signature 假阳性。思路来源:NETs 多组学范文(CytoHubba×相关×Boruta 交集)。
# Turnkey: Rscript 502_biomarker_triple_vote.R  (默认读 example_data/, 写 results/+assets/)
#          换数据: --input expr.csv --group group.csv --candidates candidate_genes.csv
# 复用 _framework/theme_pub.R;代码注释中文,图中文字英文。
# =============================================================================
suppressWarnings(suppressMessages({
  library(igraph); library(Boruta); library(Hmisc); library(ggplot2)
}))

## ---- 定位脚本目录 + 载入框架(向上搜 _framework) ----------------------------
.this <- tryCatch({ a <- commandArgs(FALSE); m <- grep("^--file=", a)
  if(length(m)) dirname(normalizePath(sub("^--file=","",a[m[1]]))) else getwd() },
  error=function(e) getwd())
.fw <- NULL; .p <- .this
for(i in 1:6){ cand_fw <- file.path(.p,"_framework","theme_pub.R"); if(file.exists(cand_fw)){ .fw <- cand_fw; break }; .p <- dirname(.p) }
if(!is.null(.fw)){ source(.fw) } else {
  theme_pub <- function(base_size=11, ...) theme_bw(base_size=base_size)
  pal_pub   <- function(n=NULL, name="npg") scales::hue_pal()(ifelse(is.null(n),6,n))
  save_fig  <- function(plot, file, width=7, height=6, dpi=300){
    ggsave(paste0(file,".pdf"), plot, width=width, height=height)
    ggsave(paste0(file,".png"), plot, width=width, height=height, dpi=dpi) }
}

set.seed(42)
TOPK <- 12   # 拓扑/相关各取 top-K 候选(投票前的单法入选数)
DIR  <- .this
DDAT <- file.path(DIR,"example_data"); DRES <- file.path(DIR,"results"); DAST <- file.path(DIR,"assets")
for(d in c(DDAT,DRES,DAST)) dir.create(d, showWarnings=FALSE, recursive=TRUE)

## ---- 1. 读 / 生成示例数据 -------------------------------------------------
# 合成:8 个真信号基因(tre 组强上调 + 彼此共表达成模块);其余为独立噪声。
# 设计使三法各有区分度且真信号被三法共识。
fE <- file.path(DDAT,"expr.csv"); fG <- file.path(DDAT,"group.csv"); fC <- file.path(DDAT,"candidate_genes.csv")
if(!file.exists(fE)){
  ng <- 60; ns <- 40
  genes <- sprintf("GENE%02d",1:ng); samp <- sprintf("S%02d",1:ns); grp <- rep(c("con","tre"), each=ns/2)
  M <- matrix(rnorm(ng*ns, 8, 1.2), ng, ns, dimnames=list(genes,samp))   # 独立噪声底
  sig <- genes[1:8]
  M[sig, grp=="tre"] <- M[sig, grp=="tre"] + 3.0                          # 真信号:强上调(表型相关)
  modk <- as.numeric(scale(rnorm(ns)))                                    # 共表达潜变量(只作用于信号)
  for(g in sig) M[g,] <- M[g,] + 1.6*modk                                 # 信号基因彼此共表达成模块(拓扑 hub)
  write.csv(round(M,3), fE); write.csv(data.frame(sample=samp, group=grp), fG, row.names=FALSE)
  write.csv(data.frame(gene=genes), fC, row.names=FALSE)
  cat("[gen] synthetic example_data: 60 genes x 40 samples, 8 true co-expressed signals\n")
}
expr <- as.matrix(read.csv(fE, row.names=1, check.names=FALSE))
grp  <- read.csv(fG); cand <- intersect(as.character(read.csv(fC)$gene), rownames(expr))
y    <- factor(grp$group[match(colnames(expr), grp$sample)]); yb <- as.integer(y) - 1

## ---- 2. 三法(各自有区分度) ------------------------------------------------
# (1) 网络拓扑:候选相关网络的【度中心性 top-K】hub(类比 CytoHubba)
cc <- cor(t(expr[cand,])); diag(cc) <- 0
gr <- graph_from_adjacency_matrix(abs(cc) > 0.5, mode="undirected", diag=FALSE)
deg <- degree(gr); topo_hub <- names(sort(deg, decreasing=TRUE))[seq_len(min(TOPK, sum(deg>0)))]
# (2) 相关中心性:与表型 point-biserial 相关最强的 top-K
pcor <- sapply(cand, function(x) abs(cor(expr[x,], yb)))
corr_hub <- names(sort(pcor, decreasing=TRUE))[seq_len(TOPK)]
# (3) 机器学习:Boruta(RF wrapper),TentativeRoughFix 兜底未决项
bdat <- as.data.frame(t(expr[cand,])); bdat$.y <- y
bor  <- TentativeRoughFix(Boruta(.y ~ ., data=bdat, doTrace=0, maxRuns=300))
ml_hub <- getSelectedAttributes(bor)

## ---- 3. 投票交集(>=2 法) --------------------------------------------------
allg <- unique(c(topo_hub, corr_hub, ml_hub))
vote <- data.frame(gene=allg,
                   topology    = as.integer(allg %in% topo_hub),
                   correlation = as.integer(allg %in% corr_hub),
                   ML_Boruta   = as.integer(allg %in% ml_hub))
vote$votes <- rowSums(vote[,2:4]); vote <- vote[order(-vote$votes, vote$gene),]
consensus  <- vote$gene[vote$votes >= 2]
write.csv(vote, file.path(DRES,"vote_table.csv"), row.names=FALSE)
write.csv(data.frame(consensus_biomarker=consensus), file.path(DRES,"consensus_biomarkers.csv"), row.names=FALSE)
cat(sprintf("[methods] topology=%d  correlation=%d  Boruta=%d\n", length(topo_hub), length(corr_hub), length(ml_hub)))
cat(sprintf("[done] %d consensus biomarkers (>=2 votes): %s\n", length(consensus), paste(consensus, collapse=", ")))

## ---- 4. 出图(矢量 PDF + 300dpi PNG) ---------------------------------------
vm <- vote[vote$votes >= 1,]
vlong <- data.frame(gene=rep(vm$gene,3),
                    method=factor(rep(c("Topology","Correlation","ML (Boruta)"), each=nrow(vm)),
                                  levels=c("Topology","Correlation","ML (Boruta)")),
                    value=c(vm$topology, vm$correlation, vm$ML_Boruta))
vlong$gene <- factor(vlong$gene, levels=rev(vm$gene))
p1 <- ggplot(vlong, aes(method, gene, fill=factor(value))) +
  geom_tile(color="white", linewidth=0.4) +
  scale_fill_manual(values=c("0"="grey90","1"=pal_pub(2)[1]), name="selected", labels=c("no","yes")) +
  labs(x="Method", y="Gene", title="Biomarker triple-vote matrix") + theme_pub(base_size=10)
save_fig(p1, file.path(DAST,"vote_matrix"), width=5, height=max(4, nrow(vm)*0.22))
if(length(consensus)){
  dd <- data.frame(gene=consensus, cor=pcor[consensus])
  dd$gene <- factor(dd$gene, levels=dd$gene[order(dd$cor)])
  p2 <- ggplot(dd, aes(cor, gene)) +                          # lollipop(顶刊优于条形)
    geom_segment(aes(x=0, xend=cor, yend=gene, colour=cor), linewidth=1.1) +
    geom_point(aes(colour=cor), size=4) +
    scale_colour_gradientn(colours=pal_pub(5), guide="none") +
    scale_x_continuous(expand=expansion(mult=c(0,0.08))) +
    labs(x="|correlation with phenotype|", y=NULL, title="Consensus biomarkers (>=2 votes)") +
    theme_pub(base_size=11)
  save_fig(p2, file.path(DAST,"consensus_bar"), width=6, height=4)
}
cat("[fig] assets/vote_matrix.{pdf,png}, consensus_bar.{pdf,png}\n")
