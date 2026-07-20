# =============================================================================
# 504 · hdWGCNA single-cell co-expression network
# -----------------------------------------------------------------------------
# 单细胞共表达网络(metacell 聚合抗稀疏),补 054(bulk WGCNA)的单细胞版。
# 思路来源:dry-AMD 范文用 hdWGCNA 在 metacell 上建网找疾病模块 hub。
# 流程:Seurat 标准 -> SetupForWGCNA -> MetacellsByGroups -> SetDatExpr ->
#       TestSoftPowers -> ConstructNetwork -> ModuleEigengenes/Connectivity -> hub
# Turnkey: Rscript 504_hdwgcna_single_cell.R  (默认读 example_data/, 写 results/+assets/)
# 复用 _framework/theme_pub.R;代码注释中文,图中文字英文。
# =============================================================================
suppressWarnings(suppressMessages({
  library(Seurat); library(hdWGCNA); library(WGCNA); library(igraph)
  library(patchwork); library(ggplot2); library(dplyr)
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

set.seed(42); enableWGCNAThreads(nThreads = 2)
DIR <- .this; DDAT <- file.path(DIR,"example_data"); DRES <- file.path(DIR,"results"); DAST <- file.path(DIR,"assets")
for(d in c(DDAT,DRES,DAST)) dir.create(d, showWarnings=FALSE, recursive=TRUE)

## ---- 1. 读 / 生成示例单细胞数据(植入 3 个共表达模块) -----------------------
fSeu <- file.path(DDAT,"sc_counts.rds")
if(!file.exists(fSeu)){
  nc <- 900; ng <- 250
  genes <- sprintf("G%03d",1:ng); cells <- sprintf("C%04d",1:nc)
  ct <- rep(c("A","B","C"), length.out=nc)
  counts <- matrix(rpois(ng*nc, lambda=1.5), ng, nc, dimnames=list(genes,cells))   # 稀疏底
  # 3 个共表达模块(各 30 基因),分别在 A/B/C 富集且模块内共变
  mods <- list(M1=1:30, M2=31:60, M3=61:90); cts <- c("A","B","C")
  for(k in 1:3){
    lat <- rnorm(nc, 0, 1) + ifelse(ct==cts[k], 3, 0)            # 模块潜变量(对应细胞型高)
    lat <- pmax(lat, 0)
    for(g in mods[[k]]) counts[g,] <- counts[g,] + rpois(nc, lambda=lat*1.5)  # 模块基因随潜变量共变
  }
  saveRDS(counts, fSeu); cat("[gen] synthetic scRNA: 900 cells x 250 genes, 3 co-expression modules\n")
}
counts <- readRDS(fSeu)
ct <- rep(c("A","B","C"), length.out=ncol(counts))

## ---- 2. Seurat 标准流程 ---------------------------------------------------
seu <- CreateSeuratObject(counts = counts)
seu$cell_type <- ct
seu <- NormalizeData(seu, verbose=FALSE)
seu <- FindVariableFeatures(seu, nfeatures=200, verbose=FALSE)
seu <- ScaleData(seu, features=rownames(seu), verbose=FALSE)
seu <- RunPCA(seu, npcs=20, verbose=FALSE)
seu <- RunUMAP(seu, dims=1:15, verbose=FALSE)

## ---- 3. hdWGCNA 流程 ------------------------------------------------------
seu <- SetupForWGCNA(seu, gene_select="fraction", fraction=0.05, wgcna_name="scWGCNA")
seu <- MetacellsByGroups(seu, group.by="cell_type", k=20, max_shared=10,
                         ident.group="cell_type", reduction="pca", min_cells=30)
seu <- NormalizeMetacells(seu)
seu <- SetDatExpr(seu, group_name=c("A","B","C"), group.by="cell_type", use_metacells=TRUE)
seu <- TestSoftPowers(seu, networkType="signed")
# 自动选 power 在合成/小数据上可能失败(scale-free fit 不达 0.8 -> Inf),显式取诊断表里 fit 最高的 power
.pt <- GetPowerTable(seu); .sp <- .pt$Power[which.max(.pt$SFT.R.sq)]
if(!is.finite(.sp) || is.na(.sp)) .sp <- 8
cat(sprintf("[hdWGCNA] using soft_power = %d (max scale-free fit)\n", .sp))
seu <- ConstructNetwork(seu, soft_power=.sp, tom_name="sc_demo", overwrite_tom=TRUE,
                        networkType="signed", minModuleSize=15, mergeCutHeight=0.2)
seu <- ModuleEigengenes(seu)
seu <- ModuleConnectivity(seu)

## ---- 4. 结果表 + hub ------------------------------------------------------
modules <- GetModules(seu)
write.csv(modules[, c("gene_name","module","color")], file.path(DRES,"modules.csv"), row.names=FALSE)
nmod <- length(setdiff(unique(modules$module), "grey"))
hubs <- tryCatch(GetHubGenes(seu, n_hubs=10), error=function(e) NULL)
if(!is.null(hubs)) write.csv(hubs, file.path(DRES,"hub_genes.csv"), row.names=FALSE)
cat(sprintf("[done] %d co-expression modules (excl. grey); hub genes -> results/\n", nmod))

## ---- 5. 出图 --------------------------------------------------------------
# (a) soft power 诊断
p_sp <- PlotSoftPowers(seu); pw <- patchwork::wrap_plots(p_sp, ncol=2)
save_fig(pw, file.path(DAST,"soft_power"), width=8, height=6)
# (b) 模块树状图(base graphics -> png/pdf 设备捕获)
png(file.path(DAST,"dendrogram.png"), width=2000, height=1200, res=300)
PlotDendrogram(seu, main="hdWGCNA dendrogram"); dev.off()
pdf(file.path(DAST,"dendrogram.pdf"), width=7, height=4.2)
PlotDendrogram(seu, main="hdWGCNA dendrogram"); dev.off()
# (c) 模块 eigengene 在 UMAP 上的特征图
p_mf <- tryCatch(ModuleFeaturePlot(seu, features="hMEs", order=TRUE), error=function(e) NULL)
if(!is.null(p_mf)) save_fig(patchwork::wrap_plots(p_mf, ncol=2), file.path(DAST,"module_featureplot"), width=8, height=6)
cat("[fig] assets/soft_power, dendrogram, module_featureplot (.pdf/.png)\n")
