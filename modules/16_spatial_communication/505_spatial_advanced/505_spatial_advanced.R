# =============================================================================
# 505 · Spatial advanced toolkit: deconvolution + niche + co-localization
# -----------------------------------------------------------------------------
# 空间高级套件(补 027/050 基础空间):
#   (1) RCTD(spacexr)参考型反卷积 -> 每 spot 细胞类型比例
#   (2) NMF(RcppML)分解比例矩阵 -> 空间生态位(niche factor)
#   (3) CellDegree(KNN-6 异型邻接,base 实现)-> 细胞类型物理共定位
#   (4) MISTy(mistyR)多视图空间依赖(可选,失败不致命)
# 思路来源:PDAC 时空范文(RCTD+MISTy+semla 邻接+NMF 多方法叠证锁定空间生态位)。
# Turnkey: Rscript 505_spatial_advanced.R  (默认读 example_data/, 写 results/+assets/)
# 复用 _framework/theme_pub.R;代码注释中文,图中文字英文。
# =============================================================================
suppressWarnings(suppressMessages({
  library(spacexr); library(RcppML); library(Matrix); library(ggplot2)
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

## ---- 1. 生成示例:单细胞参考 + 空间切片(spot=细胞类型混合,有空间结构) --------
fD <- file.path(DDAT,"spatial_demo.rds")
if(!file.exists(fD)){
  ng <- 200; nct <- 3; cts <- c("Tumor","CAF","Immune")
  genes <- sprintf("G%03d",1:ng)
  # 每类 marker(前 60 基因分 3 段各 20 为类特异高表达)
  mk <- list(Tumor=1:20, CAF=21:40, Immune=41:60)
  prof <- sapply(cts, function(ct){ v <- rep(2, ng); v[mk[[ct]]] <- 20; v })  # 类表达谱
  rownames(prof) <- genes
  # 单细胞参考:每类 200 细胞
  ncell <- 600; cell_ct <- rep(cts, each=ncell/nct)
  sc <- sapply(seq_len(ncell), function(i) rpois(ng, prof[,cell_ct[i]]))
  rownames(sc) <- genes; colnames(sc) <- sprintf("cell%04d",1:ncell)
  # 空间切片:18x18 grid = 324 spot;按区域设细胞类型比例(左=Tumor 中=CAF 右=Immune 渐变)
  gx <- 18; coords <- expand.grid(x=1:gx, y=1:gx); rownames(coords) <- sprintf("spot%03d",1:nrow(coords))
  fx <- coords$x/gx
  prop <- cbind(Tumor=pmax(0,1-1.4*fx)+0.1, CAF=1-abs(fx-0.5)*2+0.1, Immune=pmax(0,1.4*fx-0.4)+0.1)
  prop <- prop/rowSums(prop); rownames(prop) <- rownames(coords)
  # 每 spot ~10 细胞混合,counts = 比例加权类谱 * 深度
  st <- sapply(seq_len(nrow(coords)), function(i) rpois(ng, (prof %*% prop[i,]) * 5))
  rownames(st) <- genes; colnames(st) <- rownames(coords)
  saveRDS(list(sc=sc, cell_ct=cell_ct, st=st, coords=coords, true_prop=prop), fD)
  cat("[gen] synthetic: 600-cell ref (3 types) + 324-spot slice (spatial gradient)\n")
}
D <- readRDS(fD)

## ---- 2. RCTD 参考型反卷积 -------------------------------------------------
ref <- Reference(D$sc, factor(setNames(D$cell_ct, colnames(D$sc))))
puck <- SpatialRNA(D$coords, as(D$st, "dgCMatrix"))
myRCTD <- create.RCTD(puck, ref, max_cores = 2, CELL_MIN_INSTANCE = 20)
myRCTD <- run.RCTD(myRCTD, doublet_mode = "full")
W <- as.matrix(myRCTD@results$weights); W <- sweep(W, 1, rowSums(W), "/")   # 归一化为比例
W <- W[, , drop=FALSE]
write.csv(W, file.path(DRES,"RCTD_weights.csv"))
# 反卷积准确性(与真值相关)
common <- intersect(rownames(W), rownames(D$true_prop))
acc <- mean(sapply(colnames(W), function(ct) cor(W[common,ct], D$true_prop[common,ct])))
cat(sprintf("[RCTD] %d spots deconvolved; mean cor(estimated, truth)=%.2f\n", nrow(W), acc))

## ---- 3. NMF 空间生态位(RcppML) -------------------------------------------
k <- min(3, ncol(W))
nm <- RcppML::nmf(W, k = k, seed = 42)
Hfac <- nm$w                                # spot x factor (RcppML 返回 list: w/d/h)
rownames(Hfac) <- rownames(W); colnames(Hfac) <- paste0("Niche", seq_len(k))
niche <- colnames(Hfac)[max.col(Hfac)]
write.csv(data.frame(spot=rownames(W), Hfac, niche=niche), file.path(DRES,"NMF_niches.csv"), row.names=FALSE)
cat(sprintf("[NMF] %d spatial niches from cell-type composition\n", k))

## ---- 4. CellDegree:KNN-6 异型邻接(base 实现,类比 semla) -------------------
co <- D$coords[rownames(W),]
dom <- colnames(W)[max.col(W)]
dmat <- as.matrix(dist(co))
hetero <- sapply(seq_len(nrow(co)), function(i){
  nn <- order(dmat[i,])[2:7]                 # 6 近邻
  mean(dom[nn] != dom[i])                     # 异型邻居比例
})
cd <- data.frame(spot=rownames(W), dominant=dom, hetero_degree=hetero)
write.csv(cd, file.path(DRES,"CellDegree.csv"), row.names=FALSE)
cat(sprintf("[CellDegree] mean heterotypic neighbor fraction=%.2f (high=interface niche)\n", mean(hetero)))

## ---- 5. MISTy 多视图(可选,失败不致命) ------------------------------------
misty_ok <- FALSE
try({
  suppressMessages(library(mistyR)); library(dplyr)
  expr <- as.data.frame(W);
  views <- create_initial_view(expr) %>% add_paraview(co, l = 3)
  run_misty(views, results.folder = file.path(DRES,"misty"))
  misty_ok <- TRUE; cat("[MISTy] multi-view spatial dependency -> results/misty/\n")
}, silent = TRUE)
if(!misty_ok) cat("[MISTy] skipped (optional; see README to enable)\n")

## ---- 6. 出图 --------------------------------------------------------------
pdat <- data.frame(co, dominant=dom, niche=niche, hetero=hetero, Tumor=W[,"Tumor"])
# (a) 主导细胞类型空间图
p1 <- ggplot(pdat, aes(x, y, color=dominant)) + geom_point(size=2.4) +
  scale_color_manual(values=pal_pub(3), name="Dominant type") + coord_fixed() +
  labs(title="RCTD: dominant cell type per spot") + theme_pub(base_size=11)
save_fig(p1, file.path(DAST,"rctd_dominant"), width=5.5, height=4.5)
# (b) NMF 生态位空间图
p2 <- ggplot(pdat, aes(x, y, color=niche)) + geom_point(size=2.4) +
  scale_color_manual(values=pal_pub(k), name="Niche") + coord_fixed() +
  labs(title="NMF spatial niches") + theme_pub(base_size=11)
save_fig(p2, file.path(DAST,"nmf_niche"), width=5.5, height=4.5)
# (c) 异型邻接度(界面)空间图
p3 <- ggplot(pdat, aes(x, y, color=hetero)) + geom_point(size=2.4) +
  scale_color_gradientn(colours=pal_pub(5), name="Heterotypic\nneighbors") + coord_fixed() +
  labs(title="CellDegree: cell-type interface") + theme_pub(base_size=11)
save_fig(p3, file.path(DAST,"celldegree_interface"), width=5.5, height=4.5)
cat("[fig] assets/rctd_dominant, nmf_niche, celldegree_interface (.pdf/.png)\n")
