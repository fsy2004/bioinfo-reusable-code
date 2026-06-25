# =============================================================================
# 515 · 弦图 (chord diagram) — 关系/流量的环形可视化
# -----------------------------------------------------------------------------
# 用 circlize::chordDiagram 把一个"来源 × 去向"强度矩阵画成环形弦图,直观展示
# 细胞-细胞通讯强度、配体-受体流量、簇间共享基因、状态转移等有向关系。带方向箭头、
# 期刊配色扇区。比堆叠条形/邻接表更能凸显主导通路。
#
# 注:circlize 为 base 绘图,本模块自带 PDF+PNG 双设备导出(不走 ggplot 的 save_fig)。
# Turnkey: Rscript 515_chord_diagram.R   (合成通讯矩阵→results/+assets/)
#          换数据: --input matrix.csv  (方阵/矩阵,行=来源,列=去向)
# 复用 _framework/theme_pub.R 的期刊配色;无条形图。
# =============================================================================
suppressWarnings(suppressMessages({ library(circlize) }))

.this <- tryCatch({ a <- commandArgs(FALSE); m <- grep("^--file=", a)
  if(length(m)) dirname(normalizePath(sub("^--file=","",a[m[1]]))) else getwd() },
  error=function(e) getwd())
.fw <- NULL; .p <- .this
for(i in 1:6){ cand <- file.path(.p,"_framework","theme_pub.R"); if(file.exists(cand)){ .fw <- cand; break }; .p <- dirname(.p) }
if(!is.null(.fw)) source(.fw) else { pal_pub <- function(n=NULL,name="npg") grDevices::rainbow(ifelse(is.null(n),6,n)) }
args <- if(exists("bio_args")) bio_args(list(input=NULL)) else list(input=NULL)

set.seed(42)
DIR <- .this
DDAT <- file.path(DIR,"example_data"); DRES <- file.path(DIR,"results"); DAST <- file.path(DIR,"assets")
for(d in c(DDAT,DRES,DAST)) dir.create(d, showWarnings=FALSE, recursive=TRUE)

## ---- 1. 读/生成 来源×去向 强度矩阵(细胞-细胞通讯)------------------------
fM <- file.path(DDAT,"interaction_matrix.csv")
if(is.null(args$input)){
  if(!file.exists(fM)){
    ct <- c("Tumor","CAF","Endothelial","Tcell","Macrophage")
    M <- matrix(round(runif(25,0,3),1), 5, 5, dimnames=list(ct,ct))
    M["Tumor","CAF"]<-9; M["CAF","Tumor"]<-6; M["Macrophage","Tcell"]<-7  # 几条主导通讯
    M["Tumor","Macrophage"]<-5; M["Endothelial","Tcell"]<-4; diag(M)<-0
    write.csv(M, fM)
    cat("[gen] synthetic 5x5 cell-cell interaction matrix (demo only)\n")
  }
  M <- as.matrix(read.csv(fM, row.names=1, check.names=FALSE))
} else M <- as.matrix(read.csv(args$input, row.names=1, check.names=FALSE))

## ---- 2. 配色 + 绘制函数(可被 PDF/PNG 两次调用)--------------------------
sectors <- union(rownames(M), colnames(M))
grid.col <- setNames(pal_pub(length(sectors),"okabe_ito"), sectors)
draw_chord <- function(){
  circos.clear()
  # canvas 扩到 ±1.25:圈占 80% 画布,四周留白给旋转的长扇区标签(防 Endothelial 被裁)
  circos.par(gap.after=4, start.degree=90, canvas.xlim=c(-1.25,1.25), canvas.ylim=c(-1.25,1.25))
  chordDiagram(M, grid.col=grid.col, directional=1,
               direction.type=c("diffHeight","arrows"), link.arr.type="big.arrow",
               diffHeight=mm_h(2), annotationTrack="grid",
               preAllocateTracks=list(track.height=0.08))
  circos.track(track.index=1, panel.fun=function(x,y){
    circos.text(CELL_META$xcenter, CELL_META$ylim[1]+mm_y(3), CELL_META$sector.index,
                facing="clockwise", niceFacing=TRUE, adj=c(0,0.5), cex=0.8, font=2)
  }, bg.border=NA)
  title("Cell-cell communication chord diagram", cex.main=1.1)
  circos.clear()
}
## ---- 3. 双设备导出(矢量 PDF + 300dpi PNG)-------------------------------
pdf_ok <- tryCatch({ grDevices::cairo_pdf(file.path(DAST,"chord_diagram.pdf"), width=6, height=6)
  draw_chord(); dev.off(); TRUE }, error=function(e){ try(dev.off(),silent=TRUE); FALSE })
png(file.path(DAST,"chord_diagram.png"), width=6, height=6, units="in", res=300)
draw_chord(); dev.off()

## ---- 4. 落盘:主导通讯对 ---------------------------------------------------
idx <- which(M>0, arr.ind=TRUE)
flow <- data.frame(source=rownames(M)[idx[,1]], target=colnames(M)[idx[,2]],
                   strength=M[idx]); flow <- flow[order(-flow$strength),]
write.csv(flow, file.path(DRES,"interaction_flows.csv"), row.names=FALSE)
cat(sprintf("[chord] top flow: %s -> %s (%.1f); PDF=%s\n", flow$source[1], flow$target[1], flow$strength[1], pdf_ok))
cat("[fig] assets/chord_diagram.{pdf,png}\n")
sink(file.path(DRES,"sessionInfo.txt")); print(sessionInfo()); sink()
