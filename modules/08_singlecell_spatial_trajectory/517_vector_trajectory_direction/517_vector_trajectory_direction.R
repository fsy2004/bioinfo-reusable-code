# =============================================================================
# 517 · VECTOR 轨迹分化方向 (differentiation direction as a vector field)
# -----------------------------------------------------------------------------
# 复现 VECTOR(Zhang et al.)核心思想:不依赖 RNA velocity,仅用表达数据在 2D 嵌入上
# 推断分化【方向】。原理:细胞"潜能/干性"可由表达广度近似(CytoTRACE 假设:表达基因
# 数越多越接近祖细胞);把嵌入网格化,箭头沿"潜能下降"梯度方向 = 分化方向,叠成矢量场。
#
# 产出本身即高级图(quiver/矢量场 + 潜能着色 UMAP),不是条形图。
# ★实现的是方法核心(turnkey 免装 niche 包);接真数据时潜能可换 CytoTRACE2(模块082)。
# Turnkey: Rscript 517_vector_trajectory_direction.R   (合成分叉轨迹→results/+assets/)
#          换数据: --embedding umap.csv --expr expr.csv
# 复用 _framework/theme_pub.R;无条形图。
# =============================================================================
suppressWarnings(suppressMessages({ library(ggplot2) }))

.this <- tryCatch({ a <- commandArgs(FALSE); m <- grep("^--file=", a)
  if(length(m)) dirname(normalizePath(sub("^--file=","",a[m[1]]))) else getwd() },
  error=function(e) getwd())
.fw <- NULL; .p <- .this
for(i in 1:6){ cand <- file.path(.p,"_framework","theme_pub.R"); if(file.exists(cand)){ .fw <- cand; break }; .p <- dirname(.p) }
if(!is.null(.fw)) source(.fw) else stop("需要 _framework/theme_pub.R")
args <- bio_args(list(embedding=NULL, expr=NULL, ngrid=14))

set.seed(42)
DIR <- .this
DDAT <- file.path(DIR,"example_data"); DRES <- file.path(DIR,"results"); DAST <- file.path(DIR,"assets")
for(d in c(DDAT,DRES,DAST)) dir.create(d, showWarnings=FALSE, recursive=TRUE)

## ---- 1. 合成:分叉轨迹(1 根→2 分支)+ 表达(潜能=表达广度)----------------
if(is.null(args$embedding)){
  n1<-260; n2<-200; n3<-200
  # 主干 + 两分支的 2D 坐标
  t1<-seq(0,1,length.out=n1); trunk<-cbind(t1*4-2, t1*2-3 + rnorm(n1,0,0.18))
  t2<-seq(0,1,length.out=n2); b1<-cbind(2+t2*2, -1+t2*2.5 + rnorm(n2,0,0.18))
  t3<-seq(0,1,length.out=n3); b2<-cbind(2+t3*2, -1-t3*2.5 + rnorm(n3,0,0.18))
  emb<-rbind(trunk,b1,b2); colnames(emb)<-c("UMAP1","UMAP2")
  potency<-c(1-t1*0.5, 0.5-t2*0.5, 0.5-t3*0.5)            # 根部潜能高→分支末端低
  potency<-pmax(0.02, potency + rnorm(length(potency),0,0.03))
  # 表达广度:表达基因数 ∝ 潜能(CytoTRACE 假设),据此【反推】潜能(不直接用真值)
  ng<-400; ncell<-nrow(emb)
  ngene_expr<-round(60 + potency*300)                      # 每细胞表达基因数
  E<-matrix(0, ng, ncell)
  for(j in 1:ncell){ idx<-sample(ng, min(ng,ngene_expr[j])); E[idx,j]<-rpois(length(idx),3)+1 }
  write.csv(round(emb,3), file.path(DDAT,"embedding.csv"), row.names=FALSE)
  saveRDS(E, file.path(DDAT,"expr.rds"))
  cat(sprintf("[gen] synthetic branching trajectory: %d cells, root->2 branches (demo only)\n", ncell))
} else { emb<-as.matrix(read.csv(args$embedding)); E<-as.matrix(read.csv(args$expr, row.names=1)) }

## ---- 2. VECTOR 核心:由表达广度估潜能(越高越祖细胞)-----------------------
ngrid <- as.integer(args$ngrid)
potency_est <- colSums(E>0); potency_est <- (potency_est-min(potency_est))/(diff(range(potency_est)))
df <- data.frame(UMAP1=emb[,1], UMAP2=emb[,2], potency=potency_est)

## ---- 3. 网格化 + 梯度 → 矢量场(箭头指向潜能下降=分化方向)-----------------
gx<-seq(min(df$UMAP1),max(df$UMAP1),length.out=ngrid)
gy<-seq(min(df$UMAP2),max(df$UMAP2),length.out=ngrid)
cellx<-findInterval(df$UMAP1,gx); celly<-findInterval(df$UMAP2,gy)
grid_pot<-matrix(NA,ngrid,ngrid); gc<-matrix(0,ngrid,ngrid)
for(i in 1:nrow(df)){ a<-max(1,cellx[i]); b<-max(1,celly[i])
  grid_pot[a,b]<-ifelse(is.na(grid_pot[a,b]),0,grid_pot[a,b])+df$potency[i]; gc[a,b]<-gc[a,b]+1 }
grid_pot<-grid_pot/gc                                     # 每格平均潜能
arrows<-data.frame()
for(a in 2:(ngrid-1)) for(b in 2:(ngrid-1)){
  if(gc[a,b]<3) next
  # 中心差分梯度;分化方向 = 负梯度(潜能下降)
  dpx<-mean(c(grid_pot[a+1,b],grid_pot[a-1,b]),na.rm=TRUE); dpy<-mean(c(grid_pot[a,b+1],grid_pot[a,b-1]),na.rm=TRUE)
  gxv<- -( ( ifelse(is.na(grid_pot[a+1,b]),grid_pot[a,b],grid_pot[a+1,b]) -
             ifelse(is.na(grid_pot[a-1,b]),grid_pot[a,b],grid_pot[a-1,b])) )
  gyv<- -( ( ifelse(is.na(grid_pot[a,b+1]),grid_pot[a,b],grid_pot[a,b+1]) -
             ifelse(is.na(grid_pot[a,b-1]),grid_pot[a,b],grid_pot[a,b-1])) )
  mag<-sqrt(gxv^2+gyv^2); if(is.na(mag)||mag<1e-6) next
  sc<-0.6*(gx[2]-gx[1])/mag                               # 统一箭长
  arrows<-rbind(arrows, data.frame(x=gx[a],y=gy[b],xend=gx[a]+gxv*sc,yend=gy[b]+gyv*sc))
}
write.csv(df, file.path(DRES,"cell_potency.csv"), row.names=FALSE)
write.csv(arrows, file.path(DRES,"vector_field.csv"), row.names=FALSE)
cat(sprintf("[vector] %d grid arrows; potency range [%.2f,%.2f]\n", nrow(arrows), min(df$potency), max(df$potency)))

## ---- 4. 出图:潜能 UMAP + 分化方向矢量场 ---------------------------------
p <- ggplot(df, aes(UMAP1, UMAP2)) +
  geom_point(aes(color=potency), size=1.1, alpha=0.8) +
  scale_color_viridis_c(option="C", name="Potency\n(stemness)", direction=-1) +
  geom_segment(data=arrows, aes(x=x,y=y,xend=xend,yend=yend),
               arrow=arrow(length=unit(0.14,"cm"),type="closed"),
               linewidth=0.6, color="grey15", alpha=0.9) +
  labs(title="VECTOR differentiation direction",
       subtitle="arrows = inferred differentiation flow (high->low potency)") +
  theme_pub(base_size=11) + theme(axis.text=element_blank(), axis.ticks=element_blank())
save_fig(p, file.path(DAST,"vector_field"), width=5.6, height=4.6)
cat("[fig] assets/vector_field.{pdf,png}\n")
sink(file.path(DRES,"sessionInfo.txt")); print(sessionInfo()); sink()
