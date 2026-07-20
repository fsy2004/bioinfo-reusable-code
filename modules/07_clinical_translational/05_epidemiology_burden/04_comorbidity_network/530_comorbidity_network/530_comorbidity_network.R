# =============================================================================
# 530 · 共病网络 (disease-pair association → igraph network → Louvain community)
# -----------------------------------------------------------------------------
# 从"每患者-疾病"长表算疾病对共现关联(2x2 列联表 → phi / odds-ratio(Haldane 校正)/
# Jaccard / Fisher p),构建 igraph 共病网络,Louvain 社区检测与中心性,出:
#   ① 共病网络图(ggraph FR 布局,节点按社区着色,边宽=关联强度);
#   ② 疾病×疾病关联热图(log2(OR) 发散,中心=OR 1);③ 加权度 hub lollipop。
#
# 接地于真实工具代码: 21/99_external_sources/comorbidity_networks/sample_code.R
#   (graph_from_adjacency_matrix(t(mat),mode=...,diag=F,weighted=T) → delete_edges(weight<1)
#    → decompose(min.vertices=3) 的真实 igraph 写法) 及 CSB-IG_Comorbidity_Networks/Jaccard.R。
# 诚实边界(见 README):关联非因果;有向 OR 网做 Louvain 须先对称化、权重须非负;
#   零格用 Haldane 0.5 校正;OR>1 过滤会丢保护性关联(本模块同时给出 phi)。
#
# Turnkey: Rscript 530_comorbidity_network.R        (合成 患者-疾病 长表 → results/+assets/)
#          换数据: --input patients.csv  (两列: patient_id, disease)
# 复用 _framework/theme_pub.R;无条形图(network/heatmap/lollipop)。
# =============================================================================
suppressWarnings(suppressMessages({ library(ggplot2); library(dplyr); library(igraph); library(ggraph) }))

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

## ---- 1. 合成 患者-疾病 长表(3 个共病簇,簇内疾病更易共现)-----------------
cat("Step 1  合成 患者-疾病 长表 (含共病簇结构)\n")
clusters <- list(
  Cardiometabolic = c("Hypertension","Diabetes","Obesity","Dyslipidemia","CHD"),
  Respiratory     = c("COPD","Asthma","SleepApnea"),
  MentalMSK       = c("Depression","Anxiety","Osteoarthritis","BackPain","Osteoporosis"))
all_dz <- unlist(clusters, use.names=FALSE)
N <- 2500; rows <- list(); ri <- 1
for(pid in 1:N){
  k <- sample(seq_along(clusters), 1, prob=c(.45,.25,.30))         # 主簇
  for(ci in seq_along(clusters)){
    base <- if(ci==k) 0.5 else 0.08                                # 主簇内高共现
    got <- clusters[[ci]][runif(length(clusters[[ci]])) < base]
    if(length(got)) { rows[[ri]] <- data.frame(patient_id=pid, disease=got); ri <- ri+1 }
  }
}
long <- do.call(rbind, rows)
long <- long[!duplicated(long[,c("patient_id","disease")]),]
write.csv(long, file.path(DDAT,"patients.csv"), row.names=FALSE)

args <- commandArgs(TRUE); i <- match("--input", args)
if(!is.na(i) && i<length(args)) long <- read.csv(args[i+1], stringsAsFactors=FALSE)

## ---- 2. 患者×疾病 0/1 矩阵 + 疾病对 2x2 关联 (phi/OR-Haldane/Jaccard/Fisher p)
cat("Step 2  疾病对 2x2 关联 (phi / OR-Haldane / Jaccard / Fisher p)\n")
dz <- sort(unique(long$disease)); Np <- length(unique(long$patient_id))
M <- table(factor(long$patient_id), factor(long$disease, levels=dz)); M <- (M>0)*1
pair <- combn(dz, 2); res <- list()
for(j in 1:ncol(pair)){
  x <- M[,pair[1,j]]; y <- M[,pair[2,j]]
  a <- sum(x&y); b <- sum(x&!y); c <- sum(!x&y); d <- sum(!x&!y)
  OR  <- ((a+0.5)*(d+0.5))/((b+0.5)*(c+0.5))            # Haldane 0.5 校正(防零格)
  phi <- (a*d - b*c)/sqrt((a+b)*(c+d)*(a+c)*(b+d))
  jac <- a/(a+b+c)
  pfish <- fisher.test(matrix(c(a,b,c,d),2))$p.value
  res[[j]] <- data.frame(d1=pair[1,j], d2=pair[2,j], a=a, OR=OR, phi=phi, jaccard=jac, p=pfish)
}
assoc <- do.call(rbind, res)
assoc$padj <- p.adjust(assoc$p, "BH")                    # 多重检验校正(铁律)
write.csv(assoc, file.path(DRES,"pairwise_association.csv"), row.names=FALSE)

## ---- 3. 构建 igraph 共病网络 (正关联 OR>1 且 padj<0.05; 权重=log2(OR)) -------
## 写法接地 sample_code.R: 邻接 → 阈值 delete → 最大连通分量
cat("Step 3  igraph 网络 + Louvain 社区\n")
edges <- assoc %>% filter(OR>1, padj<0.05) %>% transmute(from=d1, to=d2, weight=log2(OR))
g <- graph_from_data_frame(edges, directed=FALSE, vertices=data.frame(name=dz))
g <- delete_vertices(g, V(g)[igraph::degree(g)==0])     # 去孤立点
# Louvain 需无向+非负权重(本图已满足)
cl <- cluster_louvain(g, weights=E(g)$weight)
V(g)$module <- as.factor(membership(cl))
V(g)$wdeg   <- strength(g, weights=E(g)$weight)         # 加权度
hub <- data.frame(disease=V(g)$name, module=V(g)$module, wdeg=V(g)$wdeg,
                  betw=betweenness(g, weights=1/E(g)$weight))
write.csv(hub[order(-hub$wdeg),], file.path(DRES,"node_metrics.csv"), row.names=FALSE)
cat(sprintf("   nodes=%d edges=%d modules=%d (modularity=%.2f)\n",
            vcount(g), ecount(g), length(cl), modularity(cl)))

set.seed(42)   # FR 布局随机 → 固定种子可复现
p1 <- ggraph(g, layout="fr") +
  geom_edge_link(aes(width=weight), alpha=0.35, color="grey55") +
  geom_node_point(aes(color=module, size=wdeg)) +
  geom_node_text(aes(label=name), repel=TRUE, size=2.6) +
  scale_edge_width(range=c(0.3,2.2), guide="none") +
  scale_color_pub("npg") + scale_size(range=c(2,7), guide="none") +
  labs(color="Module", title="Comorbidity network (Louvain modules)") +
  theme_void(base_size=11) + theme(plot.title=element_text(face="bold", hjust=0.5))
save_fig(p1, file.path(DAST,"comorbidity_network"), width=6.2, height=5.2)

## ---- 4. 疾病×疾病关联热图 (log2(OR) 发散, 中心=OR 1) ------------------------
cat("Step 4  关联热图 (log2 OR 发散)\n")
heat <- assoc %>% transmute(d1, d2, l2or=log2(OR))
heat2 <- rbind(heat, data.frame(d1=heat$d2, d2=heat$d1, l2or=heat$l2or))
ord <- hub$disease[order(hub$module, -hub$wdeg)]
ord <- c(ord, setdiff(dz, ord))
heat2$d1 <- factor(heat2$d1, levels=ord); heat2$d2 <- factor(heat2$d2, levels=ord)
p2 <- ggplot(heat2, aes(d1, d2, fill=l2or)) + geom_tile(color="white", linewidth=0.3) +
  scale_fill_diverge(midpoint=0, name="log2(OR)") +
  labs(x=NULL, y=NULL, title="Disease-pair association (log2 OR)") +
  theme_pub(base_size=9) + theme(axis.text.x=element_text(angle=45, hjust=1))
save_fig(p2, file.path(DAST,"association_heatmap"), width=5.6, height=5.0)

## ---- 5. 加权度 hub lollipop -----------------------------------------------
cat("Step 5  加权度 hub lollipop\n")
hb <- hub[order(hub$wdeg),]; hb$disease <- factor(hb$disease, levels=hb$disease)
p3 <- ggplot(hb, aes(wdeg, disease, color=module)) +
  geom_segment(aes(x=0, xend=wdeg, yend=disease), linewidth=0.6) +
  geom_point(size=2.6) + scale_color_pub("npg") +
  labs(x="Weighted degree (sum log2 OR)", y=NULL, color="Module",
       title="Comorbidity hubs") + theme_pub()
save_fig(p3, file.path(DAST,"hub_lollipop"), width=5.0, height=3.8)

cat("Done 530 · figures → assets/ , tables → results/\n")
