# =============================================================================
# 514 · 哑铃图 + 斜率图 (dumbbell & slopegraph) — 配对/前后变化的顶刊画法
# -----------------------------------------------------------------------------
# 两类替代"分组条形图"的配对变化图:
#   · 哑铃图 dumbbell:每个条目两点(如 before/after、control/treated)连线,
#     一眼看出每项的方向与幅度,按变化量排序;
#   · 斜率图 slopegraph:两时点间连线,凸显"谁升谁降、排名交叉"。
# 适合:处理前后的通路评分、两队列的效应、两时点的丰度等。
#
# Turnkey: Rscript 514_dumbbell_slope_plot.R   (合成配对数据→results/+assets/)
#          换数据: --input data.csv  (列: item,cond1,cond2;或长表 item,condition,value)
# 复用 _framework/theme_pub.R;无条形图。
# =============================================================================
suppressWarnings(suppressMessages({ library(ggplot2); library(ggrepel) }))

.this <- tryCatch({ a <- commandArgs(FALSE); m <- grep("^--file=", a)
  if(length(m)) dirname(normalizePath(sub("^--file=","",a[m[1]]))) else getwd() },
  error=function(e) getwd())
.fw <- NULL; .p <- .this
for(i in 1:6){ cand <- file.path(.p,"_framework","theme_pub.R"); if(file.exists(cand)){ .fw <- cand; break }; .p <- dirname(.p) }
if(!is.null(.fw)) source(.fw) else stop("需要 _framework/theme_pub.R")
args <- bio_args(list(input=NULL, c1="Control", c2="Treated"))

set.seed(42)
DIR <- .this
DDAT <- file.path(DIR,"example_data"); DRES <- file.path(DIR,"results"); DAST <- file.path(DIR,"assets")
for(d in c(DDAT,DRES,DAST)) dir.create(d, showWarnings=FALSE, recursive=TRUE)

## ---- 1. 合成:12 个条目两条件配对值 ---------------------------------------
fD <- file.path(DDAT,"paired.csv")
if(is.null(args$input)){
  if(!file.exists(fD)){
    items <- sprintf("Pathway %02d", 1:12)
    c1 <- runif(12, 1.5, 5); delta <- rnorm(12, 0.6, 1.0); c2 <- pmax(0.3, c1+delta)
    df <- data.frame(item=items, Control=round(c1,2), Treated=round(c2,2))
    write.csv(df, fD, row.names=FALSE)
    cat("[gen] synthetic: 12 items, two conditions (demo only)\n")
  }
  df <- read.csv(fD)
} else { df <- read.csv(args$input) }
names(df)[2:3] <- c("c1","c2")
df$change <- df$c2 - df$c1
df$dir <- ifelse(df$change>=0, "up", "down")
df <- df[order(df$c2),]; df$item <- factor(df$item, levels=df$item)
write.csv(df, file.path(DRES,"paired_change.csv"), row.names=FALSE)
cat(sprintf("[change] up=%d down=%d (median |Δ|=%.2f)\n", sum(df$dir=="up"), sum(df$dir=="down"), median(abs(df$change))))

col <- pal_pub(2,"npg"); names(col) <- c("up","down")
## ---- 2. 哑铃图(两点连线,按 cond2 排序,方向着色)------------------------
p1 <- ggplot(df) +
  geom_segment(aes(x=c1, xend=c2, y=item, yend=item, color=dir), linewidth=1.1, alpha=0.6) +
  geom_point(aes(c1, item), color="grey55", size=3) +
  geom_point(aes(c2, item, color=dir), size=3) +
  scale_color_manual(values=col, name="Change") +
  labs(x="Value", y=NULL, title=sprintf("Dumbbell: %s -> %s", args$c1, args$c2),
       subtitle="grey = baseline, colored = endpoint") +
  theme_pub(base_size=11)
save_fig(p1, file.path(DAST,"dumbbell"), width=5.8, height=4.2)

## ---- 3. 斜率图(两时点连线 + 端点标签)------------------------------------
long <- rbind(data.frame(item=df$item, x=args$c1, value=df$c1, dir=df$dir),
              data.frame(item=df$item, x=args$c2, value=df$c2, dir=df$dir))
long$x <- factor(long$x, levels=c(args$c1, args$c2))
p2 <- ggplot(long, aes(x, value, group=item, color=dir)) +
  geom_line(linewidth=0.9, alpha=0.7) + geom_point(size=2.4) +
  ggrepel::geom_text_repel(data=subset(long, x==args$c2), aes(label=item),
                           size=2.8, direction="y", nudge_x=0.15, segment.size=0.2) +
  scale_color_manual(values=col, name="Change") +
  labs(x=NULL, y="Value", title="Slopegraph") +
  theme_pub(base_size=11)
save_fig(p2, file.path(DAST,"slopegraph"), width=5.2, height=4.4)

cat("[fig] assets/: dumbbell, slopegraph (.pdf+.png)\n")
sink(file.path(DRES,"sessionInfo.txt")); print(sessionInfo()); sink()
