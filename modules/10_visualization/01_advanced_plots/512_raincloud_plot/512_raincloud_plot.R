# =============================================================================
# 512 · 云雨图 (raincloud plot) — 顶刊级分布可视化,替代条形图
# -----------------------------------------------------------------------------
# 云雨图 = 半小提琴(密度"云")+ 箱线(四分位)+ 抖动点("雨"),一张图同时给出
# 分布形状、中位/四分位、原始样本点,信息密度远高于条形图。适合组间评分/表达/
# 指标比较(如不同分期/细胞类型/处理)。可选叠加显著性标注。
#
# Turnkey: Rscript 512_raincloud_plot.R   (合成多组评分→results/+assets/)
#          换数据: --input data.csv  (列: group,value;可选 --value col --group col)
# 复用 _framework/theme_pub.R;依赖 ggdist(半小提琴/点);无条形图。
# =============================================================================
suppressWarnings(suppressMessages({ library(ggplot2); library(ggdist) }))

.this <- tryCatch({ a <- commandArgs(FALSE); m <- grep("^--file=", a)
  if(length(m)) dirname(normalizePath(sub("^--file=","",a[m[1]]))) else getwd() },
  error=function(e) getwd())
.fw <- NULL; .p <- .this
for(i in 1:6){ cand <- file.path(.p,"_framework","theme_pub.R"); if(file.exists(cand)){ .fw <- cand; break }; .p <- dirname(.p) }
if(!is.null(.fw)) source(.fw) else stop("需要 _framework/theme_pub.R")
args <- bio_args(list(input=NULL, value="value", group="group"))

set.seed(42)
DIR <- .this
DDAT <- file.path(DIR,"example_data"); DRES <- file.path(DIR,"results"); DAST <- file.path(DIR,"assets")
for(d in c(DDAT,DRES,DAST)) dir.create(d, showWarnings=FALSE, recursive=TRUE)

## ---- 1. 读 / 生成多组分布数据(4 组,均值与离散度不同)---------------------
fD <- file.path(DDAT,"groups_value.csv")
if(is.null(args$input)){
  if(!file.exists(fD)){
    grps <- c("Healthy","Stage I","Stage II","Stage III")
    mu <- c(2.0,3.2,4.1,5.3); sdv <- c(0.7,0.9,1.1,1.0); n <- c(80,70,65,60)
    df <- do.call(rbind, lapply(seq_along(grps), function(i)
      data.frame(group=grps[i], value=rnorm(n[i], mu[i], sdv[i]))))
    write.csv(df, fD, row.names=FALSE)
    cat("[gen] synthetic: 4 groups, graded score (demo only)\n")
  }
  df <- read.csv(fD)
} else { df <- read.csv(args$input); names(df)[match(c(args$value,args$group), names(df))] <- c("value","group") }
df$group <- factor(df$group, levels=unique(df$group))

## ---- 2. 组间统计(Kruskal-Wallis + 两两 Wilcoxon, BH 校正)------------------
kw <- kruskal.test(value~group, df)
pw <- pairwise.wilcox.test(df$value, df$group, p.adjust.method="BH")
write.csv(data.frame(group=levels(df$group),
                     median=tapply(df$value, df$group, median),
                     IQR=tapply(df$value, df$group, IQR)),
          file.path(DRES,"group_summary.csv"), row.names=FALSE)
capture.output(print(pw), file=file.path(DRES,"pairwise_wilcoxon.txt"))
cat(sprintf("[stat] Kruskal-Wallis p=%.2g\n", kw$p.value))

## ---- 3. 云雨图(half-eye 云 + box + 抖动雨;水平方向)-----------------------
col <- pal_pub(nlevels(df$group), "okabe_ito")
p <- ggplot(df, aes(group, value, fill=group, color=group)) +
  ggdist::stat_halfeye(adjust=0.7, width=0.6, .width=0, justification=-0.25,
                       point_colour=NA, alpha=0.7) +
  geom_boxplot(width=0.14, outlier.shape=NA, alpha=0.5, linewidth=0.5) +
  geom_jitter(width=0.06, size=0.7, alpha=0.4) +
  scale_fill_manual(values=col, guide="none") +
  scale_color_manual(values=col, guide="none") +
  coord_flip() +
  labs(x=NULL, y="Score", title="Raincloud plot",
       subtitle=sprintf("Kruskal-Wallis p = %.2g", kw$p.value)) +
  theme_pub(base_size=11)
save_fig(p, file.path(DAST,"raincloud"), width=6, height=4.2)
cat("[fig] assets/raincloud.{pdf,png}\n")
sink(file.path(DRES,"sessionInfo.txt")); print(sessionInfo()); sink()
