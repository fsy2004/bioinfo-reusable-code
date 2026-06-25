# =============================================================================
# 513 · 山脊图 (ridgeline / joyplot) — 多组分布的紧凑堆叠展示
# -----------------------------------------------------------------------------
# 山脊图把多组(常 5-20 组)的密度分布上下错开堆叠,适合展示"分布随分期/时间/
# 细胞类型的渐变",比并排小提琴更省空间、比条形图信息量大得多。支持按数值梯度
# 着色(viridis)凸显分布位移。
#
# Turnkey: Rscript 513_ridgeline_plot.R   (合成多组分布→results/+assets/)
#          换数据: --input data.csv  (列: group,value)
# 复用 _framework/theme_pub.R;依赖 ggridges;无条形图。
# =============================================================================
suppressWarnings(suppressMessages({ library(ggplot2); library(ggridges) }))

.this <- tryCatch({ a <- commandArgs(FALSE); m <- grep("^--file=", a)
  if(length(m)) dirname(normalizePath(sub("^--file=","",a[m[1]]))) else getwd() },
  error=function(e) getwd())
.fw <- NULL; .p <- .this
for(i in 1:6){ cand <- file.path(.p,"_framework","theme_pub.R"); if(file.exists(cand)){ .fw <- cand; break }; .p <- dirname(.p) }
if(!is.null(.fw)) source(.fw) else stop("需要 _framework/theme_pub.R")
args <- bio_args(list(input=NULL))

set.seed(42)
DIR <- .this
DDAT <- file.path(DIR,"example_data"); DRES <- file.path(DIR,"results"); DAST <- file.path(DIR,"assets")
for(d in c(DDAT,DRES,DAST)) dir.create(d, showWarnings=FALSE, recursive=TRUE)

## ---- 1. 合成:8 个时间点/细胞类型,分布渐变位移 ---------------------------
fD <- file.path(DDAT,"groups_value.csv")
if(is.null(args$input)){
  if(!file.exists(fD)){
    grps <- sprintf("Day %d", c(0,2,4,7,10,14,21,28))
    df <- do.call(rbind, lapply(seq_along(grps), function(i)
      data.frame(group=grps[i], value=rnorm(200, 2 + 0.45*i + ifelse(i>5,0.5,0), 0.8+0.03*i))))
    write.csv(df, fD, row.names=FALSE)
    cat("[gen] synthetic: 8 timepoints, drifting distributions (demo only)\n")
  }
  df <- read.csv(fD)
} else df <- read.csv(args$input)
df$group <- factor(df$group, levels=rev(unique(df$group)))   # 第一个组在顶部

## ---- 2. 出图(梯度着色山脊 + 中位线)--------------------------------------
p <- ggplot(df, aes(x=value, y=group, fill=after_stat(x))) +
  ggridges::geom_density_ridges_gradient(scale=2.3, rel_min_height=0.01,
                                         quantile_lines=TRUE, quantiles=2,
                                         color="white", linewidth=0.3) +
  scale_fill_viridis_c(option="D", name="Value") +
  labs(x="Value", y=NULL, title="Ridgeline plot (distribution over time)") +
  theme_pub(base_size=11) +
  theme(axis.text.y=element_text(vjust=0))
save_fig(p, file.path(DAST,"ridgeline"), width=6, height=4.6)

write.csv(aggregate(value~group, df, function(v) round(c(median=median(v), IQR=IQR(v)),3)),
          file.path(DRES,"group_summary.csv"), row.names=FALSE)
cat("[fig] assets/ridgeline.{pdf,png}\n")
sink(file.path(DRES,"sessionInfo.txt")); print(sessionInfo()); sink()
