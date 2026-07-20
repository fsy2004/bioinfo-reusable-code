# =============================================================================
# 529 · CHARLS 纵向队列 (wave 描述 / 趋势 / 等百分位等值化 / 混合模型 / 生存)
# -----------------------------------------------------------------------------
# CHARLS 式重复测量随访队列的标准产出:
#   ① 波次分层描述(Table 1);② 重复结局的纵向趋势(均值±SE,原始 vs 等值化);
#   ③ 跨波等百分位等值化(把后续波分数 crosswalk 到参照波刻度);
#   ④ lme4 线性混合模型轨迹(随机截距);⑤ 基线多病 → 事件(CVD)的 Cox + KM。
#
# 接地于真实工具代码: 21/99_external_sources/charls_memory_equating/scripts/*.R
#   (按波合并 → 清洗重复测量 → 校准样本 → 加权等值化 → 趋势图 的真实流程)。
# 诚实边界(见 README):equate 包未装 → 等百分位用加权 ECDF 反演自实现(标准等价做法);
#   grip-IPD 仓库仅有 00_setup_paths.R,其 Cox/multistate/joint 脚本不在盘上 → 不臆造,
#   本模块用 survival(已装)做基础 Cox+KM;关联非因果。
#
# Turnkey: Rscript 529_charls_longitudinal_cohort.R (合成 长面板 → results/+assets/)
#          换数据: --input panel.csv  (列见 README;长格式 一人一波一行)
# 复用 _framework/theme_pub.R;无条形图(line/violin/curve/forest/KM)。
# =============================================================================
suppressWarnings(suppressMessages({ library(ggplot2); library(dplyr); library(tidyr)
  library(survival); library(lme4) }))

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

## ---- 1. 合成 CHARLS 式长面板(ID×wave;认知分随龄/波下降 + 基线多病→CVD)----
cat("Step 1  合成 CHARLS 长面板\n")
waves <- 1:5; years <- c(2011,2013,2015,2018,2020); n_id <- 1200
base <- data.frame(ID=1:n_id, female=rbinom(n_id,1,.52),
                   edu=sample(0:3, n_id, replace=TRUE, prob=c(.4,.3,.2,.1)),
                   rural=rbinom(n_id,1,.55), age0=round(runif(n_id,45,75)),
                   u=rnorm(n_id,0,1))                         # 个体随机截距
mm0 <- rpois(n_id, lambda=0.5 + 0.03*(base$age0-45))          # 基线多病计数
panel <- list()
for(w in waves){
  age <- base$age0 + (years[w]-2011)
  # 认知(即时回忆 0-10):随龄/波下降,女性/教育↑;含个体随机效应 u
  mu <- 7.2 - 0.06*(age-60) + 0.4*base$female + 0.5*base$edu - 0.15*(w-1) + base$u
  score <- pmin(10, pmax(0, round(rnorm(n_id, mu, 1.4))))
  # 后续波"难度漂移":2018/2020 测验更难(整体偏低)→ 需等值化校正
  if(years[w] %in% c(2018,2020)) score <- pmin(10, pmax(0, round(score - 1.0)))
  wt <- exp(rnorm(n_id, 0, 0.3)) * 8000                       # 个体抽样权重
  panel[[w]] <- data.frame(ID=base$ID, wave=w, year=years[w], age=age,
    female=base$female, edu=base$edu, rural=base$rural, mm0=mm0, score=score, wtresp=wt)
}
panel <- do.call(rbind, panel)
# 基线多病 → incident CVD 生存(基线进入,随访至 2020)
haz <- 0.04*exp(0.35*mm0 + 0.02*(base$age0-60))
ttime <- pmin(rexp(n_id, haz), 9); event <- as.integer(ttime < 9)
surv <- data.frame(ID=base$ID, time=round(ttime,2), event=event,
                   mm_group=ifelse(mm0>=2,"Multimorbid (>=2)","0-1 condition"), age0=base$age0)
write.csv(panel, file.path(DDAT,"panel.csv"), row.names=FALSE)
write.csv(surv,  file.path(DDAT,"survival.csv"), row.names=FALSE)

args <- commandArgs(TRUE); i <- match("--input", args)
if(!is.na(i) && i<length(args)) panel <- read.csv(args[i+1], stringsAsFactors=FALSE)

## ---- 2. 波次分层描述 (Table 1) ---------------------------------------------
cat("Step 2  波次分层 Table 1\n")
tab1 <- panel %>% group_by(year) %>% summarise(n=n(), age=sprintf("%.1f (%.1f)",mean(age),sd(age)),
  female_pct=sprintf("%.1f%%",100*mean(female)), score=sprintf("%.2f (%.2f)",mean(score),sd(score)),
  .groups="drop")
write.csv(tab1, file.path(DRES,"table1_by_wave.csv"), row.names=FALSE)

## ---- 3. 等百分位等值化 (加权 ECDF 反演; 替代未装的 equate 包) ----------------
cat("Step 3  等百分位等值化 (加权ECDF反演)\n")
# 加权百分位秩 (mid-point) 与 加权分位反函数
wpctrank <- function(s, w, at){
  o <- order(s); s<-s[o]; w<-w[o]; cw <- cumsum(w)/sum(w); below <- c(0,cw[-length(cw)])
  mid <- (below + cw)/2
  approx(s, mid, xout=at, rule=2, ties=function(x) mean(x))$y
}
wquantile <- function(s, w, p){
  o <- order(s); s<-s[o]; w<-w[o]; cw <- cumsum(w)/sum(w); below <- c(0,cw[-length(cw)])
  mid <- (below + cw)/2
  approx(mid, s, xout=p, rule=2, ties=function(x) mean(x))$y
}
ref_w <- 3                                                    # 参照波 = 2015
ref <- panel %>% filter(wave==ref_w)
xwalk <- lapply(c(4,5), function(w){                          # 把 2018/2020 crosswalk 到 2015
  lat <- panel %>% filter(wave==w)
  sc  <- 0:10
  pr  <- wpctrank(lat$score, lat$wtresp, sc)                  # 后续波分数的百分位
  eq  <- pmin(10, pmax(0, wquantile(ref$score, ref$wtresp, pr)))  # 反演到参照波刻度 + 顶/底封顶
  data.frame(wave=w, raw=sc, equated=eq)
})
xwalk <- do.call(rbind, xwalk)
write.csv(xwalk, file.path(DRES,"equipercentile_crosswalk.csv"), row.names=FALSE)
# 应用 crosswalk 得到 score_eqt
panel <- panel %>% left_join(xwalk, by=c("wave","score"="raw")) %>%
  mutate(score_eqt = ifelse(is.na(equated), score, equated))

## ---- 4. 纵向趋势 (原始 vs 等值化, 均值±SE) ---------------------------------
cat("Step 4  纵向趋势 (原始 vs 等值化)\n")
trend <- bind_rows(
  panel %>% group_by(year) %>% summarise(m=mean(score), se=sd(score)/sqrt(n()), type="Original",.groups="drop"),
  panel %>% group_by(year) %>% summarise(m=mean(score_eqt), se=sd(score_eqt)/sqrt(n()), type="Equated",.groups="drop"))
p1 <- ggplot(trend, aes(year, m, color=type, fill=type)) +
  geom_ribbon(aes(ymin=m-1.96*se, ymax=m+1.96*se), alpha=0.15, color=NA) +
  geom_line(linewidth=0.8) + geom_point(size=1.8) + scale_color_pub("npg") + scale_fill_pub("npg") +
  labs(x="Year", y="Immediate recall (0-10)", color=NULL, fill=NULL,
       title="Cognitive trajectory: original vs equated") + theme_pub()
save_fig(p1, file.path(DAST,"trend_original_vs_equated"), width=5.2, height=3.6)

## ---- 5. 等值化 concordance 曲线 -------------------------------------------
xwalk$wave_lab <- c(`4`="2018", `5`="2020")[as.character(xwalk$wave)]
p2 <- ggplot(xwalk, aes(raw, equated, color=wave_lab)) +
  geom_abline(slope=1, intercept=0, linetype=2, color="grey60") +
  geom_line(linewidth=0.8) + geom_point(size=1.6) +
  scale_color_pub("npg") +
  labs(x="Raw score (later wave)", y="Equated score (2015 scale)", color="Wave",
       title="Equipercentile crosswalk") + theme_pub()
save_fig(p2, file.path(DAST,"equipercentile_concordance"), width=4.6, height=3.6)

## ---- 6. 波次分数分布 (violin + box + jitter, raincloud 风) ------------------
cat("Step 6  波次分布 (violin/raincloud)\n")
p3 <- ggplot(panel, aes(factor(year), score, fill=factor(year))) +
  geom_violin(alpha=0.5, color=NA, width=0.9) +
  geom_boxplot(width=0.16, outlier.shape=NA, alpha=0.9) +
  geom_jitter(width=0.08, size=0.25, alpha=0.18) +
  scale_fill_pub("npg", guide="none") +
  labs(x="Year", y="Immediate recall (0-10)", title="Score distribution by wave") + theme_pub()
save_fig(p3, file.path(DAST,"score_distribution_violin"), width=5.2, height=3.4)

## ---- 7. lme4 线性混合模型轨迹 (随机截距) ----------------------------------
cat("Step 7  lme4 混合模型\n")
panel$year_c <- panel$year-2011; panel$age_c <- panel$age-60
m <- lmer(score ~ year_c + age_c + female + edu + (1|ID), data=panel, REML=TRUE)
fe <- summary(m)$coefficients
fm <- data.frame(term=rownames(fe), est=fe[,1], se=fe[,2])
fm <- fm[fm$term!="(Intercept)",]; fm$lo<-fm$est-1.96*fm$se; fm$hi<-fm$est+1.96*fm$se
fm <- fm[order(fm$est),]; fm$term <- factor(fm$term, levels=fm$term)
write.csv(fm, file.path(DRES,"lmer_fixed_effects.csv"), row.names=FALSE)
p4 <- ggplot(fm, aes(est, term)) + geom_vline(xintercept=0, linetype=2, color="grey50") +
  geom_errorbarh(aes(xmin=lo, xmax=hi), height=0.16, color="#0072B5") +
  geom_point(size=2.6, color="#BC3C29") +
  labs(x="Fixed-effect estimate (95% CI)", y=NULL, title="LMM: cognitive trajectory") + theme_pub()
save_fig(p4, file.path(DAST,"lmer_forest"), width=5.0, height=3.0)

## ---- 8. 基线多病 → incident CVD: Cox + KM ----------------------------------
cat("Step 8  Cox + KM (基线多病 → CVD)\n")
cox <- coxph(Surv(time,event) ~ I(mm0) + I(age0-60), data=surv)
write.csv(broom_cox <- data.frame(term=names(coef(cox)), HR=exp(coef(cox)),
          lo=exp(confint(cox)[,1]), hi=exp(confint(cox)[,2])), file.path(DRES,"cox_hr.csv"), row.names=FALSE)
fit <- survfit(Surv(time,event) ~ mm_group, data=surv)
km <- data.frame(time=fit$time, surv=fit$surv,
                 group=rep(names(fit$strata), fit$strata))
km$group <- sub("mm_group=","",km$group)
p5 <- ggplot(km, aes(time, surv, color=group)) + geom_step(linewidth=0.9) +
  scale_color_pub("npg") + ylim(0,1) +
  labs(x="Follow-up (years)", y="CVD-free probability", color="Baseline",
       title=sprintf("Incident CVD by baseline multimorbidity (HR/condition=%.2f)", exp(coef(cox)[1]))) +
  theme_pub()
save_fig(p5, file.path(DAST,"km_incident_cvd"), width=5.2, height=3.6)

cat("Done 529 · figures → assets/ , tables → results/\n")
