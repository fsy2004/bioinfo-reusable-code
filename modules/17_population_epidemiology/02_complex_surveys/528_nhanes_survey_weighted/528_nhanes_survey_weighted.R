# =============================================================================
# 528 · NHANES 复杂抽样加权分析 (survey design / svymean / svyby / svyglm)
# -----------------------------------------------------------------------------
# 用 survey 包按 NHANES 复杂抽样设计(权重 + 分层 + PSU)做加权描述与加权回归:
#   ① 加权 vs 未加权分组均值(dumbbell,显示忽略权重的偏差);
#   ② svyglm 设计校正回归 → 系数森林图;③ 加权患病率(0/1 结局 svymean)→ lollipop。
#
# 接地于真实工具代码: 21/99_external_sources/nhanes/vignettes/UsingSurveyWeights.{rmd,R}
#   (svydesign→subset→svymean/svyby/svyquantile→svyglm→系数森林的全套真实调用)
#   及 nhanes/R/nhanes.R(SEQN 为合并键、translated 因子化的数据契约)。
# 诚实边界(见 README):横断面关联非因果;先建 design 再 subset(不可先筛行);
#   权重须配数据来源(MEC 检查项用 WTMEC2YR,访谈项用 WTINT2YR)。
#
# Turnkey: Rscript 528_nhanes_survey_weighted.R     (合成 NHANES 形状数据 → results/+assets/)
#          换数据: --input nhanes.csv  (列见 README;须含 SEQN/SDMVPSU/SDMVSTRA/WTMEC2YR/...)
# 复用 _framework/theme_pub.R;无条形图(dumbbell/forest/lollipop)。
# =============================================================================
suppressWarnings(suppressMessages({ library(ggplot2); library(dplyr); library(survey) }))

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
options(survey.lonely.psu = "adjust")   # 单PSU分层:避免 svyglm/svymean 报错(真实坑)

## ---- 1. 合成 NHANES 形状数据(SEQN/PSU/strata/MEC权重/人口学/结局)----------
cat("Step 1  合成 NHANES 复杂抽样数据\n")
n <- 4000
eth_lv <- c("Mexican American","Other Hispanic","Non-Hispanic White","Non-Hispanic Black","Other")
strata <- sample(1:15, n, replace=TRUE)
psu    <- sample(1:2, n, replace=TRUE)                 # PSU 仅在层内唯一 → nest=TRUE
age    <- round(runif(n, 18, 80))
sex    <- sample(c("Male","Female"), n, replace=TRUE)
# 少数族裔被过抽(设计如此)→ 抽样权重更大;制造加权≠未加权
eth    <- sample(eth_lv, n, replace=TRUE, prob=c(.18,.12,.30,.25,.15))
base_w <- c(`Mexican American`=2.4,`Other Hispanic`=2.6,`Non-Hispanic White`=1.0,
            `Non-Hispanic Black`=1.8,`Other`=2.0)[eth]
WTMEC2YR <- as.numeric(base_w) * 8000 * exp(rnorm(n,0,0.25))
# 结局:舒张压随年龄/性别/族裔变化(族裔效应与过抽相关→加权修偏)
eth_eff <- c(`Mexican American`=2,`Other Hispanic`=1,`Non-Hispanic White`=0,
             `Non-Hispanic Black`=6,`Other`=1)[eth]
BPXDI1 <- 70 + 0.10*(age-50) + ifelse(sex=="Male",3,0) + as.numeric(eth_eff) + rnorm(n,0,9)
HTN    <- as.integer(BPXDI1 > 80 | runif(n) < plogis(-2 + 0.04*(age-50) + 0.3*(eth=="Non-Hispanic Black")))
dat <- data.frame(SEQN=1:n, SDMVSTRA=strata, SDMVPSU=psu, WTMEC2YR=WTMEC2YR,
                  RIDAGEYR=age, RIAGENDR=sex, RIDRETH1=factor(eth, levels=eth_lv),
                  BPXDI1=round(BPXDI1,1), HTN=HTN)
write.csv(dat, file.path(DDAT,"nhanes.csv"), row.names=FALSE)

args <- commandArgs(TRUE); i <- match("--input", args)
if(!is.na(i) && i<length(args)){ dat <- read.csv(args[i+1], stringsAsFactors=FALSE)
  dat$RIDRETH1 <- factor(dat$RIDRETH1) }

## ---- 2. 复杂抽样设计对象(★先建 design,再 subset;含 weights+strata+PSU)----
cat("Step 2  svydesign (weights + strata + PSU, nest=TRUE)\n")
des  <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTMEC2YR, nest=TRUE, data=dat)
desa <- subset(des, RIDAGEYR >= 20)                   # 在 design 上 subset,保持方差结构

## ---- 3. 加权 vs 未加权分组均值 → dumbbell ----------------------------------
cat("Step 3  加权 vs 未加权均值 (dumbbell)\n")
wm <- svyby(~BPXDI1, ~RIDRETH1, desa, svymean, na.rm=TRUE)
wtd  <- data.frame(eth=as.character(wm$RIDRETH1), mean=as.numeric(wm$BPXDI1), type="Weighted")
uw   <- dat %>% filter(RIDAGEYR>=20) %>% group_by(eth=as.character(RIDRETH1)) %>%
        summarise(mean=mean(BPXDI1), .groups="drop") %>% mutate(type="Unweighted")
db <- bind_rows(wtd, uw); db$eth <- factor(db$eth, levels=eth_lv)
p1 <- ggplot(db, aes(mean, eth)) +
  geom_line(aes(group=eth), color="grey70", linewidth=1) +
  geom_point(aes(color=type), size=3) + scale_color_pub("npg") +
  labs(x="Mean diastolic BP (mmHg)", y=NULL, color=NULL,
       title="Survey-weighted vs unweighted means") + theme_pub()
save_fig(p1, file.path(DAST,"weighted_vs_unweighted_dumbbell"), width=5.4, height=3.2)

## ---- 4. svyglm 设计校正回归 → 系数森林图 -----------------------------------
cat("Step 4  svyglm 设计校正回归 (系数森林)\n")
fit <- svyglm(BPXDI1 ~ RIDAGEYR + RIAGENDR + RIDRETH1, design=desa)
co  <- summary(fit)$coefficients
cf  <- data.frame(term=rownames(co), est=co[,1], se=co[,2])
cf  <- cf[cf$term!="(Intercept)",]
cf$lo <- cf$est-1.96*cf$se; cf$hi <- cf$est+1.96*cf$se
cf$term <- gsub("RIDRETH1","Eth: ",cf$term); cf$term <- gsub("RIAGENDR","Sex: ",cf$term)
cf$term <- gsub("RIDAGEYR","Age (per year)",cf$term)
cf <- cf[order(cf$est),]; cf$term <- factor(cf$term, levels=cf$term)
write.csv(cf, file.path(DRES,"svyglm_coefficients.csv"), row.names=FALSE)
p2 <- ggplot(cf, aes(est, term)) +
  geom_vline(xintercept=0, linetype=2, color="grey50") +
  geom_errorbarh(aes(xmin=lo, xmax=hi), height=0.18, color="#0072B5", linewidth=0.7) +
  geom_point(size=2.6, color="#BC3C29") +
  labs(x="Design-adjusted coefficient (95% CI)", y=NULL,
       title="svyglm: diastolic BP associations") + theme_pub()
save_fig(p2, file.path(DAST,"svyglm_forest"), width=5.2, height=3.4)

## ---- 5. 加权患病率 (0/1 结局 svymean) by 族裔 → lollipop --------------------
cat("Step 5  加权 HTN 患病率 (lollipop + SE)\n")
pv <- svyby(~HTN, ~RIDRETH1, desa, svymean, na.rm=TRUE)
prev <- data.frame(eth=factor(as.character(pv$RIDRETH1), levels=eth_lv),
                   prev=as.numeric(pv$HTN), se=as.numeric(pv$se))
prev <- prev[order(prev$prev),]; prev$eth <- factor(as.character(prev$eth), levels=as.character(prev$eth))
write.csv(prev, file.path(DRES,"weighted_prevalence.csv"), row.names=FALSE)
p3 <- ggplot(prev, aes(prev, eth)) +
  geom_segment(aes(x=0, xend=prev, yend=eth), color="grey70", linewidth=0.7) +
  geom_errorbarh(aes(xmin=prev-1.96*se, xmax=prev+1.96*se), height=0.15, color="#0072B5") +
  geom_point(size=3, color="#BC3C29") +
  scale_x_continuous(labels=scales::percent) +
  labs(x="Weighted prevalence (95% CI)", y=NULL,
       title="Survey-weighted prevalence by subgroup") + theme_pub()
save_fig(p3, file.path(DAST,"weighted_prevalence_lollipop"), width=5.2, height=3.0)

cat("Done 528 · figures → assets/ , tables → results/\n")
