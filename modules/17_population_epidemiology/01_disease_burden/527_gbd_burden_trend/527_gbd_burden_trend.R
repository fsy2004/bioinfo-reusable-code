# =============================================================================
# 527 · GBD 疾病负担趋势 (age-standardized rate / EAPC / Das Gupta decomposition / SDI)
# -----------------------------------------------------------------------------
# 把 GBD 长表(measure/location/sex/age/metric/year/val/upper/lower)做成一套标准
# 负担分析:① 年龄标化率(ASR)时间趋势 + 95%UI 带;② EAPC(对数线性 lm(log(ASR)~year));
# ③ 年龄-性别结构(背靠背 lollipop,非金字塔条形);④ Das Gupta 三因子分解
# (老龄化/人口增长/流行病学变化);⑤ ASR–SDI 跨地区关联(Spearman + LOESS)。
#
# 接地于真实工具代码(modules/21_disease_burden_gbd/99_external_sources):
#   R-script-for-GBD/decomposition.R   → Das Gupta 公式逐字照搬
#   R-script-for-GBD/{Comparison of different types of trends, Bilateral diagram, Age group line chart}.R
#   GBD2021/.../SDI/SDI_incidence.R    → ASR~SDI Spearman+LOESS
# 诚实边界(见 README): Joinpoint=外部 NCI 软件不在 R 内复现 → 用对数线性 EAPC 替代;
#   BAPC 投影需 INLA(未装)→ 本模块不含投影;ASR 此处由年龄别率直接标化得到。
#
# Turnkey: Rscript 527_gbd_burden_trend.R           (合成 GBD 长表 → results/ + assets/)
#          换数据: --burden burden.csv --pop pop.csv --sdi sdi.csv
# 复用 _framework/theme_pub.R;无条形图(line/ribbon/lollipop/dumbbell/scatter)。
# =============================================================================
suppressWarnings(suppressMessages({ library(ggplot2); library(dplyr) }))

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

AGES <- c("<5","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49",
          "50-54","55-59","60-64","65-69","70-74","75-79","80-84","85-89","90-94","95+")
# 合成"世界标准人口"权重(随年龄递减,和=1;demo 用,真实分析换 GBD world standard)
WSTD <- { w <- exp(-(0:19)/6); w/sum(w) }
YEARS <- 1990:2021
LOCS  <- c("China", paste0("Region", sprintf("%02d",1:11)))

## ---- 1. 合成 GBD 长表(年龄别率 + 人口 → ASR 直接标化 + Number)-------------
cat("Step 1  合成 GBD 负担长表 (年龄别率/人口/标化)\n")
# 每地区一个基线 SDI(0.4~0.85),随年份缓升;率随 SDI 升而降(发病),随年龄升而升
sdi_base <- setNames(seq(0.42, 0.84, length.out=length(LOCS)), LOCS)
loc_noise <- setNames(exp(rnorm(length(LOCS), 0, 0.13)), LOCS)  # 地区级随机效应(打破完美单调)
age_grad <- seq(0.4, 3.2, length.out=20)                       # 年龄别率梯度
burden <- list(); popn <- list(); k <- 1; kp <- 1
for(loc in LOCS) for(sex in c("Male","Female")){
  sex_f <- ifelse(sex=="Male", 1.15, 0.9)
  for(yr in YEARS){
    sdi <- min(0.95, sdi_base[[loc]] + (yr-1990)*0.004)
    base_rate <- 60 * (1.1 - sdi) * sex_f * loc_noise[[loc]]   # /1e5 量级基线(含地区噪声)
    r_age <- base_rate * age_grad * exp(rnorm(20,0,0.05))      # 年龄别率 (/1e5)
    # 人口:总量随年份增长 + 老龄化(高龄占比随年份升)
    tot <- 1e6 * (1 + (yr-1990)*0.012) * ifelse(loc=="China", 12, 1)
    shape <- exp(-(0:19)/8) * (1 + (yr-1990)*0.012)^(0:19/4)   # 高龄随时间变厚
    a_prop <- shape/sum(shape); pop_age <- tot*a_prop
    asr <- sum(WSTD * r_age)                                   # 直接标化 ASR
    num_all <- sum(r_age/1e5 * pop_age)                        # All ages Number
    burden[[k]] <- data.frame(measure_name="Incidence", location_name=loc, sex_name=sex,
      age_name=c(AGES,"Age-standardized","All ages"),
      metric_name=c(rep("Rate",20),"Rate","Number"),
      year=yr, val=c(r_age, asr, num_all),
      sdi=sdi, stringsAsFactors=FALSE); k <- k+1
    popn[[kp]] <- data.frame(location_name=loc, sex_name=sex, year=yr,
      age_name=AGES, val=pop_age, stringsAsFactors=FALSE); kp <- kp+1
  }
}
burden <- do.call(rbind, burden); popn <- do.call(rbind, popn)
burden$upper <- burden$val*(1+0.07); burden$lower <- burden$val*(1-0.07)
sdi_tab <- unique(burden[,c("location_name","year","sdi")])
write.csv(burden, file.path(DDAT,"burden.csv"), row.names=FALSE)
write.csv(popn,   file.path(DDAT,"pop.csv"),    row.names=FALSE)
write.csv(sdi_tab,file.path(DDAT,"sdi.csv"),    row.names=FALSE)

# 支持换数据
args <- commandArgs(TRUE)
getarg <- function(k, d){ i <- match(k, args); if(!is.na(i) && i<length(args)) args[i+1] else d }
burden <- read.csv(getarg("--burden", file.path(DDAT,"burden.csv")), stringsAsFactors=FALSE)
popn   <- read.csv(getarg("--pop",    file.path(DDAT,"pop.csv")),    stringsAsFactors=FALSE)

## ---- 2. ASR 趋势 + 95%UI 带 (China, by sex) --------------------------------
cat("Step 2  ASR 时间趋势 + UI 带\n")
asr <- burden %>% filter(location_name=="China", age_name=="Age-standardized", metric_name=="Rate")
p1 <- ggplot(asr, aes(year, val, color=sex_name, fill=sex_name)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.15, color=NA) +
  geom_line(linewidth=0.8) + geom_point(size=1.1) +
  scale_color_pub("npg") + scale_fill_pub("npg") +
  labs(x="Year", y="ASR (per 100,000)", color="Sex", fill="Sex",
       title="Age-standardized incidence, China 1990-2021") + theme_pub()
save_fig(p1, file.path(DAST,"asr_trend"), width=5.4, height=3.8)

## ---- 3. EAPC = 100*(exp(beta)-1), lm(log(ASR) ~ year) (对数线性,非Joinpoint)-
cat("Step 3  EAPC (对数线性) + 95%CI\n")
eapc <- asr %>% group_by(sex_name, location_name=="China") %>% group_modify(~{
  m <- lm(log(val) ~ year, data=.x); b <- coef(m)["year"]; ci <- confint(m)["year",]
  data.frame(eapc=100*(exp(b)-1), lo=100*(exp(ci[1])-1), hi=100*(exp(ci[2])-1))
}) %>% ungroup() %>% select(sex_name, eapc, lo, hi)
# 跨地区 EAPC 也算(供 SDI 关联展示)
eapc_loc <- burden %>% filter(age_name=="Age-standardized", metric_name=="Rate", sex_name=="Female") %>%
  group_by(location_name) %>% group_modify(~{ m <- lm(log(val) ~ year, data=.x)
    data.frame(eapc=100*(exp(coef(m)["year"])-1)) }) %>% ungroup()
write.csv(eapc, file.path(DRES,"eapc_china.csv"), row.names=FALSE)
p2 <- ggplot(eapc, aes(eapc, sex_name, color=sex_name)) +
  geom_vline(xintercept=0, linetype=2, color="grey50") +
  geom_errorbarh(aes(xmin=lo, xmax=hi), height=0.12, linewidth=0.7) +
  geom_point(size=3) + scale_color_pub("npg") + guides(color="none") +
  labs(x="EAPC (%/year, 95% CI)", y=NULL,
       title="EAPC of ASR (log-linear), China") + theme_pub()
save_fig(p2, file.path(DAST,"eapc_forest"), width=4.6, height=2.6)

## ---- 4. 年龄-性别结构 (latest year, Number) → 背靠背 lollipop ---------------
cat("Step 4  年龄-性别结构 (背靠背 lollipop)\n")
yr_last <- max(burden$year)
ageN <- burden %>% filter(location_name=="China", year==yr_last, age_name %in% AGES, metric_name=="Rate") %>%
  left_join(popn %>% filter(location_name=="China", year==yr_last), by=c("location_name","sex_name","year","age_name")) %>%
  mutate(number = val.x/1e5 * val.y,
         age_f = factor(age_name, levels=AGES),
         signed = ifelse(sex_name=="Male", number, -number))
p3 <- ggplot(ageN, aes(signed, age_f, color=sex_name)) +
  geom_segment(aes(x=0, xend=signed, yend=age_f), linewidth=0.6) +
  geom_point(size=1.8) +
  geom_vline(xintercept=0, color="grey40") +
  scale_color_pub("npg") +
  scale_x_continuous(labels=function(x) format(abs(x), big.mark=",", scientific=FALSE)) +
  labs(x=paste0("Incident cases  (Female  |  Male), ", yr_last), y="Age group", color="Sex",
       title="Age-sex case structure, China") + theme_pub()
save_fig(p3, file.path(DAST,"agesex_lollipop"), width=5.0, height=4.4)

## ---- 5. Das Gupta 三因子分解 (老龄化/人口/流行病学), China 1990→2021 --------
## 公式逐字照搬 R-script-for-GBD/decomposition.R
cat("Step 5  Das Gupta 三因子分解 (公式照搬真实脚本)\n")
decomp <- list(); di <- 1
for(a in c("Male","Female")){
  pop90 <- popn %>% filter(sex_name==a, year==1990, location_name=="China") %>% arrange(match(age_name,AGES))
  pop21 <- popn %>% filter(sex_name==a, year==2021, location_name=="China") %>% arrange(match(age_name,AGES))
  P_1990 <- sum(pop90$val); P_2021 <- sum(pop21$val)
  a_1990 <- pop90$val/P_1990; a_2021 <- pop21$val/P_2021
  p_1990 <- P_1990; p_2021 <- P_2021
  r90 <- burden %>% filter(sex_name==a, year==1990, location_name=="China", metric_name=="Rate", age_name %in% AGES) %>% arrange(match(age_name,AGES))
  r21 <- burden %>% filter(sex_name==a, year==2021, location_name=="China", metric_name=="Rate", age_name %in% AGES) %>% arrange(match(age_name,AGES))
  r_1990 <- r90$val/1e5; r_2021 <- r21$val/1e5
  a_effect <- (sum(a_2021*p_1990*r_1990) + sum(a_2021*p_2021*r_2021))/3 +
              (sum(a_2021*p_1990*r_2021) + sum(a_2021*p_2021*r_1990))/6 -
              (sum(a_1990*p_1990*r_1990) + sum(a_1990*p_2021*r_2021))/3 -
              (sum(a_1990*p_1990*r_2021) + sum(a_1990*p_2021*r_1990))/6
  p_effect <- (sum(a_1990*p_2021*r_1990) + sum(a_2021*p_2021*r_2021))/3 +
              (sum(a_1990*p_2021*r_2021) + sum(a_2021*p_2021*r_1990))/6 -
              (sum(a_1990*p_1990*r_1990) + sum(a_2021*p_1990*r_2021))/3 -
              (sum(a_1990*p_1990*r_2021) + sum(a_2021*p_1990*r_1990))/6
  r_effect <- (sum(a_1990*p_1990*r_2021) + sum(a_2021*p_2021*r_2021))/3 +
              (sum(a_1990*p_2021*r_2021) + sum(a_2021*p_1990*r_2021))/6 -
              (sum(a_1990*p_1990*r_1990) + sum(a_2021*p_2021*r_1990))/3 -
              (sum(a_1990*p_2021*r_1990) + sum(a_2021*p_1990*r_1990))/6
  observed <- sum(a_2021*p_2021*r_2021) - sum(a_1990*p_1990*r_1990)   # All-ages Number 差
  decomp[[di]] <- data.frame(sex=a, Aging=a_effect, Population=p_effect, Epidemiology=r_effect,
                             total=a_effect+p_effect+r_effect, observed=observed); di <- di+1
}
decomp <- do.call(rbind, decomp)
# sanity: 三效应之和 ≈ 观测变化(真实脚本的 round(diff)==round(overll) 断言)
chk <- all(abs(decomp$total - decomp$observed) < 1e-6*pmax(1,abs(decomp$observed)))
cat(sprintf("   [check] sum(effects)==observed change: %s\n", chk))
write.csv(decomp, file.path(DRES,"dasgupta_decomposition.csv"), row.names=FALSE)
dlong <- decomp %>% select(sex, Aging, Population, Epidemiology) %>%
  tidyr::pivot_longer(-sex, names_to="factor", values_to="effect")
p4 <- ggplot(dlong, aes(effect, factor, color=effect>0)) +
  geom_vline(xintercept=0, color="grey40") +
  geom_segment(aes(x=0, xend=effect, yend=factor), linewidth=0.7) +
  geom_point(size=3) + facet_wrap(~sex) +
  scale_color_manual(values=c(`TRUE`="#BC3C29", `FALSE`="#0072B5"), guide="none") +
  scale_x_continuous(labels=function(x) format(x, big.mark=",", scientific=FALSE)) +
  labs(x="Contribution to change in cases (1990→2021)", y=NULL,
       title="Das Gupta decomposition, China") + theme_pub()
save_fig(p4, file.path(DAST,"decomposition_lollipop"), width=6.0, height=3.0)

## ---- 6. ASR–SDI 跨地区关联 (latest year) + Spearman + LOESS -----------------
cat("Step 6  ASR–SDI 关联 (Spearman + LOESS)\n")
sdi_asr <- burden %>% filter(age_name=="Age-standardized", metric_name=="Rate",
                             sex_name=="Female", year==yr_last) %>%
  select(location_name, val, sdi)
sp <- suppressWarnings(cor.test(sdi_asr$sdi, sdi_asr$val, method="spearman", exact=FALSE))
cat(sprintf("   Spearman rho=%.2f, p=%.3g\n", sp$estimate, sp$p.value))
p5 <- ggplot(sdi_asr, aes(sdi, val)) +
  geom_smooth(method="loess", se=TRUE, color="#0072B5", fill="#0072B5", alpha=0.15, linewidth=0.8) +
  geom_point(size=2.4, color="#BC3C29") +
  annotate("text", x=min(sdi_asr$sdi), y=max(sdi_asr$val), hjust=0, vjust=1,
           label=sprintf("Spearman rho = %.2f\np = %.3g", sp$estimate, sp$p.value), size=3) +
  labs(x="SDI", y="ASR (per 100,000)",
       title=paste0("ASR vs SDI across regions, ", yr_last)) + theme_pub()
save_fig(p5, file.path(DAST,"asr_sdi_scatter"), width=4.8, height=3.8)

cat("Done 527 · figures → assets/ , tables → results/\n")
