# =============================================================================
# 508 · 两步(网络)中介 MR — Sobel / Delta / Monte-Carlo + 中介比例
# -----------------------------------------------------------------------------
# 在孟德尔随机化框架下分解 暴露 X → 中介 M → 结局 Y 的中介效应:
#   Step1 用【X 的工具变量】IVW 估 α = X→M;
#   Step2 用【M 的工具变量】IVW 估 β = M→Y;
#   总效应 βT = X→Y(X 工具);间接(中介)效应 = α×β;直接效应 = βT − α×β。
# 间接效应显著性三法对照:Delta 法 SE、Sobel 检验、Monte-Carlo 置信区间(更稳健,
# 不依赖 α×β 正态假设);并给中介比例 PM = α×β/βT 的 MC 置信区间。
#
# 用 cis-pQTL/eQTL 作工具时即【药物靶点中介 MR(cis-MR)】的标准做法。
# ★诚实措辞:中介比例对弱总效应不稳;Delta/Sobel 假设正态,首选 MC-CI;
#   直接效应 SE 为近似(βT 与间接共享 X 工具),严谨场景请用 MVMR 估直接效应。
#
# Turnkey: Rscript 508_twostep_mediation_mr.R   (默认合成工具汇总统计→results/+assets/)
#          换数据: --x_instruments x.csv --m_instruments m.csv  (列见 README)
# 复用 _framework/theme_pub.R;图中文字英文,注释中文。无条形图(森林+路径图)。
# =============================================================================
suppressWarnings(suppressMessages({ library(ggplot2) }))

## ---- 定位脚本目录 + 载入框架(向上搜 _framework) --------------------------
.this <- tryCatch({ a <- commandArgs(FALSE); m <- grep("^--file=", a)
  if(length(m)) dirname(normalizePath(sub("^--file=","",a[m[1]]))) else getwd() },
  error=function(e) getwd())
.fw <- NULL; .p <- .this
for(i in 1:6){ cand <- file.path(.p,"_framework","theme_pub.R"); if(file.exists(cand)){ .fw <- cand; break }; .p <- dirname(.p) }
if(!is.null(.fw)) source(.fw) else {
  theme_pub <- function(base_size=11, ...) theme_bw(base_size=base_size)
  pal_pub   <- function(n=NULL, name="npg") scales::hue_pal()(ifelse(is.null(n),6,n))
  save_fig  <- function(plot, file, width=7, height=6, dpi=300){
    ggsave(paste0(file,".pdf"), plot, width=width, height=height)
    ggsave(paste0(file,".png"), plot, width=width, height=height, dpi=dpi) }
}

set.seed(42)
DIR <- .this
DDAT <- file.path(DIR,"example_data"); DRES <- file.path(DIR,"results"); DAST <- file.path(DIR,"assets")
for(d in c(DDAT,DRES,DAST)) dir.create(d, showWarnings=FALSE, recursive=TRUE)

## ---- IVW(固定效应):outcome-beta ~ exposure-beta 过原点,权重 1/se_out^2 ----
ivw <- function(b_exp, b_out, se_out){
  w <- 1/se_out^2
  est <- sum(b_exp*b_out*w)/sum(b_exp^2*w)
  se  <- sqrt(1/sum(b_exp^2*w))
  c(est=est, se=se)
}

## ---- 1. 读 / 生成合成工具汇总统计 ------------------------------------------
# 合成两组工具:X 工具(估 α 与 βT)、M 工具(估 β)。真值:α=0.40, β=0.50,
# 直接 βd=0.20 → 总 βT=0.20+0.40*0.50=0.40,真中介比例=0.20/0.40=50%。
fX <- file.path(DDAT,"x_instruments.csv"); fM <- file.path(DDAT,"m_instruments.csv")
if(!file.exists(fX)){
  aT<-0.40; bT<-0.50; bd<-0.20; betaT<-bd+aT*bT
  mk_se <- function(n) runif(n,0.02,0.05)
  nX<-30; gX<-abs(rnorm(nX,0.12,0.04)); seX<-mk_se(nX)
  # X 工具对 M、Y 的效应(含测量噪声)
  gM_x <- aT*gX + rnorm(nX,0,0.01);    se_gM_x<-mk_se(nX)
  gY_x <- betaT*gX + rnorm(nX,0,0.01); se_gY_x<-mk_se(nX)
  dX <- data.frame(SNP=sprintf("rsX%02d",1:nX),
                   beta_exposure=gX, se_exposure=seX,
                   beta_mediator=gM_x, se_mediator=se_gM_x,
                   beta_outcome=gY_x, se_outcome=se_gY_x)
  nM<-30; gM<-abs(rnorm(nM,0.12,0.04)); seM<-mk_se(nM)
  gY_m <- bT*gM + rnorm(nM,0,0.01); se_gY_m<-mk_se(nM)
  dM <- data.frame(SNP=sprintf("rsM%02d",1:nM),
                   beta_mediator=gM, se_mediator=seM,
                   beta_outcome=gY_m, se_outcome=se_gY_m)
  write.csv(dX,fX,row.names=FALSE); write.csv(dM,fM,row.names=FALSE)
  cat("[gen] synthetic instruments: X(n=30) + M(n=30); true PM=50% (for demo only)\n")
}
dX <- read.csv(fX); dM <- read.csv(fM)

## ---- 2. 两步 MR 估计 -------------------------------------------------------
# Step1 α = X→M(X 工具);Step2 β = M→Y(M 工具);总 βT = X→Y(X 工具)
A  <- ivw(dX$beta_exposure, dX$beta_mediator, dX$se_mediator)   # α
B  <- ivw(dM$beta_mediator, dM$beta_outcome, dM$se_outcome)     # β
TT <- ivw(dX$beta_exposure, dX$beta_outcome, dX$se_outcome)     # βT
alpha<-A["est"]; se_a<-A["se"]; beta<-B["est"]; se_b<-B["se"]; bT<-TT["est"]; se_T<-TT["se"]

## ---- 3. 间接效应 + 三法显著性 + 中介比例 -----------------------------------
ind <- alpha*beta                                  # 间接(中介)效应 = α×β
delta_se <- sqrt(alpha^2*se_b^2 + beta^2*se_a^2)   # Delta 法 SE(product of coefficients)
sobel_z  <- ind/delta_se; sobel_p <- 2*pnorm(-abs(sobel_z))
delta_ci <- c(ind-1.96*delta_se, ind+1.96*delta_se)
# Monte-Carlo:对 α、β、βT 抽样,得间接效应与中介比例的经验区间(不依赖正态乘积)
NMC<-2e5
sa<-rnorm(NMC,alpha,se_a); sb<-rnorm(NMC,beta,se_b); sT<-rnorm(NMC,bT,se_T)
mc_ind <- sa*sb; mc_pm <- (sa*sb)/sT
mc_ci  <- quantile(mc_ind, c(.025,.975));
direct <- bT-ind; se_direct <- sqrt(se_T^2+delta_se^2)   # 近似(见 README 警示)
direct_ci <- c(direct-1.96*se_direct, direct+1.96*se_direct)
pm <- ind/bT; pm_ci <- quantile(mc_pm[is.finite(mc_pm)], c(.025,.975))

res <- data.frame(
  effect = c("Total (X->Y)","Direct (X->Y | M)","Indirect (X->M->Y)"),
  estimate = c(bT, direct, ind),
  lci = c(bT-1.96*se_T, direct_ci[1], mc_ci[1]),
  uci = c(bT+1.96*se_T, direct_ci[2], mc_ci[2]))
write.csv(res, file.path(DRES,"mediation_effects.csv"), row.names=FALSE)
write.csv(data.frame(
  metric=c("alpha_X_to_M","beta_M_to_Y","indirect_alpha*beta","delta_SE","Sobel_z","Sobel_p",
           "MC_indirect_LCI","MC_indirect_UCI","proportion_mediated","PM_LCI","PM_UCI"),
  value=round(c(alpha,beta,ind,delta_se,sobel_z,sobel_p,mc_ci[1],mc_ci[2],pm,pm_ci[1],pm_ci[2]),4)),
  file.path(DRES,"mediation_stats.csv"), row.names=FALSE)
cat(sprintf("[mr] alpha(X->M)=%.3f  beta(M->Y)=%.3f  total=%.3f\n", alpha, beta, bT))
cat(sprintf("[mediation] indirect=%.3f  Sobel p=%.3g  MC 95%%CI=[%.3f,%.3f]\n", ind, sobel_p, mc_ci[1], mc_ci[2]))
cat(sprintf("[mediation] proportion mediated=%.1f%%  MC 95%%CI=[%.1f%%,%.1f%%]\n", 100*pm, 100*pm_ci[1], 100*pm_ci[2]))

## ---- 4. 出图(路径图 + 森林图;无条形图)----------------------------------
col <- pal_pub(3,"npg")
# Fig1: 中介路径图(X→M→Y + 直接边;箭上标效应)
nodes <- data.frame(x=c(1,3,5), y=c(1,2.4,1), lab=c("Exposure (X)","Mediator (M)","Outcome (Y)"))
edges <- data.frame(
  x=c(1.25,3,1.2), y=c(1.2,2.3,0.9), xend=c(2.78,4.78,4.8), yend=c(2.25,1.2,0.9),
  lab=c(sprintf("alpha = %.2f", alpha), sprintf("beta = %.2f", beta),
        sprintf("direct = %.2f", direct)),
  lx=c(1.7,4.25,3), ly=c(1.95,1.95,0.74),
  kind=c("path","path","direct"))
p_path <- ggplot() +
  geom_segment(data=edges, aes(x=x,y=y,xend=xend,yend=yend, color=kind),
               linewidth=1.1, arrow=arrow(length=unit(0.18,"cm"), type="closed")) +
  geom_text(data=edges, aes(lx,ly,label=lab), size=3.4, fontface="italic") +
  geom_label(data=nodes, aes(x,y,label=lab), fill="white", linewidth=0.5,
             fontface="bold", size=3.7) +
  scale_color_manual(values=c(path=col[2], direct="grey55"), guide="none") +
  coord_cartesian(xlim=c(0.4,5.7), ylim=c(0.4,2.9)) +
  labs(title=sprintf("Two-step mediation MR  (proportion mediated = %.0f%%)", 100*pm)) +
  theme_void(base_size=12) +
  theme(plot.title=element_text(face="bold", hjust=0.5, size=12))
save_fig(p_path, file.path(DAST,"mediation_path"), width=6.4, height=3.2)

# Fig2: 效应森林图(总/直接/间接,点+95%CI)
res$effect <- factor(res$effect, levels=rev(res$effect))
p_for <- ggplot(res, aes(estimate, effect, color=effect)) +
  geom_vline(xintercept=0, linetype=2, color="grey60") +
  geom_errorbar(aes(xmin=lci, xmax=uci), orientation="y", width=0.18, linewidth=0.8) +
  geom_point(size=3.4) +
  scale_color_manual(values=rev(col), guide="none") +
  labs(x="MR effect estimate (95% CI)", y=NULL,
       title="Total / direct / indirect effects") +
  theme_pub(base_size=11)
save_fig(p_for, file.path(DAST,"effects_forest"), width=6, height=2.8)

# Fig3: 间接效应三法一致性(Delta vs Sobel-implied vs Monte-Carlo)森林对照
meth <- data.frame(
  method=c("Delta (normal)","Sobel (normal)","Monte-Carlo"),
  est=c(ind, ind, mean(mc_ind)),
  lci=c(delta_ci[1], ind-1.96*delta_se, mc_ci[1]),
  uci=c(delta_ci[2], ind+1.96*delta_se, mc_ci[2]))
meth$method <- factor(meth$method, levels=rev(meth$method))
p_cmp <- ggplot(meth, aes(est, method)) +
  geom_vline(xintercept=0, linetype=2, color="grey60") +
  geom_errorbar(aes(xmin=lci, xmax=uci), orientation="y", width=0.16, linewidth=0.8, color=col[1]) +
  geom_point(size=3.2, color=col[1]) +
  labs(x="Indirect effect (95% CI)", y=NULL,
       title="Indirect-effect CI by method (MC most robust)") +
  theme_pub(base_size=11)
save_fig(p_cmp, file.path(DAST,"indirect_methods"), width=6, height=2.4)

cat("[fig] assets/: mediation_path, effects_forest, indirect_methods (.pdf+.png)\n")
sink(file.path(DRES,"sessionInfo.txt")); print(sessionInfo()); sink()
