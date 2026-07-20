# =============================================================================
# 595 · MR-EILLS — 不变性(invariance)稳健 MR,整合多个异质 GWAS summary 数据集
# -----------------------------------------------------------------------------
# 上游论文 : Hou L, Chen H, Zhou XH. MR-EILLS: an invariance-based Mendelian
#            randomization method integrating multiple heterogeneous GWAS summary
#            datasets. Nat Commun. 2025 Aug 18;16(1):7668.
#            doi:10.1038/s41467-025-62823-6 · PMID 40826001  (已核实,见 README)
# 上游代码 : https://github.com/hhoulei/MREILLS
#            R/MREILLS.R · R/BOOT_MREILLS.R · R/CHO_lambda.R
#            (本脚本的接口与损失函数逐行照读自上述 raw 源码,未臆造)
#
# 本模块做什么:
#   合成「E 个异质 GWAS summary 数据集 + 一部分无效工具(水平多效性)」→
#   跑 ①朴素基线 MVMR-IVW / MR-Egger ②MR-EILLS(装了官方包就调官方包,
#   没装就用逐行转写的本地实现)→ 对同一份数据比偏倚 → 出 4 张顶刊级图。
#
# ★ 基线是硬要求:任何「更稳健」的说法都必须和 IVW/Egger 在同一份数据上对照。
# ★ 图中文字英文;不用条形图(dot-whisker / violin+jitter / 折线 / 热图)。
# =============================================================================

suppressWarnings(suppressMessages({
  HERE <- NULL
}))

# ---- 定位与框架 --------------------------------------------------------------
.args_raw <- commandArgs(FALSE)
.m <- grep("^--file=", .args_raw)
HERE <- if (length(.m)) dirname(normalizePath(sub("^--file=", "", .args_raw[.m[1]]))) else getwd()
FRAMEWORK <- normalizePath(file.path(HERE, "..", "..", "..", "_framework", "theme_pub.R"),
                           mustWork = FALSE)
if (!file.exists(FRAMEWORK)) stop("找不到框架样式文件: ", FRAMEWORK)
source(FRAMEWORK)

set.seed(42)   # 固定种子

# ---- 参数区(默认指向 example_data/,关键参数支持 --key value 覆盖)----------
opt <- bio_args(list(
  input     = file.path(HERE, "example_data", "mreills_multi_gwas_summary.csv"),
  outdir    = file.path(HERE, "results"),
  assets    = file.path(HERE, "assets"),
  r1        = "0.1",          # 上游 MREILLS() 的 r1 (论文里的 gamma)
  lambda    = "auto",         # 上游 MREILLS() 的 lambda;auto = 数据驱动启发式(见 README)
  maxinvalid= "0.4",          # auto-lambda 假设的无效工具上限比例
  meth      = "L-BFGS-B",     # 上游 MREILLS() 的 optim 方法
  numboot   = "200",          # 上游 BOOT_MREILLS() 的 numBoot(200:SD 估计的相对误差 ~5%)
  start     = "ivw",          # optim 起点:ivw = 从 IVW 估计 warm start;zero = 上游的原点起步
  seed      = "42"
))
r1        <- as.numeric(opt$r1)
maxinvalid<- as.numeric(opt$maxinvalid)
numboot   <- as.integer(opt$numboot)
meth      <- as.character(opt$meth)
set.seed(as.integer(opt$seed))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$assets, recursive = TRUE, showWarnings = FALSE)

TRUE_THETA <- c(X1 = 0.40, X2 = 0.00)   # 合成数据的真因果效应(仅合成已知)

# =============================================================================
# Step 0 · 上游官方包入口(装上即用真 API;没装则走本地转写实现)
# -----------------------------------------------------------------------------
# 官方 API(读自 https://raw.githubusercontent.com/hhoulei/MREILLS/HEAD/R/*.R):
#   MREILLS(fdata_all, r1, meth, lambda)              -> optim 结果 list($par 为各暴露因果估计)
#   BOOT_MREILLS(fdata_all, r1, meth, lambda, numBoot)-> 各暴露估计的 bootstrap SE 向量
#   CHO_lambda(fdata_all)                             -> list(QSj=, plot=) 选 lambda 的岭线图
#   fdata_all = E 个数据集的 list,每个元素 list(betaGX, sebetaGX, betaGY, sebetaGY)
#     betaGX  : nSNP × nExposure 矩阵     sebetaGX: 同维标准误矩阵
#     betaGY  : 长度 nSNP 向量            sebetaGY: 同长标准误向量
# =============================================================================
HAS_UPSTREAM <- requireNamespace("MREILLS", quietly = TRUE)
cat("Step 0 · 上游官方包 MREILLS:", ifelse(HAS_UPSTREAM, "已安装 → 调用官方实现",
                                           "未安装 → 使用本地逐行转写实现"), "\n")
if (!HAS_UPSTREAM) {
  cat("        安装命令(需要联网 + devtools):devtools::install_github(\"hhoulei/MREILLS\")\n")
}

# ---- 本地转写:MR-EILLS 损失函数与估计 ---------------------------------------
# 逐行转写自上游 R/MREILLS.R 的 Qloss + optim。数学形式未改;唯一改动:
#   (a) 矩阵子集加 drop = FALSE,(b) 选中工具数 < Nv+2 时返回大损失,
#   防止 optim 探到退化解时因维度塌缩报错。两处均为数值稳健性护栏,非模型改动。
# 未安装官方包时无法逐位对拍,故本实现标注为「转写,未与官方包对拍验证」。
mreills_local <- function(fdata_all, r1, meth, lambda, start = NULL) {
  Nv <- ncol(fdata_all[[1]]$betaGX)
  Ng <- nrow(fdata_all[[1]]$betaGX)

  Qloss <- function(bb) {
    # --- 工具变量筛选:每个 SNP 的不变性得分 QSj,QSj < lambda 者入选 ---
    QSj <- numeric(Ng)
    for (sj in seq_len(Ng)) {
      QSj1 <- NULL; Xe <- NULL
      for (sje in seq_along(fdata_all)) {
        oncee <- fdata_all[[sje]]
        QSj1 <- c(QSj1, oncee$betaGY[sj] - matrix(oncee$betaGX[sj, ], nrow = 1) %*% bb)
        Xe   <- cbind(Xe, oncee$betaGX[sj, ])
      }
      Xe <- c(abs(Xe) %*% abs(QSj1))
      QSj[sj] <- sum(abs(QSj1)) + sum(Xe)
    }
    Sj <- (QSj < lambda)
    if (sum(Sj) < Nv + 2) return(1e6)      # 护栏(见上)

    # --- R:各数据集内的逆方差加权残差平方 ---
    Re <- unlist(lapply(fdata_all, function(x)
      mean((x$sebetaGY[Sj]^(-2) / sum(x$sebetaGY[Sj]^(-2))) *
             ((x$betaGY[Sj] - x$betaGX[Sj, , drop = FALSE] %*% bb)^2))))

    # --- J:不变性惩罚(残差与各暴露 betaGX 的相关,跨数据集加权求和)---
    Jpe <- NULL; weight <- NULL
    for (oe in seq_along(fdata_all)) {
      once <- fdata_all[[oe]]
      epi  <- once$betaGY[Sj] - once$betaGX[Sj, , drop = FALSE] %*% bb
      Jpe_once <- numeric(Nv)
      for (oj in seq_len(Nv)) Jpe_once[oj] <- mean(once$betaGX[Sj, oj] * epi)^2
      Jpe    <- cbind(Jpe, Jpe_once)
      weight <- c(weight, sum(once$sebetaGY[Sj]^(-2)))
    }
    weight <- weight / sum(weight)
    J <- sum(Jpe %*% weight)
    R <- sum(weight * Re)
    R + r1 * J
  }

  # 起点:上游写死 par = rep(0, Nv)。上游示例用 lambda = 100,在其数据尺度下
  # 相当于「不筛选」,损失面处处平滑,从 0 起步没问题。但当 lambda 收紧到真正起
  # 筛选作用时,bb = 0 处残差 = 全部因果信号 → QSj 普遍超阈 → 入选工具过少 →
  # 损失退化为常数,L-BFGS-B 的数值梯度为 0,估计会卡死在原点(本模块实测过)。
  # 因此默认用 IVW 估计 warm start(--start zero 可切回上游的原点起步)。
  if (is.null(start)) start <- rep(0, Nv)
  optim(par = start, fn = Qloss, method = meth)
}

# 逐行转写自上游 R/BOOT_MREILLS.R:对 SNP 有放回重抽,取估计的 sd 作为 SE
boot_mreills_local <- function(fdata_all, r1, meth, lambda, numBoot, engine, start = NULL) {
  estboot <- NULL
  for (ojo in seq_len(numBoot)) {
    locboot <- sample(seq_len(nrow(fdata_all[[1]]$betaGX)), replace = TRUE)
    fdata_boot <- lapply(fdata_all, function(f) list(
      betaGX   = f$betaGX[locboot, , drop = FALSE],
      sebetaGX = f$sebetaGX[locboot, , drop = FALSE],
      betaGY   = f$betaGY[locboot],
      sebetaGY = f$sebetaGY[locboot]))
    x.inv <- try(engine(fdata_boot, r1, meth, lambda, start), silent = TRUE)
    if (inherits(x.inv, "try-error")) next
    estboot <- cbind(estboot, x.inv$par)
  }
  apply(estboot, 1, function(x) stats::sd(x, na.rm = TRUE))
}

# 统一调度:有官方包走官方,否则走转写
# 注意:官方 MREILLS()/BOOT_MREILLS() 的形参里没有起点参数(起点写死 rep(0, Nv)),
# 所以走官方包时无法 warm start —— 这是官方 API 的事实,不做假封装。
.warned_no_warmstart <- FALSE
mreills_run  <- function(f, r1, meth, lambda, start = NULL) {
  if (HAS_UPSTREAM) {
    if (!is.null(start) && any(start != 0) && !.warned_no_warmstart) {
      cat("        [注意] 官方 MREILLS() 起点写死 rep(0, Nv),--start ivw 对官方路径无效;",
          "lambda 偏紧时估计可能卡在原点(见 README §②)\n")
      .warned_no_warmstart <<- TRUE
    }
    MREILLS::MREILLS(f, r1 = r1, meth = meth, lambda = lambda)
  } else mreills_local(f, r1, meth, lambda, start)
}
mreills_boot <- function(f, r1, meth, lambda, numBoot, start = NULL) {
  if (HAS_UPSTREAM) MREILLS::BOOT_MREILLS(f, r1 = r1, meth = meth, lambda = lambda,
                                          numBoot = numBoot)
  else boot_mreills_local(f, r1, meth, lambda, numBoot, mreills_local, start)
}
# 转写自上游 R/CHO_lambda.R:在 bb 处算各 SNP 的 QSj。两处明示偏离:
#   (a) 上游把 bb 写死为 rep(0, Nv),这里把 bb 提成参数(本模块要在 IVW / EILLS 解处取值);
#   (b) 上游按 rownames 取 SNP 的并集逐个匹配(允许各数据集 SNP 不同),这里按行位置索引 ——
#       等价的前提是各数据集 SNP 顺序一致,该前提在 Step 1 已用 stopifnot 强制检查。
qsj_score <- function(fdata_all, bb) {
  Ng <- nrow(fdata_all[[1]]$betaGX)
  vapply(seq_len(Ng), function(sj) {
    QSj1 <- NULL; QSj2 <- NULL
    for (je in seq_along(fdata_all)) {
      oncee <- fdata_all[[je]]
      epi1 <- c(oncee$betaGY[sj] - matrix(oncee$betaGX[sj, ], nrow = 1) %*% bb)
      QSj1 <- c(QSj1, abs(epi1))
      QSj2 <- c(QSj2, sum(abs(epi1 * oncee$betaGX[sj, ])))
    }
    sum(QSj1) + sum(QSj2)
  }, numeric(1))
}

# =============================================================================
# Step 1 · 合成示例数据(synthetic, for demo only)
# =============================================================================
make_synthetic <- function(path, g = 150, E = 3, p_invalid = 0.30) {
  rows <- list()
  gamma_base <- cbind(runif(g, 0.06, 0.30), runif(g, 0.06, 0.30))   # 共享的 G→X 强度
  invalid    <- rep(FALSE, g); invalid[seq_len(round(g * p_invalid))] <- TRUE
  alpha_base <- ifelse(invalid, runif(g, 0.02, 0.06), 0)            # 水平多效性(违反排他)
  for (e in seq_len(E)) {
    # 异质性:各数据集的 G→X 强度与多效性强度都不同(这正是 MR-EILLS 针对的场景)
    gamma_e <- gamma_base + matrix(rnorm(g * 2, 0, 0.02), ncol = 2)
    alpha_e <- alpha_base * runif(1, 0.8, 1.2) + ifelse(invalid, rnorm(g, 0, 0.005), 0)
    beta_Y  <- as.vector(gamma_e %*% TRUE_THETA) + alpha_e + rnorm(g, 0, 0.008)
    rows[[e]] <- data.frame(
      dataset  = paste0("GWAS_", e),
      SNP      = sprintf("rs%04d", seq_len(g)),
      beta_X1  = gamma_e[, 1], se_X1 = runif(g, 0.01, 0.04),
      beta_X2  = gamma_e[, 2], se_X2 = runif(g, 0.01, 0.04),
      beta_Y   = beta_Y,       se_Y  = runif(g, 0.01, 0.04),
      true_invalid_IV = as.integer(invalid),
      stringsAsFactors = FALSE)
  }
  df <- do.call(rbind, rows)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(c(
    "# synthetic, for demo only — 非真实 GWAS。3 个异质 GWAS summary 数据集,150 SNP,2 个暴露。",
    sprintf("# 真因果效应 theta = (X1=%.2f, X2=%.2f);前 30%% SNP 为无效工具(有水平多效性直接效应)。",
            TRUE_THETA[1], TRUE_THETA[2]),
    "# true_invalid_IV 仅合成数据可知,真实分析中不存在,只用于评估工具筛选表现。"),
    path)
  utils::write.table(df, path, sep = ",", row.names = FALSE, append = TRUE,
                     col.names = TRUE, quote = FALSE)
  invisible(df)
}

cat("Step 1 · 读取输入数据\n")
if (!file.exists(opt$input)) {
  cat("        输入不存在,合成示例数据 →", opt$input, "\n")
  make_synthetic(opt$input)
}
dat <- utils::read.csv(opt$input, comment.char = "#", stringsAsFactors = FALSE)
xcols <- grep("^beta_X", names(dat), value = TRUE)
scols <- sub("^beta_", "se_", xcols)
stopifnot(all(scols %in% names(dat)), all(c("beta_Y", "se_Y", "dataset", "SNP") %in% names(dat)))
exposures <- sub("^beta_", "", xcols)
datasets  <- unique(dat$dataset)
cat(sprintf("        %d 个数据集 × %d SNP × %d 个暴露:%s\n", length(datasets),
            sum(dat$dataset == datasets[1]), length(exposures),
            paste(exposures, collapse = ", ")))

# 组装上游要求的 fdata_all 结构
fdata_all <- lapply(datasets, function(d) {
  s <- dat[dat$dataset == d, ]
  bx <- as.matrix(s[, xcols, drop = FALSE]); rownames(bx) <- s$SNP; colnames(bx) <- exposures
  sx <- as.matrix(s[, scols, drop = FALSE]); rownames(sx) <- s$SNP
  by <- stats::setNames(s$beta_Y, s$SNP); sy <- stats::setNames(s$se_Y, s$SNP)
  list(betaGX = bx, sebetaGX = sx, betaGY = by, sebetaGY = sy)
})
# 上游 CHO_lambda() 按 rownames 匹配 SNP,本模块按行位置索引 → 必须先确认各数据集 SNP 完全对齐
snp_ref <- rownames(fdata_all[[1]]$betaGX)
stopifnot("各数据集必须是同一批 SNP 且顺序一致(见 README §①)" =
            all(vapply(fdata_all, function(f) identical(rownames(f$betaGX), snp_ref), logical(1))))

truth_invalid <- if ("true_invalid_IV" %in% names(dat))
  dat$true_invalid_IV[dat$dataset == datasets[1]] == 1 else NULL

# =============================================================================
# Step 2 · ★诚实基线:MVMR-IVW / MR-Egger(加权最小二乘,逆方差权重)
# =============================================================================
# 标准做法(Burgess et al.):以 1/se_Y^2 为权重,把 beta_Y 对 betaGX 做加权回归;
# 不含截距 = IVW,含截距 = MR-Egger(截距即定向多效性)。
wls_mr <- function(bx, by, sey, intercept = FALSE) {
  X <- if (intercept) cbind(`(Egger intercept)` = 1, bx) else bx
  w <- 1 / sey^2
  XtWX <- t(X) %*% (X * w)
  bhat <- solve(XtWX, t(X) %*% (by * w))
  res  <- by - X %*% bhat
  s2   <- max(sum(w * res^2) / (nrow(X) - ncol(X)), 1e-12)
  se   <- sqrt(diag(solve(XtWX)) * s2)
  data.frame(term = colnames(X), estimate = as.vector(bhat), se = as.vector(se),
             stringsAsFactors = FALSE)
}

cat("Step 2 · 基线:MVMR-IVW / MR-Egger\n")
BX_all  <- do.call(rbind, lapply(fdata_all, `[[`, "betaGX"))
BY_all  <- unlist(lapply(fdata_all, `[[`, "betaGY"))
SEY_all <- unlist(lapply(fdata_all, `[[`, "sebetaGY"))

ivw_pooled   <- wls_mr(BX_all, BY_all, SEY_all, intercept = FALSE)
egger_pooled <- wls_mr(BX_all, BY_all, SEY_all, intercept = TRUE)
per_ds <- do.call(rbind, lapply(seq_along(fdata_all), function(i) {
  r <- wls_mr(fdata_all[[i]]$betaGX, fdata_all[[i]]$betaGY, fdata_all[[i]]$sebetaGY)
  r$dataset <- datasets[i]; r
}))
# 固定效应逆方差 meta:先各数据集单独 IVW,再合并(常见的「多队列」朴素做法)
ivw_meta <- do.call(rbind, lapply(exposures, function(v) {
  s <- per_ds[per_ds$term == v, ]
  w <- 1 / s$se^2
  data.frame(term = v, estimate = sum(w * s$estimate) / sum(w), se = sqrt(1 / sum(w)))
}))
for (v in exposures) cat(sprintf("        IVW(pooled) %s = %+.3f   (truth %+.2f)\n",
                                 v, ivw_pooled$estimate[ivw_pooled$term == v], TRUE_THETA[[v]]))

# =============================================================================
# Step 3 · MR-EILLS
# =============================================================================
cat("Step 3 · MR-EILLS\n")
bb_ivw <- ivw_pooled$estimate[match(exposures, ivw_pooled$term)]
# ---- lambda 选择 ------------------------------------------------------------
# 上游 CHO_lambda() 只出 QSj 的岭线密度图,由使用者肉眼在谷底选 lambda;
# 这里给一个可复现的自动起点(非上游功能,README 已注明):在 IVW 估计处算 QSj,
# 取 (1 - maxinvalid) 分位数作为阈值,并额外输出 lambda 敏感性曲线供复核。
qsj_at_ivw <- qsj_score(fdata_all, bb_ivw)
if (identical(as.character(opt$lambda), "auto")) {
  lambda <- as.numeric(stats::quantile(qsj_at_ivw, 1 - maxinvalid))
  cat(sprintf("        lambda = auto → %.4f (QSj 的 %.0f%% 分位)\n", lambda, 100 * (1 - maxinvalid)))
} else {
  lambda <- as.numeric(opt$lambda)
  cat(sprintf("        lambda = %.4f (用户指定)\n", lambda))
}

start_par <- if (identical(as.character(opt$start), "zero")) rep(0, length(exposures)) else bb_ivw
cat(sprintf("        optim 起点 = %s (%s)\n", opt$start,
            paste(sprintf("%+.3f", start_par), collapse = ", ")))
fit <- mreills_run(fdata_all, r1 = r1, meth = meth, lambda = lambda, start = start_par)
eills_est <- as.vector(fit$par)
cat(sprintf("        optim convergence = %s\n", fit$convergence))
cat("        bootstrap SE(numBoot =", numboot, ") ...\n")
# ★ 重置种子:bootstrap 之前脚本消耗的随机数是路径依赖的(例如 example_data 缺失时
#   会先跑 make_synthetic() 抽一大批随机数),会把 RNG 流推到不同位置 → 同一份数据
#   两次运行得到不同的 bootstrap SE(实测 ±0.015 vs ±0.025)。在此固定独立种子,
#   使 SE 只取决于数据与超参,与「示例数据是否需要重建」无关。
set.seed(as.integer(opt$seed) + 1L)
eills_se <- mreills_boot(fdata_all, r1 = r1, meth = meth, lambda = lambda,
                         numBoot = numboot, start = start_par)
for (i in seq_along(exposures))
  cat(sprintf("        MR-EILLS %s = %+.3f ± %.3f   (truth %+.2f)\n",
              exposures[i], eills_est[i], eills_se[i], TRUE_THETA[[exposures[i]]]))

# 最终解处的工具筛选结果(用于诊断图)
qsj_final <- qsj_score(fdata_all, eills_est)
selected  <- qsj_final < lambda
cat(sprintf("        选中工具 %d / %d\n", sum(selected), length(selected)))
if (!is.null(truth_invalid)) {
  sens <- mean(!selected[truth_invalid])        # 无效工具被正确剔除的比例
  spec <- mean(selected[!truth_invalid])        # 有效工具被正确保留的比例
  cat(sprintf("        工具筛选:剔除无效 IV %.0f%% · 保留有效 IV %.0f%%\n", 100 * sens, 100 * spec))
}

# =============================================================================
# Step 4 · lambda 敏感性路径
# =============================================================================
cat("Step 4 · lambda 敏感性扫描\n")
lam_grid <- stats::quantile(qsj_at_ivw, seq(0.2, 0.98, length.out = 14))
path <- do.call(rbind, lapply(seq_along(lam_grid), function(k) {
  lam <- as.numeric(lam_grid[k])
  f <- try(mreills_run(fdata_all, r1 = r1, meth = meth, lambda = lam, start = start_par),
           silent = TRUE)
  if (inherits(f, "try-error")) return(NULL)
  data.frame(lambda = lam, quantile = as.numeric(sub("%", "", names(lam_grid)[k])) / 100,
             exposure = exposures, estimate = as.vector(f$par),
             n_selected = sum(qsj_score(fdata_all, as.vector(f$par)) < lam),
             stringsAsFactors = FALSE)
}))

# =============================================================================
# Step 5 · 汇总结果表
# =============================================================================
res <- rbind(
  data.frame(method = "IVW (pooled SNPs)", exposure = ivw_pooled$term,
             estimate = ivw_pooled$estimate, se = ivw_pooled$se),
  data.frame(method = "IVW (per-dataset meta)", exposure = ivw_meta$term,
             estimate = ivw_meta$estimate, se = ivw_meta$se),
  data.frame(method = "MR-Egger (pooled SNPs)", exposure = egger_pooled$term,
             estimate = egger_pooled$estimate, se = egger_pooled$se),
  data.frame(method = "MR-EILLS", exposure = exposures,
             estimate = eills_est, se = eills_se))
res$lo <- res$estimate - 1.96 * res$se
res$hi <- res$estimate + 1.96 * res$se
res$truth <- TRUE_THETA[res$exposure]
res$bias  <- res$estimate - res$truth
res$engine <- ifelse(res$method == "MR-EILLS",
                     ifelse(HAS_UPSTREAM, "MREILLS package", "local transcription"), "base R WLS")
utils::write.csv(res, file.path(opt$outdir, "MR_estimates_eills_vs_baseline.csv"), row.names = FALSE)
sel_df <- data.frame(SNP = rownames(fdata_all[[1]]$betaGX), QSj = qsj_final,
                     selected = selected,
                     true_invalid_IV = if (is.null(truth_invalid)) NA else truth_invalid)
utils::write.csv(sel_df, file.path(opt$outdir, "IV_selection_QSj.csv"), row.names = FALSE)
if (!is.null(path)) utils::write.csv(path, file.path(opt$outdir, "lambda_path.csv"), row.names = FALSE)

# =============================================================================
# Step 6 · 出图(顶刊风格;无条形图)
# =============================================================================
cat("Step 6 · 出图\n")
pal4 <- pal_pub(4, "npg")
# Egger 截距是多效性诊断量,不是因果效应,留在结果表里但不进因果效应图
egger_int <- res[res$exposure == "(Egger intercept)", ]
cat(sprintf("        MR-Egger 截距(定向多效性检验) = %+.4f ± %.4f\n",
            egger_int$estimate[1], egger_int$se[1]))
res <- res[res$exposure != "(Egger intercept)", ]
res$method <- factor(res$method, levels = c("IVW (pooled SNPs)", "IVW (per-dataset meta)",
                                            "MR-Egger (pooled SNPs)", "MR-EILLS"))
tru <- data.frame(exposure = names(TRUE_THETA), truth = as.numeric(TRUE_THETA))

# --- Fig 1 · dot-and-whisker:各方法估计 vs 真值 ---
p1 <- ggplot(res, aes(x = estimate, y = method, colour = method)) +
  geom_vline(data = tru, aes(xintercept = truth), linetype = "dashed",
             colour = "grey35", linewidth = 0.5) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0, linewidth = 0.7) +
  geom_point(size = 3) +
  facet_wrap(~ exposure, scales = "free_x",
             labeller = as_labeller(function(x) paste0("Exposure ", x,
               "  (true = ", sprintf("%.2f", TRUE_THETA[x]), ")"))) +
  scale_colour_manual(values = pal4, guide = "none") +
  labs(x = "Causal effect estimate (95% CI)", y = NULL,
       title = "MR-EILLS vs standard MR baselines",
       subtitle = "Synthetic data: 3 heterogeneous GWAS datasets, 30% invalid instruments") +
  theme_pub(11)
save_fig(p1, file.path(opt$assets, "fig1_estimates_dotwhisker"), width = 8.4, height = 3.6)

# --- Fig 2 · violin + jitter:不变性得分 QSj 按真实工具有效性分层 ---
if (!is.null(truth_invalid)) {
  d2 <- data.frame(QSj = qsj_final,
                   status = factor(ifelse(truth_invalid, "Invalid IV (pleiotropic)", "Valid IV"),
                                   levels = c("Valid IV", "Invalid IV (pleiotropic)")),
                   selected = selected)
  p2 <- ggplot(d2, aes(x = status, y = QSj)) +
    geom_violin(aes(fill = status), colour = NA, alpha = 0.30, width = 0.85, trim = FALSE) +
    geom_jitter(aes(colour = selected), width = 0.14, height = 0, size = 1.7, alpha = 0.85) +
    geom_hline(yintercept = lambda, linetype = "dashed", colour = "grey25", linewidth = 0.5) +
    annotate("text", x = 0.62, y = lambda, vjust = -0.6, hjust = 0, size = 3.2,
             colour = "grey25", label = sprintf("lambda = %.3f", lambda)) +
    scale_fill_manual(values = pal4[c(2, 1)], guide = "none") +
    scale_colour_manual(values = c(`TRUE` = pal4[3], `FALSE` = "grey60"),
                        labels = c(`TRUE` = "retained", `FALSE` = "screened out"),
                        name = "MR-EILLS") +
    labs(x = NULL, y = "Invariance score  QSj",
         title = "Instrument screening by the MR-EILLS invariance score",
         subtitle = "QSj evaluated at the MR-EILLS solution; SNPs with QSj < lambda enter the loss") +
    theme_pub(11)
  save_fig(p2, file.path(opt$assets, "fig2_iv_screening_violin"), width = 6.6, height = 4.4)
}

# --- Fig 3 · lambda 敏感性路径(折线 + 点) ---
if (!is.null(path)) {
  p3 <- ggplot(path, aes(x = lambda, y = estimate, colour = exposure)) +
    geom_hline(data = tru, aes(yintercept = truth, colour = exposure),
               linetype = "dashed", linewidth = 0.45, show.legend = FALSE) +
    geom_line(linewidth = 0.8) +
    geom_point(aes(size = n_selected), alpha = 0.9) +
    scale_colour_manual(values = pal4[c(1, 4)], name = "Exposure") +
    scale_size_continuous(range = c(1.2, 4.5), name = "IVs retained") +
    labs(x = "lambda (instrument-screening threshold)", y = "Causal effect estimate",
         title = "Sensitivity of MR-EILLS to the screening threshold",
         subtitle = "Dashed lines = true causal effects; point size = number of instruments retained") +
    theme_pub(11)
  save_fig(p3, file.path(opt$assets, "fig3_lambda_sensitivity"), width = 7.2, height = 4.4)
}

# --- Fig 4 · 偏倚热图(RdBu 发散) ---
p4 <- ggplot(res, aes(x = exposure, y = method, fill = bias)) +
  geom_tile(colour = "white", linewidth = 1.2) +
  geom_text(aes(label = sprintf("%+.3f", bias)), size = 3.4, colour = "black") +
  scale_fill_diverge(midpoint = 0, name = "Bias\n(est - truth)") +
  labs(x = "Exposure", y = NULL, title = "Signed bias against the known truth",
       subtitle = "Same synthetic data for every method") +
  theme_pub(11) + theme(axis.line = element_blank(), axis.ticks = element_blank())
save_fig(p4, file.path(opt$assets, "fig4_bias_heatmap"), width = 6.4, height = 3.6)

# =============================================================================
utils::capture.output(utils::sessionInfo(), file = file.path(opt$outdir, "sessionInfo.txt"))
cat("完成 · 结果 →", opt$outdir, " · 图 →", opt$assets, "\n")
