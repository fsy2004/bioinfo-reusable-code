# =============================================================================
# 593 · CASE — Cell-type-specific And Shared EQTL fine-mapping
# -----------------------------------------------------------------------------
# 用途   : 多细胞类型 eQTL 精细定位。把「跨细胞类型共享的 eQTL 效应」与
#          「只在某一种细胞类型里出现的特异效应」分开,同时解开 LD 混淆。
# 上游   : CASE (R 包) · Lin C, Lin Y, Li W, Xu L, Zhang X, Zhao H.
#          "Leveraging cell-type specificity and similarity improves single-cell
#          eQTL fine-mapping." Nat Commun 2026;17:5591.
#          doi:10.1038/s41467-026-72176-3 · PMID 42020412  (二者均已核实)
#          Repo: https://github.com/leaffur/CASE
#
# 本模块结构(库规矩:任何「更好」的方法都必须带一个本机能跑的朴素对照):
#   ① 基线 A「独立单细胞类型精细定位」 —— 每个细胞类型各跑一遍单效应贝叶斯回归
#      (Wakefield ABF / SuSiE 内部的 SER 步),完全不借用跨细胞类型信息。
#   ② 基线 B「完全共享(固定效应 meta)」 —— 先把各细胞类型的 z 合并成一个 meta z,
#      再定位,等于假设效应在所有细胞类型间完全一致。
#   ③ CASE 路径 —— 守卫式:装了 CASE 才跑,签名取自上游 man/CASE.Rd(见下)。
#   基线 A / B 是 CASE 想要超越的两个极端(全特异 vs 全共享),纯 base R 实现,
#   不装任何包即可跑完并出图。
#
# 依赖   : R base + ggplot2(框架 theme_pub.R)。CASE 为可选。
# 约定   : 图中文字英文,注释中文;相对路径;固定随机种子。
# =============================================================================

suppressWarnings(suppressMessages(library(ggplot2)))

# ---- 定位脚本目录 & 载入统一顶刊主题 ----------------------------------------
.this_dir <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  if (length(m)) return(dirname(normalizePath(sub("^--file=", "", a[m[1]]))))
  if (!is.null(sys.frames()[[1]]$ofile)) return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  getwd()
}
HERE <- .this_dir()
FRAMEWORK <- normalizePath(file.path(HERE, "..", "..", "..", "_framework", "theme_pub.R"),
                           mustWork = FALSE)
if (file.exists(FRAMEWORK)) {
  source(FRAMEWORK)
} else {
  stop("找不到框架文件 theme_pub.R:", FRAMEWORK)
}

# ---- 参数区(默认全部指向 example_data/,关键参数支持 --key value 覆盖)-----
opt <- bio_args(list(
  geno        = file.path(HERE, "example_data", "genotypes.csv"),   # 样本 × SNP 剂量
  expr        = file.path(HERE, "example_data", "expression.csv"),  # 样本 × 细胞类型 表达
  outdir      = file.path(HERE, "results"),
  assets      = file.path(HERE, "assets"),
  prior_var   = "0.04",    # 单效应先验方差 W(标准化尺度;0.2^2)
  coverage    = "0.95",    # 可信集覆盖度,对齐 CASE::get_credible_sets 默认 coverage_thres
  cor_min     = "0.5",     # 可信集纯度阈值,对齐 CASE::get_credible_sets 默认 cor.min
  ruled_out   = "1e-4",    # PIP 低于此值不进可信集,对齐 CASE::get_credible_sets 默认
  run_case    = FALSE,     # --run_case 时尝试真 CASE 路径(需已安装 CASE)
  seed        = "593"
))
PRIOR_VAR <- as.numeric(opt$prior_var)
COVERAGE  <- as.numeric(opt$coverage)
COR_MIN   <- as.numeric(opt$cor_min)
RULED_OUT <- as.numeric(opt$ruled_out)
SEED      <- as.integer(opt$seed)
set.seed(SEED)
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$assets, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Step 0 · 合成示例数据(缺失时自动生成,保证 turnkey)
# =============================================================================
make_example <- function(dir, n = 400, m = 200, C = 3, block = 20, maf = 0.3, seed = 593) {
  set.seed(seed)
  # LD:每 block 内 AR(1) 相关的两条单倍型,阈值化成 0/1 等位,相加得 0/1/2 剂量
  hap <- function() {
    H <- matrix(0, n, m)
    for (b in seq(1, m, by = block)) {
      idx <- b:min(b + block - 1, m); k <- length(idx)
      S <- 0.85 ^ abs(outer(seq_len(k), seq_len(k), "-"))   # AR(1) 相关阵
      L <- chol(S)
      H[, idx] <- matrix(rnorm(n * k), n, k) %*% L
    }
    (H > qnorm(1 - maf)) * 1
  }
  X <- hap() + hap()
  colnames(X) <- sprintf("rs%04d", seq_len(m))
  rownames(X) <- sprintf("S%03d", seq_len(n))

  Xs <- scale(X); Xs[is.na(Xs)] <- 0
  # 真值:idx_shared 在全部 3 个细胞类型有效应(shared);idx_spec 只在细胞类型 1-2(specific)
  idx_shared <- 30; idx_spec <- 150
  B <- matrix(0, m, C, dimnames = list(colnames(X), sprintf("CellType_%d", seq_len(C))))
  B[idx_shared, ] <- c(0.30, 0.28, 0.26)
  B[idx_spec, 1:2] <- c(0.26, 0.24)
  Y <- Xs %*% B + matrix(rnorm(n * C), n, C)
  colnames(Y) <- colnames(B); rownames(Y) <- rownames(X)

  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  hdr <- c("# synthetic, for demo only — 593 CASE module (seed 593)",
           "# 400 samples x 200 SNPs x 3 cell types; true eQTL: rs0030 (shared, 3 cell types), rs0150 (specific, cell types 1-2)")
  wr <- function(mat, path, id) {
    con <- file(path, "w"); writeLines(hdr, con); close(con)
    df <- data.frame(X = rownames(mat), mat, check.names = FALSE)
    names(df)[1] <- id
    suppressWarnings(utils::write.table(df, path, sep = ",", row.names = FALSE, quote = FALSE,
                                        append = TRUE, col.names = TRUE))
  }
  wr(X, file.path(dir, "genotypes.csv"), "sample")
  wr(Y, file.path(dir, "expression.csv"), "sample")
  truth <- data.frame(snp = c("rs0030", "rs0150"), index = c(idx_shared, idx_spec),
                      pattern = c("shared_all3", "specific_ct1_ct2"))
  utils::write.csv(truth, file.path(dir, "true_eqtl.csv"), row.names = FALSE)
  invisible(truth)
}

if (!file.exists(opt$geno) || !file.exists(opt$expr)) {
  cat("Step 0 · example_data 缺失,自动生成合成数据\n")
  make_example(file.path(HERE, "example_data"))
}

read_mat <- function(path) {
  df <- utils::read.csv(path, comment.char = "#", check.names = FALSE, stringsAsFactors = FALSE)
  rn <- df[[1]]; df <- df[, -1, drop = FALSE]
  M <- as.matrix(df); rownames(M) <- rn; M
}

cat("Step 1 · 读入基因型与分细胞类型表达\n")
X <- read_mat(opt$geno)
Y <- read_mat(opt$expr)
stopifnot(nrow(X) == nrow(Y))
N <- nrow(X); M <- ncol(X); C <- ncol(Y)
cat(sprintf("       N = %d samples · M = %d SNPs · C = %d cell types\n", N, M, C))

truth_file <- file.path(HERE, "example_data", "true_eqtl.csv")
truth <- if (file.exists(truth_file)) utils::read.csv(truth_file, stringsAsFactors = FALSE) else NULL

# =============================================================================
# Step 2 · 边际关联统计:每个 SNP × 每个细胞类型的 z 分数 + LD 矩阵 R
#          (与 CASE 官方 vignette 的输入构造一致:Z 为 M×C,R 为 M×M)
# =============================================================================
cat("Step 2 · 计算边际 z 分数矩阵与 LD 矩阵\n")
Xs <- scale(X); Xs[is.na(Xs)] <- 0
Ys <- scale(Y)
Z <- matrix(0, M, C, dimnames = list(colnames(X), colnames(Y)))
for (j in seq_len(M)) {
  for (c in seq_len(C)) {
    fit <- summary(stats::lm(Ys[, c] ~ Xs[, j]))
    Z[j, c] <- fit$coefficients[2, 3]     # t 值 ≈ z
  }
}
R <- stats::cor(Xs)
R[is.na(R)] <- 0

# =============================================================================
# Step 3 · 基线核心:单效应贝叶斯回归 (SER, Wakefield ABF) + 前向条件化
#   (a) SER:V = se^2 ≈ 1/N(标准化尺度),W = 先验效应方差
#       log BF_j = 0.5*log(V/(V+W)) + 0.5*z^2*W/(V+W)
#       alpha_j  = BF_j / sum_k BF_k   (SNP 先验均匀,假设该层只有一个因果位点)
#       这一层就是 SuSiE 内部的 SER 步。
#   (b) 前向条件化(stepwise conditional fine-mapping,GCTA-COJO 式经典做法):
#       第 1 层取 top SNP,把它放进回归作为协变量,重算所有 SNP 的条件 z,
#       再跑第 2 层 SER。因为我们有个体水平数据,条件 z 直接由 lm 精确算出,
#       不需要用 LD 近似。
#   (c) 汇总 PIP = 1 - prod_l (1 - alpha_l)  (SuSiE 的 PIP 定义)
#   逐细胞类型独立跑这套 = 「完全特异」极端(基线 A)。
# =============================================================================
ser_alpha <- function(z, N, W = PRIOR_VAR) {
  V <- 1 / N
  lbf <- 0.5 * log(V / (V + W)) + 0.5 * z^2 * W / (V + W)
  lbf <- lbf - max(lbf)
  bf <- exp(lbf)
  bf / sum(bf)
}

# 给定表型向量与条件集 S,返回所有 SNP 的条件 t 值(S 内的 SNP 记 0)
cond_z <- function(y, Xs, S = integer(0)) {
  m <- ncol(Xs); out <- numeric(m)
  if (!length(S)) {
    for (j in seq_len(m)) out[j] <- summary(stats::lm(y ~ Xs[, j]))$coefficients[2, 3]
  } else {
    Cov <- Xs[, S, drop = FALSE]
    for (j in seq_len(m)) {
      if (j %in% S) next
      cf <- summary(stats::lm(y ~ Cov + Xs[, j]))$coefficients
      out[j] <- cf[nrow(cf), 3]
    }
  }
  out
}

#' 前向条件化精细定位:返回 L 层的 alpha、可信集,以及汇总 PIP
#' @param z_of 函数(S) -> 长度 M 的条件 z 向量(封装单细胞类型或 meta 两种情形)
finemap_forward <- function(z_of, R, L = 2, N) {
  S <- integer(0); alphas <- list(); sets <- list()
  for (l in seq_len(L)) {
    z <- z_of(S)
    a <- ser_alpha(z, N = N)
    cs <- credible_set(a, R)
    alphas[[l]] <- a; sets[[l]] <- cs
    top <- which.max(a)
    S <- c(S, top)
  }
  pip <- 1 - Reduce(`*`, lapply(alphas, function(a) 1 - a))
  list(pip = pip, alphas = alphas, sets = sets, selected = S)
}

# 可信集:按 PIP 降序累加到 coverage,再用 cor.min 纯度过滤
# (阈值默认对齐上游 CASE::get_credible_sets 的 coverage_thres / cor.min / ruled_out)
credible_set <- function(pip, R, coverage = COVERAGE, cor.min = COR_MIN, ruled_out = RULED_OUT) {
  ord <- order(pip, decreasing = TRUE)
  keep <- ord[pip[ord] > ruled_out]
  if (!length(keep)) return(integer(0))
  cs <- keep[seq_len(which(cumsum(pip[keep]) >= coverage)[1])]
  if (is.na(cs[1])) cs <- keep
  if (length(cs) > 1) {
    purity <- min(abs(R[cs, cs]))
    if (purity < cor.min) return(integer(0))    # 纯度不足 → 不报告该可信集
  }
  sort(cs)
}

L_SIG <- 2   # 每个位点最多找几个独立信号(合成数据真值为 2)

cat("Step 3 · 基线 A:每个细胞类型独立精细定位(完全特异极端)\n")
fits_indep <- lapply(seq_len(C), function(c) {
  finemap_forward(function(S) cond_z(Ys[, c], Xs, S), R = R, L = L_SIG, N = N)
})
names(fits_indep) <- colnames(Z)
pip_indep <- sapply(fits_indep, `[[`, "pip")
dimnames(pip_indep) <- dimnames(Z)
cs_indep <- lapply(fits_indep, function(f) sort(unique(unlist(f$sets))))

cat("Step 4 · 基线 B:固定效应 meta 合并后定位(完全共享极端)\n")
# 每层都在同一条件集下把各细胞类型的条件 z 做固定效应合并(误差近似独立)
fit_meta <- finemap_forward(
  function(S) rowSums(sapply(seq_len(C), function(c) cond_z(Ys[, c], Xs, S))) / sqrt(C),
  R = R, L = L_SIG, N = N)
pip_meta_v <- fit_meta$pip
pip_meta <- matrix(pip_meta_v, M, C, dimnames = dimnames(Z))   # 共享 → 各细胞类型同一套 PIP
cs_meta <- sort(unique(unlist(fit_meta$sets)))

for (c in seq_len(C)) {
  sel <- rownames(Z)[fits_indep[[c]]$selected]
  cat(sprintf("       %-12s indep CS(total) = %3d SNPs | signals: %s\n",
              colnames(Z)[c], length(cs_indep[[c]]), paste(sel, collapse = ", ")))
}
cat(sprintf("       %-12s meta  CS(total) = %3d SNPs | signals: %s\n", "SHARED-meta",
            length(cs_meta), paste(rownames(Z)[fit_meta$selected], collapse = ", ")))

# =============================================================================
# Step 5 · CASE 路径(守卫式引用封装)
#   真实签名取自上游 man/CASE.Rd(2026-07 读取):
#     CASE(Z = NULL, R, hatB = NULL, hatS = NULL, N, V = NULL, cs = TRUE,
#          verbose = TRUE, ...)
#       Z    : M*C z 分数矩阵
#       R    : M*M LD 矩阵
#       N    : 长度 C 的样本量向量,或 C*C 矩阵(对角=样本量,非对角=重叠)
#       V    : (可选) C*C 细胞类型间噪声相关阵,默认单位阵
#     返回 "CASE" 对象:pi / U / V / pip (M*C) / post_mean (M*C);vignette 另用 fit$sets
#   其余导出函数(NAMESPACE):CASE_train, CASE_test, get_credible_sets
#     get_credible_sets(pips, R, verbose = TRUE, cor.min = 0.5,
#                       coverage_thres = 0.95, ruled_out = 1e-04)
#   ★ 未安装 CASE 时本步跳过并打印真实安装命令,绝不伪造结果。
# =============================================================================
case_res <- NULL
run_case_path <- function() {
  if (!requireNamespace("CASE", quietly = TRUE)) {
    cat("       CASE 未安装 → 跳过。安装:  devtools::install_github(\"leaffur/CASE\")\n")
    cat("       教程: https://github.com/leaffur/CASE  (vignette: Introduction_to_CASE)\n")
    return(NULL)
  }
  cat("       CASE 已安装,按 man/CASE.Rd 的签名调用\n")
  fit <- CASE::CASE(Z = Z, R = R, N = rep(N, C))
  list(pip = fit$pip, sets = fit$sets, pi = fit$pi, post_mean = fit$post_mean)
}
cat("Step 5 · CASE 路径(可选)\n")
if (isTRUE(opt$run_case) || identical(opt$run_case, "TRUE")) {
  case_res <- tryCatch(run_case_path(), error = function(e) {
    cat("       CASE 调用失败:", conditionMessage(e), "\n"); NULL })
} else {
  cat("       未加 --run_case,仅跑基线(CASE 结果不会被伪造)\n")
}

# =============================================================================
# Step 6 · 结果落盘
# =============================================================================
cat("Step 6 · 写出结果表\n")
long <- do.call(rbind, lapply(seq_len(C), function(c) data.frame(
  snp = rownames(Z), cell_type = colnames(Z)[c], z = Z[, c],
  pip_independent = pip_indep[, c], pip_shared_meta = pip_meta[, c],
  stringsAsFactors = FALSE)))
if (!is.null(case_res)) {
  cp <- as.matrix(case_res$pip)
  long$pip_CASE <- as.vector(cp[, seq_len(C)])
}
utils::write.csv(long, file.path(opt$outdir, "593_pip_table.csv"), row.names = FALSE)
utils::write.csv(Z, file.path(opt$outdir, "593_zscores.csv"))

cs_tab <- rbind(
  data.frame(method = "Independent (per cell type)", cell_type = names(cs_indep),
             cs_size = sapply(cs_indep, length),
             cs_snps = sapply(cs_indep, function(i) paste(rownames(Z)[i], collapse = ";")),
             stringsAsFactors = FALSE),
  data.frame(method = "Shared meta (fixed effect)", cell_type = "ALL",
             cs_size = length(cs_meta),
             cs_snps = paste(rownames(Z)[cs_meta], collapse = ";"), stringsAsFactors = FALSE))
utils::write.csv(cs_tab, file.path(opt$outdir, "593_credible_sets.csv"), row.names = FALSE)

# 真值处的 PIP(合成数据自带 ground truth,便于比较两个极端基线)
if (!is.null(truth)) {
  ev <- do.call(rbind, lapply(truth$snp, function(s) data.frame(
    snp = s, cell_type = colnames(Z),
    pip_independent = pip_indep[s, ], pip_shared_meta = pip_meta[s, ],
    stringsAsFactors = FALSE)))
  utils::write.csv(ev, file.path(opt$outdir, "593_truth_pip.csv"), row.names = FALSE)
  cat("       真 eQTL 处 PIP:\n")
  print(ev, row.names = FALSE)
}

# =============================================================================
# Step 7 · 出图(框架样式;不用条形图)
# =============================================================================
cat("Step 7 · 绘图\n")
snp_pos <- setNames(seq_len(M), rownames(Z))

# --- 图1:z 分数 scatter(按细胞类型分面),真 eQTL 圈出 -----------------------
d1 <- long[, c("snp", "cell_type", "z")]
d1$pos <- snp_pos[d1$snp]
d1$is_true <- FALSE
if (!is.null(truth)) {
  d1$is_true <- (d1$snp == "rs0030") |
                (d1$snp == "rs0150" & d1$cell_type %in% colnames(Z)[1:2])
}
bonf <- stats::qnorm(1 - 0.05 / (2 * M))
p1 <- ggplot(d1, aes(pos, z)) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.4) +
  geom_hline(yintercept = c(bonf, -bonf), linetype = "dashed",
             colour = pal_pub(name = "npg")[1], linewidth = 0.4) +
  geom_point(aes(colour = abs(z)), size = 1.5, alpha = 0.85) +
  geom_point(data = subset(d1, is_true), shape = 21, size = 4,
             colour = "black", fill = NA, stroke = 0.8) +
  scale_color_cont(name = "|z|") +
  facet_wrap(~ cell_type) +
  labs(x = "SNP position (index)", y = "Marginal z score",
       title = "Marginal eQTL association by cell type",
       subtitle = sprintf("Dashed line = Bonferroni threshold (|z| = %.2f); circles = true causal eQTL", bonf)) +
  theme_pub(base_size = 10)
save_fig(p1, file.path(opt$assets, "fig1_zscore_scatter"), width = 9, height = 3.6)

# --- 图2:PIP 热图(候选区域 × 细胞类型,两种基线并排)----------------------
top_idx <- order(pmax(apply(pip_indep, 1, max), pip_meta_v), decreasing = TRUE)[1:30]
top_snps <- rownames(Z)[sort(top_idx)]
d2 <- rbind(
  data.frame(snp = rep(top_snps, C), cell_type = rep(colnames(Z), each = length(top_snps)),
             pip = as.vector(pip_indep[top_snps, ]), method = "Independent per cell type"),
  data.frame(snp = rep(top_snps, C), cell_type = rep(colnames(Z), each = length(top_snps)),
             pip = as.vector(pip_meta[top_snps, ]), method = "Shared meta (fixed effect)"))
d2$snp <- factor(d2$snp, levels = top_snps)
p2 <- ggplot(d2, aes(cell_type, snp, fill = pip)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_cont(name = "PIP", limits = c(0, 1)) +
  facet_wrap(~ method) +
  labs(x = NULL, y = "SNP (top 30 by PIP)",
       title = "Posterior inclusion probability across cell types",
       subtitle = "Left: no information shared. Right: complete sharing forced. CASE models the middle ground.") +
  theme_pub(base_size = 9, border = TRUE) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_fig(p2, file.path(opt$assets, "fig2_pip_heatmap"), width = 7.5, height = 6.5)

# --- 图3:slopegraph —— 真 eQTL 处 PIP 在两种极端假设下的移动 ---------------
if (!is.null(truth)) {
  d3 <- do.call(rbind, lapply(truth$snp, function(s) data.frame(
    snp = s, cell_type = colnames(Z),
    Independent = pip_indep[s, ], `Shared meta` = pip_meta[s, ],
    check.names = FALSE, stringsAsFactors = FALSE)))
  d3l <- reshape(d3, direction = "long", varying = c("Independent", "Shared meta"),
                 v.names = "pip", timevar = "assumption",
                 times = c("Independent", "Shared meta"), idvar = c("snp", "cell_type"))
  d3l$assumption <- factor(d3l$assumption, levels = c("Independent", "Shared meta"))
  d3l$label <- paste0(d3l$snp, " · ",
                      ifelse(d3l$snp == "rs0030", "shared eQTL", "cell-type-specific eQTL"))
  p3 <- ggplot(d3l, aes(assumption, pip, group = interaction(snp, cell_type), colour = cell_type)) +
    geom_line(linewidth = 0.9, alpha = 0.85) +
    geom_point(size = 3) +
    scale_colour_manual(values = pal_pub(name = "npg"), name = "Cell type") +
    facet_wrap(~ label) +
    ylim(0, 1) +
    labs(x = NULL, y = "PIP at the true causal SNP",
         title = "What each sharing assumption costs you",
         subtitle = "Forcing complete sharing helps a truly shared eQTL and misleads a cell-type-specific one") +
    theme_pub(base_size = 10)
  save_fig(p3, file.path(opt$assets, "fig3_sharing_slopegraph"), width = 7.5, height = 4.2)
}

# --- 图4:可信集大小 dot plot(分辨率:越小越好)-----------------------------
d4 <- data.frame(
  cell_type = c(names(cs_indep), "Shared meta"),
  cs_size = c(sapply(cs_indep, length), length(cs_meta)),
  method = c(rep("Independent per cell type", C), "Shared meta (fixed effect)"),
  stringsAsFactors = FALSE)
d4$cell_type <- factor(d4$cell_type, levels = rev(d4$cell_type))
p4 <- ggplot(d4, aes(cs_size, cell_type, colour = method)) +
  geom_segment(aes(x = 0, xend = cs_size, yend = cell_type), linewidth = 0.5, colour = "grey75") +
  geom_point(size = 4.5) +
  scale_colour_manual(values = pal_pub(name = "npg")[c(4, 1)], name = NULL) +
  labs(x = "95% credible set size (SNPs, all signals pooled)", y = NULL,
       title = "Fine-mapping resolution",
       subtitle = "Smaller credible set = sharper localisation of the causal variant") +
  theme_pub(base_size = 10, legend = "bottom")
save_fig(p4, file.path(opt$assets, "fig4_credible_set_size"), width = 6.5, height = 3.6)

# ---- 依赖快照(可复现:记录 R 与包版本、随机种子)---------------------------
si <- utils::capture.output(utils::sessionInfo())
writeLines(c(sprintf("# 593 CASE module · seed = %d · prior_var = %g · L = %d", SEED, PRIOR_VAR, L_SIG),
             sprintf("# CASE installed: %s", requireNamespace("CASE", quietly = TRUE)), si),
           file.path(opt$outdir, "593_sessionInfo.txt"))

cat("Step 8 · 完成。结果 →", opt$outdir, " 图 →", opt$assets, "\n")
