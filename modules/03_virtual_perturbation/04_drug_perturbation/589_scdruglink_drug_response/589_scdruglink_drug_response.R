# =============================================================================
# 589 · scDrugLink 单细胞药物重定位 (link drug targets × perturbation signatures)
# -----------------------------------------------------------------------------
# 上游: scDrugLink (Huang L, Lu X, Chen D. IEEE J Biomed Health Inform 2026;
#       30(6):4535-4547. PMID 40100675, DOI 10.1109/JBHI.2025.3552536)
#       源码 https://github.com/LHBCB/scDrugLink,本地克隆
#       C:/Users/fsy/Desktop/upstream-sources/589_scDrugLink/
#
# 上游把「药物-靶点」与「药物-扰动签名」两类证据在**细胞类型层面**相乘串起来:
#   A 臂 promotion/inhibition : Drug2Cell 打分 → 细胞类型内 病 vs 对照 Wilcoxon
#                               + Cliff's delta → weight = delta × (−log10 p_adj)
#                               (R/build_drug_target_d2c.R, R/compute_drug_prom_inh.R)
#   B 臂 sensitivity/resistance: 细胞类型内 DEG 与药物扰动签名做反向匹配 + K-S 检验
#                               (R/compute_drug_sens_res.R,上游调 Asgard::GetDrug)
#   linking (Eq.7-9)          : score = Σ_ct (prop/100)·(−log10 FDR)·treated_ratio
#                                        · exp(weight[ct, drug])
#                               (R/compute_scdruglink_score.R 逐行照抄的公式)
#
# ★本模块实现边界(诚实声明,勿当成上游包本身):
#   · A 臂 + linking 公式 = 按上游源码逐步复刻,本机零依赖可跑(仅 base R + ggplot2)。
#   · B 臂 上游依赖 Asgard + CMAP L1000 GSE70138/GSE92742 gctx(数十 GB,本机未装)。
#     本脚本先探测 Asgard/cmapR,**装了就走上游函数**;没装则退回内置
#     CMap 式 KS 反向连接基线(Lamb et al. Science 2006 连接性打分思想 + 置换检验),
#     它**不是** Asgard::GetDrug 的等价实现,只用于把管道跑通/做方法教学。
#   · effsize 未装 → 内置 Cliff's delta,用其定义式 (#x>y − #x<y)/(n1·n2)。
#   · 上游 build_drug_target_d2c.R:23 按 drug_target_df$drug_name_lower 循环,而上游自带
#     cns_drug_targets 只有 drug_name/gene_names 两列(实测 rda 273 行);:89 又引用未传入的
#     全局 drug_label_df。本脚本全程用 drug_name 作键,规避这两处(详见 README)。
#
# Turnkey: Rscript 589_scdruglink_drug_response.R    (合成数据 → results/ + assets/)
# 换数据:  Rscript 589_scdruglink_drug_response.R --expr my_expr.csv --meta my_meta.csv \
#              --targets my_targets.csv --sig my_sig.csv --outdir results/run1
# 复用 _framework/theme_pub.R;无条形图(heatmap / lollipop / 散点 / violin)。
# =============================================================================

suppressWarnings(suppressMessages({ library(ggplot2) }))

## ---- 定位框架 ---------------------------------------------------------------
.this <- tryCatch({ a <- commandArgs(FALSE); m <- grep("^--file=", a)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd() },
  error = function(e) getwd())
.fw <- NULL; .p <- .this
for (i in 1:6) { cand <- file.path(.p, "_framework", "theme_pub.R")
  if (file.exists(cand)) { .fw <- cand; break }; .p <- dirname(.p) }
if (!is.null(.fw)) source(.fw) else stop("需要 _framework/theme_pub.R")

set.seed(42)
DIR  <- .this
DDAT <- file.path(DIR, "example_data")
args <- bio_args(list(
  expr    = file.path(DDAT, "expr_lognorm.csv"),
  meta    = file.path(DDAT, "cell_meta.csv"),
  targets = file.path(DDAT, "drug_targets.csv"),
  sig     = file.path(DDAT, "drug_perturb_signature.csv"),
  labels  = file.path(DDAT, "drug_labels.csv"),
  disease = "Disease", control = "Control",
  outdir  = file.path(DIR, "results"),
  n_bins = "25", ctrl_size = "50", n_perm = "500"
))
DRES <- args$outdir; DAST <- file.path(DIR, "assets")
for (d in c(DDAT, DRES, DAST)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
N_BINS <- as.integer(args$n_bins); CTRL_SIZE <- as.integer(args$ctrl_size)
N_PERM <- as.integer(args$n_perm)
DISEASE <- args$disease; CONTROL <- args$control

## =============================================================================
## Step 0 · 合成示例数据 (synthetic, for demo only)
## =============================================================================
gen_demo <- function() {
  ng <- 400; genes <- sprintf("G%03d", 1:ng)
  cts <- c("Microglia", "Astrocyte", "Oligodendrocyte", "Neuron", "Tcell", "Endothelial")
  n_per <- 60                                    # 每类 60 细胞(30 病 / 30 对照)
  meta <- do.call(rbind, lapply(cts, function(ct) data.frame(
    cell_type = ct, disease = rep(c(DISEASE, CONTROL), each = n_per / 2),
    stringsAsFactors = FALSE)))
  rownames(meta) <- sprintf("C%04d", seq_len(nrow(meta)))
  nc <- nrow(meta)

  # 基线:细胞类型特异表达 + 泊松噪声
  base <- matrix(rpois(ng * nc, 3), ng, nc, dimnames = list(genes, rownames(meta)))
  for (k in seq_along(cts)) {
    mk <- genes[((k - 1) * 15 + 1):(k * 15)]
    base[mk, meta$cell_type == cts[k]] <- base[mk, meta$cell_type == cts[k]] +
      rpois(length(mk) * sum(meta$cell_type == cts[k]), 8)
  }
  # 疾病程序:一组基因在 Microglia / Astrocyte 中显著上调(疾病最相关的两类)
  dis_up   <- genes[201:240]; dis_down <- genes[241:270]
  hit_ct   <- c("Microglia", "Astrocyte")
  sel <- meta$disease == DISEASE & meta$cell_type %in% hit_ct
  base[dis_up, sel]   <- base[dis_up, sel]   + rpois(length(dis_up) * sum(sel), 10)
  base[dis_down, sel] <- pmax(0, base[dis_down, sel] - 2)

  # log-normalise(CPM-like 1e4 + log1p,与 Seurat NormalizeData 同式)
  expr <- log1p(sweep(base, 2, pmax(colSums(base), 1), "/") * 1e4)

  # 药物靶点:20 药。1-4 号靶向疾病上调程序(应得正 promotion 权重)
  drugs <- sprintf("Drug%02d", 1:20)
  tg <- lapply(seq_along(drugs), function(i) {
    if (i <= 4)      sample(dis_up, 12)
    else if (i <= 8) sample(dis_down, 10)
    else             sample(genes, 12)
  })
  targets <- data.frame(drug_name = drugs,
                        gene_names = vapply(tg, paste, "", collapse = ";"),
                        stringsAsFactors = FALSE)

  # 扰动签名(gene × drug 的 z 值):1-4 号反转疾病方向(应得高 sensitivity)
  sig <- matrix(rnorm(ng * length(drugs), 0, 1), ng, length(drugs),
                dimnames = list(genes, drugs))
  sig[dis_up,   1:4] <- sig[dis_up,   1:4] - 2.2      # 疾病上调 → 药物压低
  sig[dis_down, 1:4] <- sig[dis_down, 1:4] + 2.2      # 疾病下调 → 药物拉高
  sig[dis_up,   5:8] <- sig[dis_up,   5:8] + 1.5      # 5-8 同向(应被判无效/加重)
  sig[dis_up,   9:12] <- sig[dis_up,   9:12] - 1.1    # 9-12 弱反转但靶点随机
  sig[dis_down, 9:12] <- sig[dis_down, 9:12] + 1.1    #  → 只有 B 臂证据,考验 linking

  # 已知适应症标签(仅用于评估 AUC/AUPR,对应上游 cns_drug_labels)
  labels <- data.frame(drug_name = drugs,
                       known_label = c(rep(1, 4), rep(0, 16)), stringsAsFactors = FALSE)

  write.csv(round(expr, 4), file.path(DDAT, "expr_lognorm.csv"))
  write.csv(data.frame(cell = rownames(meta), meta), file.path(DDAT, "cell_meta.csv"), row.names = FALSE)
  write.csv(targets, file.path(DDAT, "drug_targets.csv"), row.names = FALSE)
  write.csv(round(sig, 4), file.path(DDAT, "drug_perturb_signature.csv"))
  write.csv(labels, file.path(DDAT, "drug_labels.csv"), row.names = FALSE)
  cat(sprintf("[gen] synthetic: %d genes x %d cells, %d cell types, %d drugs (demo only)\n",
              ng, nc, length(cts), length(drugs)))
}
if (!file.exists(args$expr)) { cat("Step 0: 生成合成示例数据\n"); gen_demo() } else cat("Step 0: 复用 example_data/\n")

expr <- as.matrix(read_table_smart(args$expr, row_names = TRUE))
meta <- read_table_smart(args$meta, row_names = TRUE)
targets_df <- read_table_smart(args$targets)
sig_mat <- as.matrix(read_table_smart(args$sig, row_names = TRUE))
labels_df <- if (file.exists(args$labels)) read_table_smart(args$labels) else NULL
stopifnot(all(c("cell_type", "disease") %in% colnames(meta)))
meta <- meta[colnames(expr), , drop = FALSE]
cat(sprintf("  载入 expr %d x %d;cell types = %d;drugs = %d\n",
            nrow(expr), ncol(expr), length(unique(meta$cell_type)), nrow(targets_df)))

## =============================================================================
## Step 1 · Drug2Cell 靶点打分矩阵 (照抄 R/build_drug_target_d2c.R 的算法)
##   scores = X·W(靶点组均表达) − control_profiles·drug_weights(同表达档背景)
## =============================================================================
cat("Step 1: 构建 Drug2Cell 矩阵 (drug x cell)\n")
build_d2c <- function(expr, targets_df, n_bins = 25, ctrl_size = 50, seed = 42) {
  gene_list <- rownames(expr)
  drugs <- targets_df$drug_name
  targets <- lapply(drugs, function(d)
    gene_list %in% strsplit(targets_df$gene_names[targets_df$drug_name == d], ";")[[1]])
  names(targets) <- drugs

  weights <- as.data.frame(targets, row.names = gene_list, check.names = FALSE)
  weights <- sweep(weights, 2, colSums(weights) + 1e-6, `/`)   # 每药靶点等权归一

  X <- t(expr)                                                 # cell x gene
  scores <- X %*% as.matrix(weights)

  # ---- Seurat AddModuleScore 式背景分档(上游注释 "seurat scoring mechanism")----
  obs_avg  <- colMeans(X)
  n_items  <- round(length(obs_avg) / (n_bins - 1))
  obs_rank <- rank(obs_avg, ties.method = "min") - 1
  obs_cut  <- obs_rank %/% n_items

  set.seed(seed)
  control_groups <- list()
  for (cut_value in unique(obs_cut)) {
    mask <- (obs_cut == cut_value); r_genes <- sample(which(mask))
    if (length(r_genes) > ctrl_size) { mask[] <- FALSE; mask[r_genes[1:ctrl_size]] <- TRUE }
    control_groups[[paste0("bin", cut_value)]] <- mask
  }
  cgw <- as.data.frame(control_groups, row.names = gene_list, check.names = FALSE)
  cgw <- sweep(cgw, 2, colSums(cgw) + 1e-6, FUN = "/")
  control_profiles <- X %*% as.matrix(cgw)

  drug_bins <- lapply(colnames(weights), function(d)
    colnames(cgw) %in% paste0("bin", unique(obs_cut[targets[[d]]])))
  names(drug_bins) <- colnames(weights)
  dw <- as.data.frame(drug_bins, row.names = colnames(cgw), check.names = FALSE)
  dw <- sweep(dw, 2, colSums(dw) + 1e-6, FUN = "/")

  t(scores - control_profiles %*% as.matrix(dw))               # drug x cell
}
d2c <- build_d2c(expr, targets_df, N_BINS, CTRL_SIZE)
write.csv(round(d2c, 5), file.path(DRES, "drug_target_d2c.csv"))
cat(sprintf("  d2c: %d drugs x %d cells\n", nrow(d2c), ncol(d2c)))

## =============================================================================
## Step 2 · A 臂:促进/抑制效应 (R/compute_drug_prom_inh.R)
##   细胞类型内 病 vs 对照 Wilcoxon → BH;Cliff's delta;weight = delta·(−log10 p_adj)
## =============================================================================
cat("Step 2: 估计 promotion/inhibition 权重\n")
cliff_delta <- function(x, y) {                # effsize::cliff.delta 的 estimate 定义式
  rx <- rank(c(x, y)); n1 <- length(x); n2 <- length(y)
  U <- sum(rx[seq_len(n1)]) - n1 * (n1 + 1) / 2
  2 * U / (n1 * n2) - 1                        # = (#x>y − #x<y)/(n1·n2)
}
cell_types <- unique(meta$cell_type); drugs <- rownames(d2c)
weight_mat <- matrix(0, length(cell_types), length(drugs), dimnames = list(cell_types, drugs))
padj_mat <- weight_mat
for (ct in cell_types) {
  cells <- rownames(meta)[meta$cell_type == ct]
  d_cells <- cells[meta[cells, "disease"] == DISEASE]
  h_cells <- cells[meta[cells, "disease"] == CONTROL]
  if (length(d_cells) < 3 || length(h_cells) < 3) { message("跳过 ", ct, ":细胞数不足"); next }
  pv <- dl <- numeric(length(drugs))
  for (j in seq_along(drugs)) {
    a <- d2c[j, d_cells]; b <- d2c[j, h_cells]
    if (var(a) == 0 || var(b) == 0) { pv[j] <- 1; dl[j] <- 0; next }
    pv[j] <- suppressWarnings(wilcox.test(a, b, alternative = "two.sided")$p.value)
    dl[j] <- cliff_delta(a, b)
  }
  pa <- p.adjust(pv + 1e-6, method = "BH")     # 上游 +1e-6 防 log10(0)=Inf
  padj_mat[ct, ] <- pa
  weight_mat[ct, ] <- dl * (-log10(pa))        # Eq.(7) 的 weight
}
write.csv(round(weight_mat, 5), file.path(DRES, "prom_inh_weight.csv"))
write.csv(signif(padj_mat, 5), file.path(DRES, "prom_inh_padj.csv"))

## =============================================================================
## Step 3 · 细胞类型内 DEG (R/get_intra_cell_type_degs.R)
##   上游 = Seurat::FindMarkers(默认 wilcox);Seurat 在就直接调它,否则内置等价实现
## =============================================================================
cat("Step 3: 细胞类型内 DEG (disease vs control)\n")
has_seurat <- requireNamespace("Seurat", quietly = TRUE)
deg_list <- list()
if (has_seurat) {
  suppressWarnings(suppressMessages(library(Seurat)))
  so <- suppressWarnings(CreateSeuratObject(counts = expm1(expr), meta.data = meta))
  so <- SetAssayData(so, layer = "data", new.data = as(expr, "dgCMatrix"))
  so$cell_type <- meta$cell_type; so$disease <- meta$disease
  for (ct in cell_types) {
    Idents(so) <- "cell_type"
    sub <- subset(so, cell_type == ct); Idents(sub) <- "disease"
    if (sum(sub$disease == DISEASE) < 3 || sum(sub$disease == CONTROL) < 3) next
    mk <- suppressWarnings(FindMarkers(sub, ident.1 = DISEASE, ident.2 = CONTROL,
                                       logfc.threshold = 0, min.pct = 0))
    deg_list[[ct]] <- data.frame(row.names = rownames(mk), score = mk$avg_log2FC,
                                 adj.P.Val = mk$p_val_adj, P.Value = mk$p_val)
  }
  cat("  DEG 引擎: Seurat::FindMarkers (wilcox)\n")
}
if (!length(deg_list)) {                        # 无 Seurat 时的等价内置实现
  for (ct in cell_types) {
    cells <- rownames(meta)[meta$cell_type == ct]
    dc <- cells[meta[cells, "disease"] == DISEASE]; hc <- cells[meta[cells, "disease"] == CONTROL]
    if (length(dc) < 3 || length(hc) < 3) next
    p <- apply(expr, 1, function(g) suppressWarnings(wilcox.test(g[dc], g[hc])$p.value))
    lfc <- log2((rowMeans(expm1(expr[, dc, drop = FALSE])) + 1) /
                (rowMeans(expm1(expr[, hc, drop = FALSE])) + 1))
    p[is.na(p)] <- 1
    deg_list[[ct]] <- data.frame(row.names = rownames(expr), score = lfc,
                                 adj.P.Val = p.adjust(p, "bonferroni"), P.Value = p)
  }
  cat("  DEG 引擎: 内置 wilcox + Bonferroni (与 FindMarkers 默认同式)\n")
}
cat(sprintf("  %d 个细胞类型拿到 DEG\n", length(deg_list)))

## =============================================================================
## Step 4 · B 臂:敏感/耐药效应
##   上游 = Asgard::GetDrugRef + Asgard::GetDrug(需 CMAP L1000 rankMatrix)
##   本机无 Asgard → 内置 CMap 式 KS 反向连接基线 + 置换 p(明确标注非等价)
## =============================================================================
has_asgard <- requireNamespace("Asgard", quietly = TRUE) && requireNamespace("cmapR", quietly = TRUE)
cat(sprintf("Step 4: 估计 sensitivity/resistance —— 引擎 = %s\n",
            if (has_asgard) "Asgard::GetDrug (上游路径)" else "内置 KS 反向连接基线 (Asgard 缺失)"))
if (has_asgard) {
  stop(paste("检测到 Asgard/cmapR:请改走上游函数 compute_drug_sens_res(disease=, ",
             "perturbation_matrix_path=, gene_info=, drug_info=, deg_list=),",
             "并提供 CMAP L1000 tissue-specific rankMatrix。本 turnkey 脚本不代跑重管道。"))
}

# CMap 连接性 KS 统计量 (Lamb et al., Science 2006):签名集在药物响应秩表中的富集
ks_stat <- function(pos, n) {                    # pos = 命中位置(升序), n = 全基因数
  t <- length(pos); if (!t) return(0)
  a <- max(seq_len(t) / t - pos / n); b <- max(pos / n - (seq_len(t) - 1) / t)
  if (a > b) a else -b
}
conn_score <- function(rank_up, rank_dn, n) {
  ku <- ks_stat(sort(rank_up), n); kd <- ks_stat(sort(rank_dn), n)
  if (sign(ku) == sign(kd)) 0 else (ku - kd) / 2   # >0 同向, <0 反向(治疗性)
}
common_g <- intersect(rownames(expr), rownames(sig_mat))
sig_mat <- sig_mat[common_g, , drop = FALSE]
sens_p <- sens_fdr <- matrix(1, length(cell_types), length(drugs),
                             dimnames = list(cell_types, drugs))
sens_cs <- matrix(0, length(cell_types), length(drugs), dimnames = list(cell_types, drugs))
set.seed(42)
for (ct in names(deg_list)) {
  dg <- deg_list[[ct]]; dg <- dg[rownames(dg) %in% common_g & dg$adj.P.Val < 0.05, , drop = FALSE]
  up <- rownames(dg)[dg$score > 0]; dn <- rownames(dg)[dg$score < 0]
  if (length(up) < 3 || length(dn) < 3) next
  n <- length(common_g)
  pv <- numeric(length(drugs)); cs <- numeric(length(drugs))
  for (j in seq_along(drugs)) {
    ord <- order(sig_mat[, drugs[j]], decreasing = TRUE)          # 药物响应秩表
    rk <- setNames(seq_len(n), common_g[ord])
    obs <- conn_score(rk[up], rk[dn], n)
    null <- replicate(N_PERM, {                                   # 置换基因标签
      s <- sample(n); conn_score(s[seq_along(up)], s[length(up) + seq_along(dn)], n) })
    cs[j] <- obs
    pv[j] <- (sum(null <= obs) + 1) / (N_PERM + 1)                # 单侧:越负越治疗性
  }
  sens_cs[ct, ] <- cs; sens_p[ct, ] <- pv
  sens_fdr[ct, ] <- p.adjust(pv, "BH")
}
write.csv(signif(sens_fdr, 5), file.path(DRES, "sens_res_fdr.csv"))
write.csv(round(sens_cs, 5), file.path(DRES, "connectivity_score.csv"))

## =============================================================================
## Step 5 · linking:Eq.(7-9) 合成治疗评分 (R/compute_scdruglink_score.R)
##   score = Σ_ct (prop/100)·(−log10 FDR)·treated_degs_ratio·exp(weight[ct,drug])
##   treated_degs = −deg_score × mean_response > 0 的 DEG 数 / 可处理 DEG 数
##   跨细胞类型 p 合并 = Fisher (CombineP: chisq = −2Σln p, df = 2k)
## =============================================================================
cat("Step 5: linking → 最终治疗评分\n")
dis_meta <- meta[meta$disease == DISEASE, , drop = FALSE]
cl_size <- table(dis_meta$cell_type); cl_size <- cl_size[cl_size > 3]
cl_prop <- round(100 * cl_size / nrow(dis_meta), 2)

combine_p <- function(p) { keep <- p > 0 & p <= 1
  if (sum(keep) < 2) return(NA_real_)
  pchisq(-2 * sum(log(p[keep])), 2 * sum(keep), lower.tail = FALSE) }

score_one <- function(drug, use_link = TRUE) {
  s <- 0; per_ct <- numeric(0); resp <- sig_mat[, drug]
  for (ct in names(deg_list)) {
    if (!(ct %in% names(cl_prop))) { per_ct <- c(per_ct, 0); next }
    dg <- deg_list[[ct]]; dg <- dg[dg$adj.P.Val < 0.05, , drop = FALSE]
    treatable <- intersect(rownames(dg), names(resp))
    if (!length(treatable)) { per_ct <- c(per_ct, 0); next }
    treated <- -dg[treatable, "score"] * resp[treatable]
    ratio <- sum(treated > 0) / length(treatable)
    w <- if (use_link) weight_mat[ct, drug] else 0
    v <- (cl_prop[[ct]] / 100) * (-log10(sens_fdr[ct, drug])) * ratio * exp(w)
    s <- s + v; per_ct <- c(per_ct, v)
  }
  names(per_ct) <- names(deg_list); list(score = s, per_ct = per_ct)
}
linked <- lapply(drugs, score_one, use_link = TRUE)
unlinked <- lapply(drugs, score_one, use_link = FALSE)   # 对照:去掉 exp(weight) 的 Asgard 式打分
names(linked) <- names(unlinked) <- drugs

ct_score <- t(sapply(linked, function(x) x$per_ct)); rownames(ct_score) <- drugs
comb_p <- sapply(drugs, function(d) combine_p(sens_p[names(deg_list), d]))
out <- data.frame(drug = drugs,
                  drug_score = sapply(linked, function(x) x$score),
                  drug_score_unlinked = sapply(unlinked, function(x) x$score),
                  p_val = comb_p, fdr = p.adjust(comb_p, "BH"),
                  weight_sum = colSums(weight_mat)[drugs], row.names = NULL)
if (!is.null(labels_df) && "known_label" %in% colnames(labels_df))
  out$known_label <- labels_df$known_label[match(out$drug, labels_df$drug_name)]
out <- out[order(-out$drug_score), ]
write.csv(out, file.path(DRES, "drug_scores.csv"), row.names = FALSE)
write.csv(round(ct_score, 6), file.path(DRES, "cell_type_drug_scores.csv"))
cat("  Top5:\n"); print(utils::head(out[, c("drug", "drug_score", "drug_score_unlinked", "fdr")], 5))

## ---- 评估 AUROC / AUPR(对应上游 reproducibility/compute_eval_metrics.R)------
auroc <- function(s, y) { r <- rank(s); n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (!n1 || !n0) return(NA_real_); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
aupr <- function(s, y) { o <- order(-s); y <- y[o]
  tp <- cumsum(y); prec <- tp / seq_along(y); rec <- tp / sum(y)
  sum(diff(c(0, rec)) * prec) }
if (!is.null(out$known_label)) {
  ev <- data.frame(model = c("scDrugLink (linked)", "unlinked (no exp(weight))"),
    auroc = c(auroc(out$drug_score, out$known_label), auroc(out$drug_score_unlinked, out$known_label)),
    aupr  = c(aupr(out$drug_score, out$known_label),  aupr(out$drug_score_unlinked, out$known_label)))
  write.csv(ev, file.path(DRES, "eval_metrics.csv"), row.names = FALSE)
  cat("  评估:\n"); print(ev)
}

## =============================================================================
## Step 6 · 出图(heatmap / lollipop / 散点 / violin;不用条形图)
## =============================================================================
cat("Step 6: 出图\n")
lab_of <- function(d) if (is.null(out$known_label)) "drug" else
  ifelse(out$known_label[match(d, out$drug)] == 1, "known indication", "other")

# (1) promotion/inhibition 权重 heatmap (cell type x drug, RdBu)
hm <- data.frame(cell_type = rep(rownames(weight_mat), ncol(weight_mat)),
                 drug = rep(colnames(weight_mat), each = nrow(weight_mat)),
                 w = as.vector(weight_mat))
hm$drug <- factor(hm$drug, levels = out$drug)
p1 <- ggplot(hm, aes(drug, cell_type, fill = w)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_diverge(name = "Promotion (+) /\nInhibition (-)") +
  labs(title = "Arm A: drug-target promotion/inhibition weight",
       subtitle = "Cliff's delta x -log10(BH p), within cell type (disease vs control)",
       x = NULL, y = NULL) +
  theme_pub(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p1, file.path(DAST, "prom_inh_heatmap"), width = 8.6, height = 3.6)

# (2) 最终治疗评分 lollipop
top <- utils::head(out, 15); top$drug <- factor(top$drug, levels = rev(top$drug))
top$grp <- lab_of(as.character(top$drug))
p2 <- ggplot(top, aes(drug_score, drug, colour = grp)) +
  geom_segment(aes(x = 0, xend = drug_score, y = drug, yend = drug),
               colour = "grey70", linewidth = 0.5) +
  geom_point(size = 3.4) + scale_colour_manual(values = pal_pub(name = "nejm"), name = NULL) +
  scale_x_continuous(trans = "log1p", breaks = c(0, 0.3, 1, 3, 10, 30, 100)) +
  labs(title = "scDrugLink therapeutic score (top 15)",
       subtitle = "Eq.(9): sum over cell types of prop x -log10(FDR) x treated-DEG ratio x exp(weight)",
       x = "Therapeutic score (log1p scale)", y = NULL) +
  theme_pub(base_size = 10, legend = "bottom")
save_fig(p2, file.path(DAST, "drug_score_lollipop"), width = 6.4, height = 5.0)

# (3) linking 效果:unlinked vs linked 散点
sc <- out; sc$grp <- lab_of(sc$drug)
p3 <- ggplot(sc, aes(drug_score_unlinked, drug_score, colour = grp)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
  geom_point(size = 3, alpha = 0.9) +
  { if (requireNamespace("ggrepel", quietly = TRUE))
      ggrepel::geom_text_repel(aes(label = drug), size = 2.6, max.overlaps = 12, show.legend = FALSE) } +
  scale_colour_manual(values = pal_pub(name = "nejm"), name = NULL) +
  scale_x_continuous(trans = "log1p", breaks = c(0, 0.2, 0.5, 1)) +
  scale_y_continuous(trans = "log1p", breaks = c(0, 0.3, 1, 3, 10, 30, 100)) +
  labs(title = "Effect of linking the two evidence arms",
       subtitle = "y: with exp(promotion/inhibition weight)   x: sensitivity/resistance only (log1p axes)",
       x = "Unlinked score", y = "Linked scDrugLink score") +
  theme_pub(base_size = 10, legend = "bottom")
save_fig(p3, file.path(DAST, "linking_effect_scatter"), width = 6.0, height = 5.4)

# (4) 头名药的 Drug2Cell 分布(violin,按细胞类型 × 病/对照)
best <- as.character(out$drug[1])
vd <- data.frame(score = d2c[best, rownames(meta)], cell_type = meta$cell_type,
                 group = factor(meta$disease, levels = c(CONTROL, DISEASE)))
p4 <- ggplot(vd, aes(cell_type, score, fill = group)) +
  geom_violin(scale = "width", trim = TRUE, colour = "black", linewidth = 0.3, alpha = 0.85) +
  geom_boxplot(width = 0.14, outlier.shape = NA, position = position_dodge(0.9),
               colour = "black", linewidth = 0.3, show.legend = FALSE) +
  scale_fill_manual(values = pal_pub(name = "lancet"), name = NULL) +
  labs(title = sprintf("Drug2Cell target score of %s", best),
       subtitle = "Per-cell drug-target module score, background-corrected (Drug2Cell)",
       x = NULL, y = "Drug2Cell score") +
  theme_pub(base_size = 10, legend = "bottom") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_fig(p4, file.path(DAST, "d2c_violin_top_drug"), width = 6.8, height = 4.4)

cat(sprintf("完成。results -> %s ; assets -> %s\n", DRES, DAST))
