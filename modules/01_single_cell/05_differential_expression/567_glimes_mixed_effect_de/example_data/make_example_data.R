# 生成 567 模块的合成示例数据 (synthetic, for demo only)
# 设计:12 个供体(6 ctrl / 6 stim),每供体 24 个细胞;150 个基因。
#   - 30 个「真 DE」基因:条件固定效应
#   - 30 个「供体驱动」空基因:无条件效应,但供体随机截距方差大 → 细胞级 t 检验会假阳
#   - 90 个普通空基因
set.seed(2026)

OUT <- local({ a <- commandArgs(FALSE); m <- grep("^--file=", a); if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd() })

n_donor <- 12; n_cell_per <- 24; n_gene <- 150
donors <- sprintf("D%02d", 1:n_donor)
cond   <- rep(c("ctrl", "stim"), each = n_donor / 2)      # 供体层面的条件
batch  <- rep(c("B1", "B2"), times = n_donor / 2)          # 与条件交错,避免完全混杂

meta <- data.frame(
  cell      = sprintf("C%04d", seq_len(n_donor * n_cell_per)),
  donor     = rep(donors, each = n_cell_per),
  condition = rep(cond,   each = n_cell_per),
  exp_batch = rep(batch,  each = n_cell_per),
  stringsAsFactors = FALSE
)
n_cell <- nrow(meta)

gene_id <- sprintf("G%03d", seq_len(n_gene))
kind <- c(rep("true_de", 30), rep("donor_driven_null", 30), rep("plain_null", 90))

# 每个细胞的测序深度因子(log 尺度)
size_factor <- exp(rnorm(n_cell, 0, 0.25))

counts <- matrix(0L, nrow = n_gene, ncol = n_cell,
                 dimnames = list(gene_id, meta$cell))

base_log_mu <- runif(n_gene, log(0.4), log(6))   # 基线表达跨度
for (g in seq_len(n_gene)) {
  # 条件效应
  beta <- if (kind[g] == "true_de") sample(c(-1, 1), 1) * runif(1, 0.7, 1.6) else 0
  # 供体随机截距标准差
  sd_donor <- switch(kind[g],
                     donor_driven_null = runif(1, 0.7, 1.1),
                     true_de           = runif(1, 0.15, 0.35),
                     runif(1, 0.10, 0.25))
  u <- setNames(rnorm(n_donor, 0, sd_donor), donors)
  eta <- base_log_mu[g] + beta * (meta$condition == "stim") + u[meta$donor] + log(size_factor)
  # Poisson-lognormal:再叠一点细胞级过离散
  mu <- exp(eta + rnorm(n_cell, 0, 0.3))
  counts[g, ] <- rpois(n_cell, mu)
}

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

hdr <- c("# synthetic, for demo only -- 由 make_example_data.R 生成,不是真实实验数据")

cdf <- data.frame(gene = gene_id, counts, check.names = FALSE)
writeLines(hdr, file.path(OUT, "counts.csv"))
write.table(cdf, file.path(OUT, "counts.csv"), sep = ",", row.names = FALSE,
            quote = FALSE, append = TRUE, col.names = TRUE)

writeLines(hdr, file.path(OUT, "cell_meta.csv"))
write.table(meta, file.path(OUT, "cell_meta.csv"), sep = ",", row.names = FALSE,
            quote = FALSE, append = TRUE, col.names = TRUE)

truth <- data.frame(gene = gene_id, gene_class = kind, is_de = as.integer(kind == "true_de"))
writeLines(hdr, file.path(OUT, "truth.csv"))
write.table(truth, file.path(OUT, "truth.csv"), sep = ",", row.names = FALSE,
            quote = FALSE, append = TRUE, col.names = TRUE)

cat("done:", n_gene, "genes x", n_cell, "cells\n")
cat("detection rate range:", round(range(rowMeans(counts > 0)), 3), "\n")
