# =============================================================================
# 编号       : 060
# 脚本名     : 单基因免疫双蝴蝶相关图 (turnkey + 顶刊图)
# 分类       : 12_tcga_prognosis
# 用途       : 计算目标基因与免疫细胞、免疫检查点基因的相关,绘成两侧发散"蝴蝶"相关图。
# 方法/包    : Spearman 相关 + ggplot2;主题 theme_pub.R
# 结果图     : Immune_butterfly(左:免疫细胞;右:检查点基因)
# -----------------------------------------------------------------------------
# 运行(示例): Rscript 060_immune_butterfly.R
# 运行(自己): Rscript 060_immune_butterfly.R --expr data/expr.csv --immune data/immune.csv --gene TYMS
# 可选参数 : --checkpoints data/checkpoint_genes.txt(检查点基因列表;缺省用内置常见检查点)
# 输入规格 : --expr 表达矩阵(首列基因含目标基因与检查点基因);--immune 免疫比例(首列 Sample);列对齐。
# 整理日期 : 2026-06-23(turnkey 重构;linkET→纯 ggplot 发散蝴蝶图)
# =============================================================================

.find_fw <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  d <- if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
  for (i in 0:6) { cand <- file.path(d, paste(rep("..", i), collapse = "/"), "_framework", "theme_pub.R")
    if (file.exists(cand)) return(normalizePath(cand)) }
  stop("未找到 _framework/theme_pub.R")
}
source(.find_fw())
suppressWarnings(suppressMessages(library(ggplot2)))

SCRIPT_DIR <- bio_script_dir()
args <- bio_args(list(expr = file.path(SCRIPT_DIR, "example_data", "expression.csv"),
                      immune = file.path(SCRIPT_DIR, "example_data", "immune_fraction.csv"),
                      checkpoints = file.path(SCRIPT_DIR, "example_data", "checkpoint_genes.txt"),
                      outdir = file.path(SCRIPT_DIR, "results"), gene = "TYMS"))
ASSETS <- file.path(SCRIPT_DIR, "assets")
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE); dir.create(ASSETS, recursive = TRUE, showWarnings = FALSE)

cat("Step 1/3: 读取数据...\n")
expr <- read_table_smart(args$expr, row_names = TRUE)
imm <- read_table_smart(args$immune, row_names = TRUE)
if (!args$gene %in% rownames(expr)) stop("目标基因不在表达矩阵: ", args$gene)
samples <- intersect(colnames(expr), rownames(imm)); if (length(samples) < 5) stop("表达/免疫样本无法对齐。")
tg <- as.numeric(expr[args$gene, samples])
ckpt <- if (!is.null(args$checkpoints) && file.exists(args$checkpoints)) intersect(trimws(readLines(args$checkpoints)), rownames(expr)) else
  intersect(c("PDCD1", "CD274", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "IDO1"), rownames(expr))
ckpt <- setdiff(ckpt, args$gene)
cat("  目标:", args$gene, "· 免疫细胞", ncol(imm), "· 检查点", length(ckpt), "\n")

cat("Step 2/3: 相关计算...\n")
cor_one <- function(x, panel, name) { ct <- suppressWarnings(cor.test(tg, x, method = "spearman"))
  data.frame(Panel = panel, Var = name, r = unname(ct$estimate), p = ct$p.value) }
res <- rbind(
  do.call(rbind, lapply(colnames(imm), function(c) cor_one(imm[samples, c], "Immune cells", c))),
  do.call(rbind, lapply(ckpt, function(g) cor_one(as.numeric(expr[g, samples]), "Checkpoint genes", g))))
res$sig <- cut(res$p, c(-Inf, .001, .01, .05, Inf), c("***", "**", "*", ""))
res$signed_r <- ifelse(res$Panel == "Immune cells", -res$r, res$r)  # 左侧取负实现蝴蝶发散
write.csv(res[, c("Panel", "Var", "r", "p")], file.path(args$outdir, "correlation.csv"), row.names = FALSE)

cat("Step 3/3: 绘制双蝴蝶图...\n")
res <- res[order(res$Panel, res$r), ]; res$Var <- factor(res$Var, levels = unique(res$Var))
p <- ggplot(res, aes(signed_r, Var, fill = r)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.6) +
  geom_text(aes(label = paste0(sprintf("%.2f", r), sig),
                hjust = ifelse(Panel == "Immune cells", 1.1, -0.1)), size = 2.9) +
  scale_fill_gradient2(low = "#3C5488", mid = "white", high = "#E64B35", midpoint = 0, name = "Spearman r") +
  scale_x_continuous(labels = function(x) abs(x), limits = c(-1, 1), expand = expansion(mult = .15)) +
  facet_grid(Panel ~ ., scales = "free_y", space = "free_y") +
  labs(title = paste0(args$gene, " — immune correlation (butterfly)"),
       subtitle = "← Immune cells          Checkpoint genes →", x = "|Spearman correlation|", y = NULL) +
  theme_pub(base_size = 11, border = TRUE) + theme(axis.text.y = element_text(face = "italic"), strip.text = element_text(face = "bold"))
save_fig(p, file.path(ASSETS, "Immune_butterfly"), 7, 6); save_fig(p, file.path(args$outdir, "Immune_butterfly"), 7, 6)
cat("完成。蝴蝶图见", normalizePath(ASSETS), "\n")
