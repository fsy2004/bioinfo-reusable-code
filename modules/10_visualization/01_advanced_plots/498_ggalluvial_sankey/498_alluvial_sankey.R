# =============================================================================
# 编号       : 498
# 脚本名     : ggalluvial 桑基/冲积图 (turnkey + 顶刊图)
# 分类       : 17_advanced_figures
# 用途       : 把多层流向关系(药物→hub→通路 / 配体→受体→细胞 等)绘成冲积/桑基图。
# 方法/包    : ggalluvial;主题 theme_pub.R
# 结果图     : Alluvial(多层冲积图)
# -----------------------------------------------------------------------------
# 运行(示例): Rscript 498_alluvial_sankey.R
# 运行(自己): Rscript 498_alluvial_sankey.R --input data/flow_table.csv
# 输入规格 : CSV 长表,前若干列=各层(2-4 层,如 Drug/Hub/Pathway),最后一列 Freq(流宽,
#            缺省则每行计 1)。
# 整理日期 : 2026-06-23(由代码片段补全为 turnkey 模块)
# =============================================================================

.find_fw <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  d <- if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
  for (i in 0:6) { cand <- file.path(d, paste(rep("..", i), collapse = "/"), "_framework", "theme_pub.R")
    if (file.exists(cand)) return(normalizePath(cand)) }
  stop("未找到 _framework/theme_pub.R")
}
source(.find_fw())
suppressWarnings(suppressMessages({ library(ggplot2); library(ggalluvial) }))

SCRIPT_DIR <- bio_script_dir()
args <- bio_args(list(input = file.path(SCRIPT_DIR, "example_data", "flow_table.csv"),
                      outdir = file.path(SCRIPT_DIR, "results")))
ASSETS <- file.path(SCRIPT_DIR, "assets")
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE); dir.create(ASSETS, recursive = TRUE, showWarnings = FALSE)

cat("Step 1/2: 读取流向表...\n")
df <- read_table_smart(args$input)
freq_col <- if (tolower(names(df)[ncol(df)]) %in% c("freq", "frequency", "weight", "count", "n")) names(df)[ncol(df)] else NULL
axes <- setdiff(names(df), freq_col)
if (length(axes) < 2) stop("需 ≥2 层(列)。")
if (is.null(freq_col)) { df$Freq <- 1; freq_col <- "Freq" }
df <- to_lodes_form(df, axes = axes, id = "flow")
cat("  ", length(axes), "层:", paste(axes, collapse = " → "), "\n")

cat("Step 2/2: 绘制冲积图...\n")
n_fill <- length(unique(df$stratum[df$x == axes[1]]))
p <- ggplot(df, aes(x = x, stratum = stratum, alluvium = flow, y = .data[[freq_col]], fill = stratum)) +
  geom_flow(width = 1/3, alpha = 0.55, curve_type = "sigmoid") +
  geom_stratum(width = 1/3, colour = "grey30", linewidth = 0.3) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3, fontface = "bold") +
  scale_fill_manual(values = pal_pub(length(unique(df$stratum)), "npg"), guide = "none") +
  scale_x_discrete(limits = axes, expand = c(.08, .08)) +
  labs(title = "Alluvial flow", x = NULL, y = "Frequency") +
  theme_pub(base_size = 12) + theme(panel.grid = element_blank(), axis.text.x = element_text(face = "bold", size = 12))
save_fig(p, file.path(ASSETS, "Alluvial"), 8, 6); save_fig(p, file.path(args$outdir, "Alluvial"), 8, 6)
cat("完成。冲积图见", normalizePath(ASSETS), "\n")
