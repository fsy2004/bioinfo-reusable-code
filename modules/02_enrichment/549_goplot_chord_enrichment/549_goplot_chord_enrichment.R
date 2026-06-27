# =============================================================================
# 编号   : 549
# 脚本名 : GOplot 富集环形 / 和弦高级图(代替富集条形图)
# 分类   : 02_enrichment
# 用途   : 把一份「GO/通路富集结果 + 基因差异表达(logFC)」可视化为顶刊级
#          多对多关系图,代替平凡的富集条形图:
#            · GOChord  基因 × 通路「和弦图」——一眼看出哪个基因同时落在哪些通路、
#                       并叠加 logFC 上/下调(条形图永远做不到的多对多关系)
#            · GOCircle 富集「圈图」——外圈散点显示每条通路成员基因的上/下调分布,
#                       内圈 z-score 条 + 旁附通路表(信息密度远高于条形)
#            · GOHeat   基因 × 通路「成员热图」——以颜色编码每个基因参与的通路计数
# ★诚实基线(必须对照,不可只报好看图):
#            内置「富集条形图 vs 升级版 lollipop(棒棒糖图)」对照。条形图是本库
#            绘图铁律明确反对的图型(feedback_avoid_bar_charts):墨多信息少、
#            只能显示一维 -log10(p)、无法表达基因-通路多对多关系。本模块用
#            lollipop 作为「同等信息但更克制」的替代,并进一步用 GOChord/GOHeat
#            展示条形图根本无法表达的二维关系。基线对照图 fig0_baseline_bar_vs_lollipop
#            会同时画出两者,直观证明升级的必要性,而非空口宣称。
# 依赖   : GOplot(自带 EC 示例)· ggplot2 · (framework: theme_pub.R)
# 运行   : Rscript 549_goplot_chord_enrichment.R                       # 零改动跑内置 EC 示例
#          Rscript 549_goplot_chord_enrichment.R --david my_david.csv --genelist my_deg.csv
# 输入   : 见 README ①。三张表(均 GOplot 内置 EC,缺则自动落盘 built-in demo only):
#          (1) david    富集结果:   列 Category,ID,Term,Genes(逗号分隔),adj_pval
#                                    → circle_dat 展开 + 算 z-score + 圈图/lollipop
#          (2) genelist 全转录组 DE: 列 ID(基因名),logFC
#                                    → 给每条通路成员基因上色、算通路激活方向 z-score
#          (3) genes    重点 DEG 子集:列 ID(基因名),logFC(关注的差异基因,通常几十个)
#                                    → GOChord/GOHeat 的「基因×通路」关系矩阵成员
#          说明:GOChord/GOHeat 是「基因-通路多对多关系图」,只在可读数量(~几十)的
#                重点基因上有意义,故需 genes 子集;circle_dat/GOCircle 用全 genelist。
#                这正是 GOplot 官方 vignette 的用法(EC$genelist 全表 + EC$genes 子集)。
# =============================================================================

## ---- 0. 载入框架顶刊主题 + 真包 ---------------------------------------------
.find_fw <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  d <- if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
  for (i in 0:6) { cand <- file.path(d, paste(rep("..", i), collapse = "/"), "_framework", "theme_pub.R")
    if (file.exists(cand)) return(normalizePath(cand)) }
  stop("未找到 _framework/theme_pub.R")
}
source(.find_fw())
suppressWarnings(suppressMessages({ library(GOplot); library(ggplot2); library(grid) }))
set.seed(42)

SCRIPT_DIR <- bio_script_dir()
DDAT   <- file.path(SCRIPT_DIR, "example_data")
ASSETS <- file.path(SCRIPT_DIR, "assets")
args <- bio_args(list(
  david    = file.path(DDAT, "david_enrichment.csv"),  # 富集结果
  genelist = file.path(DDAT, "genelist_logFC.csv"),     # 全转录组差异表达 logFC
  genes    = file.path(DDAT, "genes_subset_logFC.csv"), # 重点 DEG 子集(关系图用)
  outdir   = file.path(SCRIPT_DIR, "results"),
  n_proc   = 7,    # GOChord/GOHeat 取前 N 条显著通路(逗号分隔关系图用,过多则不可读)
  n_circ   = 10))  # GOCircle 圈图展示的通路条数
for (k in c("n_proc","n_circ")) args[[k]] <- as.integer(args[[k]])
for (d in c(DDAT, ASSETS, args$outdir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# 兼容 grob/gtable(GOCircle 返回 gtable 而非 ggplot)的双格式导出。
# 强制白底:GOplot 几个函数默认不设背景,透明背景在多数 PNG 查看器里显示为黑底
# 并吞掉黑色文字/图例 → 这里统一铺白底 + ggsave bg="white"。
save_any <- function(obj, file, width = 7, height = 6, dpi = 300) {
  stem <- sub("\\.(pdf|png)$", "", file)
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  is_gg <- inherits(obj, "ggplot") || inherits(obj, "gg")
  if (is_gg) {
    obj <- obj + theme(plot.background  = element_rect(fill = "white", colour = NA),
                       panel.background = element_rect(fill = "white", colour = NA))
    ggsave(paste0(stem, ".pdf"), obj, width = width, height = height,
           device = grDevices::cairo_pdf, bg = "white")
    ggsave(paste0(stem, ".png"), obj, width = width, height = height, dpi = dpi, bg = "white")
    return(invisible())
  }
  # grob/gtable: 用基础设备绘制,先铺白底
  draw_white <- function() {
    grid::grid.newpage()
    grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))
    grid::grid.draw(obj)
  }
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = width, height = height); draw_white(); grDevices::dev.off()
  grDevices::png(paste0(stem, ".png"), width = width, height = height, units = "in", res = dpi, bg = "white")
    draw_white(); grDevices::dev.off()
  invisible()
}

## ---- 1. 准备输入(无文件 → 落盘 GOplot 内置 EC 示例) -----------------------
cat("Step 1: 准备 david 富集表 + genelist(全) + genes(子集) logFC 表...\n")
if (!(file.exists(args$david) && file.exists(args$genelist) && file.exists(args$genes))) {
  data(EC)  # GOplot 内置:心内皮细胞 GO 富集 + 差异表达(built-in demo only)
  # david: Category/ID/Term/Genes/adj_pval ;Genes 列已是逗号分隔字符串
  write.csv(EC$david, args$david, row.names = FALSE)
  # genelist: 全转录组 logFC(circle_dat 用它给每条通路成员上色 + 算 z-score)
  write.csv(EC$genelist[, c("ID", "logFC")], args$genelist, row.names = FALSE)
  # genes: 重点 DEG 子集(GOChord/GOHeat 关系矩阵成员;官方 vignette 即用 EC$genes)
  write.csv(EC$genes[, c("ID", "logFC")], args$genes, row.names = FALSE)
  cat(sprintf("  [gen] 已写内置 EC 示例 → david %d 通路 · genelist %d 基因 · genes %d DEG (built-in demo only)\n",
              nrow(EC$david), nrow(EC$genelist), nrow(EC$genes)))
}
david    <- read.csv(args$david,    check.names = FALSE, stringsAsFactors = FALSE)
genelist <- read.csv(args$genelist, check.names = FALSE, stringsAsFactors = FALSE)
genes    <- read.csv(args$genes,    check.names = FALSE, stringsAsFactors = FALSE)
cat(sprintf("  读入: %d 条富集通路 · %d 全基因 logFC · %d 重点 DEG\n",
            nrow(david), nrow(genelist), nrow(genes)))

## ---- 2. circle_dat:把富集表展开为「基因-通路」长表 + z-score ---------------
# 真包核心函数:GOplot::circle_dat(terms=david, genes=genelist)
# 返回长表 category/ID/term/count/genes/logFC/adj_pval/zscore
#   zscore = (上调基因数 - 下调基因数)/sqrt(count),反映通路整体激活方向。
cat("Step 2: circle_dat 展开基因-通路长表 + 计算 z-score...\n")
circ <- circle_dat(david, genelist)
write.csv(circ, file.path(args$outdir, "circle_dat_long.csv"), row.names = FALSE)
# 每条通路一行的汇总(画基线条形/lollipop 用)
term_tab <- unique(circ[, c("ID", "term", "adj_pval", "zscore", "count")])
term_tab <- term_tab[order(term_tab$adj_pval), ]
term_tab$neglog10p <- -log10(term_tab$adj_pval)
write.csv(term_tab, file.path(args$outdir, "term_summary.csv"), row.names = FALSE)
cat(sprintf("  展开 %d 行(基因×通路);共 %d 条通路;最显著: %s (p=%.1e)\n",
            nrow(circ), nrow(term_tab), term_tab$term[1], term_tab$adj_pval[1]))

## ---- 3. ★诚实基线对照:富集条形图 vs 升级版 lollipop -------------------------
# 铁律(feedback_avoid_bar_charts):顶刊很少用条形图。这里如实把「条形图」与
# 「lollipop」并排画出 —— 同样的 -log10(p) 信息,但 lollipop 墨水更少、可读性更高,
# 且基线条形图根本无法表达后面 GOChord/GOHeat 的多对多关系。不空喊「条形不好」,
# 用对照图自证。
cat("Step 3: ★诚实基线对照(bar vs lollipop)...\n")
topb <- head(term_tab, min(args$n_circ, nrow(term_tab)))
topb$term <- factor(topb$term, levels = rev(topb$term))   # 由上到下显著性递减

p_bar <- ggplot(topb, aes(term, neglog10p)) +
  geom_col(fill = pal_pub(1, "npg"), width = 0.7) +
  coord_flip() +
  labs(title = "Baseline: enrichment bar chart",
       subtitle = "1-D only; ink-heavy; cannot show gene-pathway links",
       x = NULL, y = expression(-log[10]~adj.~italic(p))) +
  theme_pub(base_size = 11)

p_lol <- ggplot(topb, aes(neglog10p, term)) +
  geom_segment(aes(x = 0, xend = neglog10p, y = term, yend = term),
               colour = "grey70", linewidth = 0.7) +
  geom_point(aes(size = count, colour = zscore)) +
  scale_colour_diverge(midpoint = 0, name = "z-score") +
  scale_size_continuous(range = c(2.5, 7), name = "Gene count") +
  labs(title = "Upgraded: lollipop (preferred)",
       subtitle = "Same significance, less ink; encodes count (size) & direction (color)",
       x = expression(-log[10]~adj.~italic(p)), y = NULL) +
  theme_pub(base_size = 11)

p_base <- compose_panels(list(p_bar, p_lol), ncol = 2, tag = "A")
save_any(p_base, file.path(ASSETS, "fig0_baseline_bar_vs_lollipop"), width = 13, height = 6)
cat("  → fig0_baseline_bar_vs_lollipop(.pdf/.png):条形 vs lollipop 自证对照\n")

## ---- 4. GOChord:基因 × 通路「和弦图」(条形图无法表达的多对多关系)----------
cat("Step 4: GOChord 基因-通路和弦图...\n")
proc_sel <- as.character(term_tab$term[seq_len(min(args$n_proc, nrow(term_tab)))])
# chord_dat: 0/1 关系矩阵(基因×通路)+ 末列 logFC;真包函数。
# genes 用「重点 DEG 子集」(几十个),保证关系图可读;全 genelist 会过密不可读。
chord <- chord_dat(data = circ, genes = genes, process = proc_sel)
# ★清洗:GOChord/GOHeat 在「成员数为 0 的通路列」或「不参与任何选中通路的基因行」上会
#   内部建弧出错(replacement rows mismatch,GOplot 已知 bug)。这里剔空列空行后再画,
#   底层关系数据不变,仅去掉无信息的空行/空列。
.np   <- ncol(chord) - 1L                       # 末列是 logFC
.pcol <- colSums(chord[, seq_len(.np), drop = FALSE]) > 0   # 有成员的通路
chord <- chord[, c(which(.pcol), ncol(chord)), drop = FALSE]
.np   <- ncol(chord) - 1L
.grow <- rowSums(chord[, seq_len(.np), drop = FALSE]) > 0   # 至少落在 1 条通路的基因
chord <- chord[.grow, , drop = FALSE]
if (sum(.pcol) < length(proc_sel) || !all(.grow))
  cat(sprintf("  [clean] 剔除 %d 个空通路列 + %d 个无连接基因行(GOplot 空弧 bug 防护)\n",
              length(proc_sel) - sum(.pcol), sum(!.grow)))
write.csv(chord, file.path(args$outdir, "chord_matrix.csv"))
n_proc_plot <- ncol(chord) - 1L   # 清洗后实际通路数
p_chord <- GOChord(chord, space = 0.02, gene.order = "logFC",
                   gene.space = 0.25, gene.size = 3,
                   lfc.col = c("#B2182B", "#F7F7F7", "#2166AC"),  # 上调红→下调蓝(与 RdBu 一致)
                   ribbon.col = pal_pub(n_proc_plot, "npg"))
save_any(p_chord, file.path(ASSETS, "fig1_GOChord_gene_pathway"), width = 9, height = 9)
cat(sprintf("  → fig1_GOChord(%d 基因 × %d 通路);带 logFC 上下调着色\n",
            nrow(chord), n_proc_plot))

## ---- 5. GOCircle:富集圈图(外散点=成员上/下调,内条=z-score,旁附通路表)----
cat("Step 5: GOCircle 富集圈图...\n")
# GOCircle 返回 gtable(非 ggplot);save_any 已兼容。
p_circle <- GOCircle(circ, nsub = min(args$n_circ, nrow(term_tab)),
                     label.size = 4, rad1 = 2, rad2 = 3, table.legend = TRUE)
save_any(p_circle, file.path(ASSETS, "fig2_GOCircle_enrichment"), width = 11, height = 7)
cat("  → fig2_GOCircle:圈图 + z-score + 通路表\n")

## ---- 6. GOHeat:基因 × 通路「成员热图」-------------------------------------
cat("Step 6: GOHeat 基因-通路成员热图...\n")
# GOHeat 输入需 chord 同构的 0/1 矩阵;nlfc=0 → 颜色编码每基因参与通路计数。
heat_mat <- chord[, -ncol(chord), drop = FALSE]   # 去掉末列 logFC,只留通路关系
p_heat <- GOHeat(heat_mat, nlfc = 0)
save_any(p_heat, file.path(ASSETS, "fig3_GOHeat_membership"), width = 10, height = 6)
cat("  → fig3_GOHeat:基因×通路成员热图(颜色=参与通路计数)\n")

cat("\n完成。结果表见", normalizePath(args$outdir), ";展示图见 assets/\n")
cat("★诚实基线结论:lollipop 以更少墨水承载同等 -log10(p),并额外编码 count(点大小)\n")
cat("  与方向(z-score 颜色);GOChord/GOHeat 进一步表达条形图根本无法表达的基因-通路\n")
cat("  多对多关系 —— 对照图 fig0 已直观自证升级必要性。\n")

sink(file.path(args$outdir, "sessionInfo.txt")); print(sessionInfo()); sink()  # 依赖版本快照(铁律6)
