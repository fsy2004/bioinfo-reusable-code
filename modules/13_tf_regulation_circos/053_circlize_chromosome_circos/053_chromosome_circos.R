# =============================================================================
# 编号       : 053
# 脚本名     : circlize 基因染色体圈图 (turnkey + 顶刊图)
# 分类       : 13_tf_regulation_circos
# 用途       : 把一组目标基因按基因组坐标画到染色体圈图(ideogram)上,带基因标签。
# 方法/包    : circlize(circos.initializeWithIdeogram / genomicLabels)
# 结果图     : Chromosome_circos(染色体圈图 + 基因标签)
# -----------------------------------------------------------------------------
# 运行(示例): Rscript 053_chromosome_circos.R
# 运行(自己): Rscript 053_chromosome_circos.R --input data/gene_positions.csv --genome hg38
# 输入规格 : CSV,列 Gene, Chr(如 chr17), Start, End(基因组坐标)。
# 备注     : 首次运行会联网获取该基因组 cytoband(circlize 缓存);无网络时改 --genome hg19(内置)。
# 整理日期 : 2026-06-23(turnkey 重构;合并双输入为单坐标表)
# =============================================================================

.find_fw <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  d <- if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
  for (i in 0:6) { cand <- file.path(d, paste(rep("..", i), collapse = "/"), "_framework", "theme_pub.R")
    if (file.exists(cand)) return(normalizePath(cand)) }
  stop("未找到 _framework/theme_pub.R")
}
source(.find_fw())
suppressWarnings(suppressMessages(library(circlize)))

SCRIPT_DIR <- bio_script_dir()
args <- bio_args(list(input = file.path(SCRIPT_DIR, "example_data", "gene_positions.csv"),
                      outdir = file.path(SCRIPT_DIR, "results"), genome = "hg38"))
ASSETS <- file.path(SCRIPT_DIR, "assets")
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE); dir.create(ASSETS, recursive = TRUE, showWarnings = FALSE)

cat("Step 1/2: 读取基因坐标...\n")
gp <- read_table_smart(args$input)
names(gp)[1:4] <- c("genename", "chr", "start", "end")
bed <- gp[, c("chr", "start", "end", "genename")]
bed$start <- as.numeric(bed$start); bed$end <- as.numeric(bed$end)
cat("  ", nrow(bed), "个基因\n")

cat("Step 2/2: 绘制染色体圈图...\n")
cyto <- tryCatch(read.cytoband(species = args$genome)$df, error = function(e) NULL)
pal <- colorRampPalette(pal_pub(name = "npg"))(24)
draw_circos <- function() {
  circos.clear()
  circos.par(start.degree = 90, gap.degree = 1)
  if (!is.null(cyto)) circos.initializeWithIdeogram(cytoband = cyto, plotType = NULL)
  else circos.initializeWithIdeogram(species = args$genome, plotType = NULL)
  circos.track(ylim = c(0, 1), track.height = 0.12, bg.border = NA, panel.fun = function(x, y) {
    chr <- CELL_META$sector.index; xlim <- CELL_META$xlim
    idx <- suppressWarnings(as.numeric(gsub("chr", "", chr)))
    idx <- if (is.na(idx)) ifelse(chr == "chrX", 23, ifelse(chr == "chrY", 24, 1)) else idx
    circos.rect(xlim[1], 0, xlim[2], 1, col = pal[min(idx, 24)], border = NA)
    circos.text(mean(xlim), 0.5, gsub("chr", "", chr), cex = 0.55, col = "white", facing = "inside", niceFacing = TRUE)
  })
  circos.genomicIdeogram(track.height = mm_h(4))
  circos.genomicLabels(bed, labels.column = 4, side = "inside", cex = 0.75,
                       col = "black", line_col = "grey50")
  circos.clear()
}
for (dest in c(file.path(ASSETS, "Chromosome_circos"), file.path(args$outdir, "Chromosome_circos"))) {
  grDevices::cairo_pdf(paste0(dest, ".pdf"), width = 7, height = 7); draw_circos(); dev.off()
  grDevices::png(paste0(dest, ".png"), width = 7, height = 7, units = "in", res = 300); draw_circos(); dev.off()
}
write.csv(gp[, c("genename", "chr", "start", "end")], file.path(args$outdir, "gene_position_info.csv"), row.names = FALSE)
cat("完成。圈图见", normalizePath(ASSETS), "\n")
