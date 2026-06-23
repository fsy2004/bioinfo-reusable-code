# =============================================================================
# 编号       : 054
# 脚本名     : WGCNA 共表达网络分析 (turnkey + 顶刊图)
# 分类       : 11_WGCNA共表达网络
# 用途       : 对表达矩阵做 WGCNA:软阈值选择、模块识别、模块-性状相关、hub 基因,
#              输出软阈值图、模块树状图、模块-性状相关热图。
# 方法/包    : WGCNA(pickSoftThreshold/blockwiseModules/moduleEigengenes);主题 theme_pub.R
# 结果图     : SoftThreshold;Module_dendrogram;Module_trait_heatmap
# -----------------------------------------------------------------------------
# 运行(示例): Rscript 054_WGCNA_coexpression.R
# 运行(自己): Rscript 054_WGCNA_coexpression.R --input data/expr.csv --traits data/traits.csv
# 可选参数 : --power 0(0=自动选) --minmodule 30 --mergecut 0.25
# 输入规格 : --input 表达矩阵 CSV(首列基因,样本列);--traits 性状表 CSV(首列 Sample,其余数值性状)。
# 整理日期 : 2026-06-23(turnkey 重构;保留 WGCNA 标准流程)
# =============================================================================

.find_fw <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  d <- if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
  for (i in 0:6) { cand <- file.path(d, paste(rep("..", i), collapse = "/"), "_framework", "theme_pub.R")
    if (file.exists(cand)) return(normalizePath(cand)) }
  stop("未找到 _framework/theme_pub.R")
}
source(.find_fw())
suppressWarnings(suppressMessages({ library(WGCNA); library(ggplot2); library(ComplexHeatmap); library(circlize) }))
options(stringsAsFactors = FALSE)

SCRIPT_DIR <- bio_script_dir()
args <- bio_args(list(input = file.path(SCRIPT_DIR, "example_data", "expr_matrix.csv"),
                      traits = file.path(SCRIPT_DIR, "example_data", "traits.csv"),
                      outdir = file.path(SCRIPT_DIR, "results"), power = "0", minmodule = "30", mergecut = "0.25"))
ASSETS <- file.path(SCRIPT_DIR, "assets")
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE); dir.create(ASSETS, recursive = TRUE, showWarnings = FALSE)

cat("Step 1/4: 读取数据...\n")
expr <- read_table_smart(args$input, row_names = TRUE)
datExpr <- as.data.frame(t(as.matrix(expr)))   # 样本×基因
gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
cat("  ", nrow(datExpr), "样本 x", ncol(datExpr), "基因\n")

# ---- Step 2. 软阈值 ----
cat("Step 2/4: 软阈值选择...\n")
powers <- 1:20
sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "unsigned", verbose = 0)
pw <- if (as.integer(args$power) > 0) as.integer(args$power) else if (!is.na(sft$powerEstimate)) sft$powerEstimate else 6
cat("  选用软阈值 power =", pw, "\n")
sdf <- data.frame(Power = powers, R2 = -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq, MeanK = sft$fitIndices$mean.k.)
p_sft <- ggplot(sdf, aes(Power, R2)) +
  geom_hline(yintercept = 0.85, linetype = "dashed", colour = "#E64B35") +
  geom_text(aes(label = Power), size = 3, colour = "#3C5488") +
  geom_vline(xintercept = pw, linetype = "dotted", colour = "grey50") +
  labs(title = "Scale-free topology fit", x = "Soft threshold (power)", y = expression(Scale-free~R^2)) +
  theme_pub(base_size = 12, border = TRUE)
save_fig(p_sft, file.path(ASSETS, "SoftThreshold"), 6, 5); save_fig(p_sft, file.path(args$outdir, "SoftThreshold"), 6, 5)

# ---- Step 3. 模块识别 + 树状图 ----
cat("Step 3/4: 模块识别...\n")
net <- blockwiseModules(datExpr, power = pw, networkType = "unsigned", TOMType = "unsigned",
                        minModuleSize = as.integer(args$minmodule), mergeCutHeight = as.numeric(args$mergecut),
                        numericLabels = TRUE, saveTOMs = FALSE, verbose = 0, maxBlockSize = ncol(datExpr))
moduleColors <- labels2colors(net$colors)
cat("  识别模块数:", length(unique(moduleColors)), "\n")
write.csv(data.frame(Gene = colnames(datExpr), Module = moduleColors), file.path(args$outdir, "gene_modules.csv"), row.names = FALSE)
for (dest in c(file.path(ASSETS, "Module_dendrogram"), file.path(args$outdir, "Module_dendrogram"))) {
  grDevices::cairo_pdf(paste0(dest, ".pdf"), width = 9, height = 5.5)
  plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]], "Module", dendroLabels = FALSE,
                      hang = 0.03, addGuide = TRUE, main = "Gene clustering dendrogram & modules"); dev.off()
  grDevices::png(paste0(dest, ".png"), width = 9, height = 5.5, units = "in", res = 300)
  plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]], "Module", dendroLabels = FALSE,
                      hang = 0.03, addGuide = TRUE, main = "Gene clustering dendrogram & modules"); dev.off()
}

# ---- Step 4. 模块-性状相关热图 ----
cat("Step 4/4: 模块-性状相关...\n")
tr <- read_table_smart(args$traits); rownames(tr) <- tr[[1]]; tr <- tr[rownames(datExpr), -1, drop = FALSE]
tr <- as.data.frame(lapply(tr, as.numeric)); rownames(tr) <- rownames(datExpr)
MEs <- orderMEs(moduleEigengenes(datExpr, moduleColors)$eigengenes)
mtcor <- cor(MEs, tr, use = "p"); mtp <- corPvalueStudent(mtcor, nrow(datExpr))
txt <- matrix(paste0(sprintf("%.2f", mtcor), "\n(", sprintf("%.0e", mtp), ")"), nrow = nrow(mtcor))
ht <- Heatmap(mtcor, name = "correlation", col = colorRamp2(c(-1, 0, 1), c("#3C5488", "white", "#E64B35")),
              rect_gp = gpar(col = "grey80"), cluster_rows = FALSE, cluster_columns = FALSE,
              column_title = "Module-trait relationships", column_title_gp = gpar(fontsize = 12, fontface = "bold"),
              row_names_gp = gpar(fontsize = 9), column_names_gp = gpar(fontsize = 10),
              cell_fun = function(j, i, x, y, w, h, fill) grid.text(txt[i, j], x, y, gp = gpar(fontsize = 7)))
for (dest in c(file.path(ASSETS, "Module_trait_heatmap"), file.path(args$outdir, "Module_trait_heatmap"))) {
  grDevices::cairo_pdf(paste0(dest, ".pdf"), width = max(4, ncol(tr) * 1.5 + 2), height = max(4, nrow(mtcor) * 0.55 + 1.5)); draw(ht); dev.off()
  grDevices::png(paste0(dest, ".png"), width = max(4, ncol(tr) * 1.5 + 2), height = max(4, nrow(mtcor) * 0.55 + 1.5), units = "in", res = 300); draw(ht); dev.off()
}
cat("完成。WGCNA 图/表见", normalizePath(args$outdir), "\n")
