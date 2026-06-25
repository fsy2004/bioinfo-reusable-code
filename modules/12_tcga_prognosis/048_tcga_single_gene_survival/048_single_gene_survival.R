# =============================================================================
# 编号       : 048
# 脚本名     : TCGA 单基因多终点生存曲线 (turnkey + 顶刊图)
# 分类       : 12_tcga_prognosis
# 用途       : 对单个基因按高/低表达分组,绘制 OS/DSS/DFI/PFI 各终点的 KM 生存曲线
#              (含 HR/95%CI/p)。
# 方法/包    : survival + survminer;主题 theme_pub.R
# 结果图     : KM_OS;KM_DSS;KM_DFI;KM_PFI(存在的终点各一张)
# -----------------------------------------------------------------------------
# 运行(示例): Rscript 048_single_gene_survival.R
# 运行(自己): Rscript 048_single_gene_survival.R --input data/gene_survival.csv --gene JPH3
# 可选参数 : --gene 基因列名(默认取除生存列外第一数值列) --cutoff median|optimal
# 输入规格 : CSV,含基因表达列 + 各终点的 <EP>.time 与 <EP>(0/1)成对列(EP∈OS/DSS/DFI/PFI)。
# 整理日期 : 2026-06-23(turnkey 重构;合并表达+生存为单表)
# =============================================================================

.find_fw <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  d <- if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
  for (i in 0:6) { cand <- file.path(d, paste(rep("..", i), collapse = "/"), "_framework", "theme_pub.R")
    if (file.exists(cand)) return(normalizePath(cand)) }
  stop("未找到 _framework/theme_pub.R")
}
source(.find_fw())
suppressWarnings(suppressMessages({ library(survival); library(survminer) }))

SCRIPT_DIR <- bio_script_dir()
args <- bio_args(list(input = file.path(SCRIPT_DIR, "example_data", "gene_survival.csv"),
                      outdir = file.path(SCRIPT_DIR, "results"), gene = "", cutoff = "median"))
ASSETS <- file.path(SCRIPT_DIR, "assets")
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE); dir.create(ASSETS, recursive = TRUE, showWarnings = FALSE)

cat("Step 1/2: 读取数据...\n")
d <- read_table_smart(args$input)
eps <- c("OS", "DSS", "DFI", "PFI"); present <- eps[sapply(eps, function(e) all(c(e, paste0(e, ".time")) %in% names(d)))]
surv_cols <- c(unlist(lapply(eps, function(e) c(e, paste0(e, ".time")))), "sample", "id")
gene <- if (nzchar(args$gene)) args$gene else setdiff(names(d)[sapply(d, is.numeric)], surv_cols)[1]
if (!gene %in% names(d)) stop("未找到基因列: ", gene)
cut <- if (args$cutoff == "median") median(d[[gene]], na.rm = TRUE) else median(d[[gene]], na.rm = TRUE)
d$Group <- factor(ifelse(d[[gene]] > cut, "High", "Low"), levels = c("Low", "High"))
cat("  基因", gene, "· 终点:", paste(present, collapse = "/"), "· High", sum(d$Group == "High"), "/ Low", sum(d$Group == "Low"), "\n")

cat("Step 2/2: 各终点 KM 曲线...\n")
res <- data.frame()
for (e in present) {
  dd <- d[!is.na(d[[paste0(e, ".time")]]), ]; dd$.t <- dd[[paste0(e, ".time")]] / 365; dd$.s <- dd[[e]]
  fit <- survfit(Surv(.t, .s) ~ Group, data = dd)
  cox <- summary(coxph(Surv(.t, .s) ~ Group, data = dd))
  hr <- sprintf("HR = %.2f (%.2f-%.2f)\n%s", cox$conf.int[1], cox$conf.int[3], cox$conf.int[4],
                if (cox$coefficients[1, 5] < 0.001) "p < 0.001" else paste0("p = ", sprintf("%.3f", cox$coefficients[1, 5])))
  res <- rbind(res, data.frame(Endpoint = e, HR = cox$conf.int[1], p = cox$coefficients[1, 5]))
  km <- ggsurvplot(fit, data = dd, conf.int = TRUE, pval = hr, pval.size = 4, risk.table = TRUE,
                   legend.labs = c(paste0(gene, " Low"), paste0(gene, " High")), legend.title = "",
                   xlab = "Time (years)", palette = c("#0072B5", "#BC3C29"), risk.table.height = .26,
                   title = paste0(e, " — ", gene), ggtheme = theme_pub(base_size = 12, border = TRUE))
  for (dest in c(file.path(ASSETS, paste0("KM_", e)), file.path(args$outdir, paste0("KM_", e)))) {
    grDevices::cairo_pdf(paste0(dest, ".pdf"), width = 6, height = 6.2, onefile = FALSE); print(km); dev.off()
    grDevices::png(paste0(dest, ".png"), width = 6, height = 6.2, units = "in", res = 300); print(km); dev.off()
  }
  cat("   ", e, "HR=", sprintf("%.2f", cox$conf.int[1]), "p=", sprintf("%.3f", cox$coefficients[1, 5]), "\n")
}
write.csv(res, file.path(args$outdir, "survival_summary.csv"), row.names = FALSE)
cat("完成。各终点 KM 见", normalizePath(ASSETS), "\n")
