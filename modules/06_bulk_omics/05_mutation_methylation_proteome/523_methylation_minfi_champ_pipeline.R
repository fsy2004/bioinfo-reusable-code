#!/usr/bin/env Rscript

# Methylation differential analysis template for beta-value matrices.

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list()
  i <- 1
  while (i <= length(args)) {
    key <- sub("^--", "", args[[i]])
    val <- if (i + 1 <= length(args) && !grepl("^--", args[[i + 1]])) args[[i + 1]] else TRUE
    out[[key]] <- val
    i <- i + if (isTRUE(val)) 1 else 2
  }
  out
}

args <- parse_args()
if (is.null(args$beta_matrix) || is.null(args$metadata) || is.null(args$group_col) || is.null(args$outdir)) {
  stop("Usage: Rscript methylation_minfi_champ_pipeline.R --beta_matrix beta.tsv --metadata meta.tsv --group_col group --outdir results/methylation")
}
if (!requireNamespace("limma", quietly = TRUE)) stop("Package 'limma' is required.")

dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
sep1 <- if (grepl("\\.csv$", args$beta_matrix, ignore.case = TRUE)) "," else "\t"
sep2 <- if (grepl("\\.csv$", args$metadata, ignore.case = TRUE)) "," else "\t"
beta <- as.matrix(read.table(args$beta_matrix, header = TRUE, row.names = 1, sep = sep1, check.names = FALSE))
meta <- read.table(args$metadata, header = TRUE, row.names = 1, sep = sep2, check.names = FALSE)
common <- intersect(colnames(beta), rownames(meta))
beta <- beta[, common, drop = FALSE]
meta <- meta[common, , drop = FALSE]
group <- factor(meta[[args$group_col]])
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
fit <- limma::lmFit(beta, design)
contrast <- limma::makeContrasts(contrasts = paste(levels(group)[2], levels(group)[1], sep = "-"), levels = design)
fit2 <- limma::eBayes(limma::contrasts.fit(fit, contrast))
dmp <- limma::topTable(fit2, number = Inf, adjust.method = "BH")
write.csv(dmp, file.path(args$outdir, "differential_methylation_probes.csv"))
