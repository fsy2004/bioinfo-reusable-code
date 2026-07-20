#!/usr/bin/env Rscript

# Metabolomics differential matrix template.

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
if (is.null(args$metabolite_matrix) || is.null(args$metadata) || is.null(args$group_col) || is.null(args$outdir)) {
  stop("Usage: Rscript metabolomics_metaboanalystR_pipeline.R --metabolite_matrix metabolite.tsv --metadata meta.tsv --group_col group --outdir results/metabolomics")
}

dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
sep1 <- if (grepl("\\.csv$", args$metabolite_matrix, ignore.case = TRUE)) "," else "\t"
sep2 <- if (grepl("\\.csv$", args$metadata, ignore.case = TRUE)) "," else "\t"
mat <- as.matrix(read.table(args$metabolite_matrix, header = TRUE, row.names = 1, sep = sep1, check.names = FALSE))
meta <- read.table(args$metadata, header = TRUE, row.names = 1, sep = sep2, check.names = FALSE)
common <- intersect(colnames(mat), rownames(meta))
mat <- mat[, common, drop = FALSE]
meta <- meta[common, , drop = FALSE]
group <- factor(meta[[args$group_col]])
if (nlevels(group) != 2) stop("This minimal template currently expects exactly two groups.")

pvals <- apply(mat, 1, function(v) t.test(v[group == levels(group)[2]], v[group == levels(group)[1]])$p.value)
logfc <- rowMeans(mat[, group == levels(group)[2], drop = FALSE], na.rm = TRUE) -
  rowMeans(mat[, group == levels(group)[1], drop = FALSE], na.rm = TRUE)
res <- data.frame(metabolite = rownames(mat), logFC = logfc, pvalue = pvals, padj = p.adjust(pvals, "BH"))
res <- res[order(res$padj), ]
write.csv(res, file.path(args$outdir, "differential_metabolites.csv"), row.names = FALSE)
