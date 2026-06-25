#!/usr/bin/env Rscript

# Protein matrix differential analysis template with limma.

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
if (is.null(args$protein_matrix) || is.null(args$metadata) || is.null(args$group_col) || is.null(args$outdir)) {
  stop("Usage: Rscript proteomics_limma_msstats_pipeline.R --protein_matrix protein.tsv --metadata meta.tsv --group_col group --outdir results/proteomics")
}
if (!requireNamespace("limma", quietly = TRUE)) stop("Package 'limma' is required.")

dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
sep1 <- if (grepl("\\.csv$", args$protein_matrix, ignore.case = TRUE)) "," else "\t"
sep2 <- if (grepl("\\.csv$", args$metadata, ignore.case = TRUE)) "," else "\t"
mat <- as.matrix(read.table(args$protein_matrix, header = TRUE, row.names = 1, sep = sep1, check.names = FALSE))
meta <- read.table(args$metadata, header = TRUE, row.names = 1, sep = sep2, check.names = FALSE)
common <- intersect(colnames(mat), rownames(meta))
mat <- mat[, common, drop = FALSE]
meta <- meta[common, , drop = FALSE]
group <- factor(meta[[args$group_col]])
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
fit <- limma::lmFit(mat, design)
contrast <- limma::makeContrasts(contrasts = paste(levels(group)[2], levels(group)[1], sep = "-"), levels = design)
fit2 <- limma::eBayes(limma::contrasts.fit(fit, contrast))
dep <- limma::topTable(fit2, number = Inf, adjust.method = "BH")
write.csv(dep, file.path(args$outdir, "differential_proteins.csv"))
