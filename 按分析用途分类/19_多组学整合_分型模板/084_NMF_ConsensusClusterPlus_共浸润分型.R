#!/usr/bin/env Rscript

# NMF and ConsensusClusterPlus wrapper for expression, immune-score or niche-pattern matrices.

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
if (is.null(args$matrix) || is.null(args$outdir)) {
  stop("Usage: Rscript 084_NMF_ConsensusClusterPlus_共浸润分型.R --matrix feature_by_sample.tsv --outdir results/pattern --k_min 2 --k_max 6")
}

dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
k_min <- if (!is.null(args$k_min)) as.integer(args$k_min) else 2
k_max <- if (!is.null(args$k_max)) as.integer(args$k_max) else 6
sep <- if (grepl("\\.csv$", args$matrix, ignore.case = TRUE)) "," else "\t"
x <- as.matrix(read.table(args$matrix, header = TRUE, row.names = 1, sep = sep, check.names = FALSE))
x <- x[apply(x, 1, var, na.rm = TRUE) > 0, , drop = FALSE]
x <- x - min(x, na.rm = TRUE)

if (requireNamespace("NMF", quietly = TRUE)) {
  ranks <- k_min:k_max
  nmf_res <- NMF::nmf(x, rank = ranks, nrun = 30, seed = 1)
  saveRDS(nmf_res, file.path(args$outdir, "nmf_rank_survey.rds"))
  best_rank <- ranks[which.max(NMF::summary(nmf_res)$cophenetic)]
  best <- NMF::nmf(x, rank = best_rank, nrun = 50, seed = 1)
  write.csv(NMF::basis(best), file.path(args$outdir, "nmf_basis_features.csv"))
  write.csv(NMF::coef(best), file.path(args$outdir, "nmf_sample_coefficients.csv"))
  write.csv(data.frame(sample = colnames(x), subtype = NMF::predict(best)),
            file.path(args$outdir, "nmf_sample_subtypes.csv"), row.names = FALSE)
}

if (requireNamespace("ConsensusClusterPlus", quietly = TRUE)) {
  ConsensusClusterPlus::ConsensusClusterPlus(
    x,
    maxK = k_max,
    reps = 100,
    pItem = 0.8,
    pFeature = 0.8,
    title = file.path(args$outdir, "consensus"),
    clusterAlg = "hc",
    distance = "pearson",
    seed = 1,
    plot = "pdf"
  )
}

write.csv(data.frame(features = nrow(x), samples = ncol(x), k_min = k_min, k_max = k_max),
          file.path(args$outdir, "pattern_run_summary.csv"), row.names = FALSE)
