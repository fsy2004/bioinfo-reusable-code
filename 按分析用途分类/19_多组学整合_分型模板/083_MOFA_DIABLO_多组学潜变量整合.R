#!/usr/bin/env Rscript

# MOFA2 / mixOmics DIABLO multi-omics integration wrapper.

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

read_view <- function(spec) {
  parts <- strsplit(spec, "=", fixed = TRUE)[[1]]
  if (length(parts) != 2) stop("View spec must be name=path: ", spec)
  sep <- if (grepl("\\.csv$", parts[2], ignore.case = TRUE)) "," else "\t"
  x <- read.table(parts[2], header = TRUE, row.names = 1, sep = sep, check.names = FALSE)
  list(name = parts[1], data = as.matrix(x))
}

args <- parse_args()
if (is.null(args$views) || is.null(args$outdir) || is.null(args$mode)) {
  stop(paste(
    "Usage:",
    "Rscript 083_MOFA_DIABLO_多组学潜变量整合.R",
    "--mode mofa --views rna=rna.tsv,protein=protein.tsv --outdir results/multiomics",
    "or --mode diablo --metadata meta.tsv --group_col group",
    sep = " "
  ))
}

dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
views <- lapply(strsplit(args$views, ",")[[1]], read_view)
names(views) <- vapply(views, `[[`, character(1), "name")
data_list <- lapply(views, `[[`, "data")

if (args$mode == "mofa") {
  if (!requireNamespace("MOFA2", quietly = TRUE)) stop("Package 'MOFA2' is required for --mode mofa.")
  mofa <- MOFA2::create_mofa(data_list)
  mofa <- MOFA2::prepare_mofa(mofa)
  model <- MOFA2::run_mofa(mofa, outfile = file.path(args$outdir, "mofa_model.hdf5"))
  factors <- MOFA2::get_factors(model, factors = "all", as.data.frame = TRUE)
  weights <- MOFA2::get_weights(model, views = "all", factors = "all", as.data.frame = TRUE)
  write.csv(factors, file.path(args$outdir, "mofa_factors.csv"), row.names = FALSE)
  write.csv(weights, file.path(args$outdir, "mofa_weights.csv"), row.names = FALSE)
  saveRDS(model, file.path(args$outdir, "mofa_model.rds"))
}

if (args$mode == "diablo") {
  if (is.null(args$metadata) || is.null(args$group_col)) stop("--metadata and --group_col are required for DIABLO.")
  if (!requireNamespace("mixOmics", quietly = TRUE)) stop("Package 'mixOmics' is required for --mode diablo.")
  sep <- if (grepl("\\.csv$", args$metadata, ignore.case = TRUE)) "," else "\t"
  meta <- read.table(args$metadata, header = TRUE, row.names = 1, sep = sep, check.names = FALSE)
  common <- Reduce(intersect, c(list(rownames(meta)), lapply(data_list, colnames)))
  if (length(common) < 3) stop("Too few common samples across omics views and metadata.")
  x_block <- lapply(data_list, function(m) t(m[, common, drop = FALSE]))
  y <- factor(meta[common, args$group_col])
  design <- matrix(0.1, ncol = length(x_block), nrow = length(x_block), dimnames = list(names(x_block), names(x_block)))
  diag(design) <- 0
  fit <- mixOmics::block.splsda(X = x_block, Y = y, ncomp = 2, design = design)
  saveRDS(fit, file.path(args$outdir, "diablo_block_splsda.rds"))
  variates <- do.call(cbind, fit$variates)
  write.csv(variates, file.path(args$outdir, "diablo_sample_variates.csv"))
}
