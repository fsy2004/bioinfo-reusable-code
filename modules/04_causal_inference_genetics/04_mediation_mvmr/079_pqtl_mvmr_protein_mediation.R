#!/usr/bin/env Rscript

# pQTL/MVMR and two-step mediation MR wrapper.
# For full MVMR, provide one row per SNP with beta/se for each exposure and outcome.

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

read_table_auto <- function(path) {
  sep <- if (grepl("\\.csv$", path, ignore.case = TRUE)) "," else "\t"
  read.table(path, header = TRUE, sep = sep, stringsAsFactors = FALSE, check.names = FALSE)
}

args <- parse_args()
if (is.null(args$input) || is.null(args$outdir) || is.null(args$exposure_betas) ||
    is.null(args$exposure_ses) || is.null(args$outcome_beta) || is.null(args$outcome_se)) {
  stop(paste(
    "Usage:",
    "Rscript 079_pqtl_mvmr_protein_mediation.R",
    "--input harmonised_mvmr.tsv --outdir results/mvmr",
    "--exposure_betas bx_gene,bx_protein --exposure_ses bxse_gene,bxse_protein",
    "--outcome_beta by --outcome_se byse",
    "[--total_effect 0.20 --direct_effect 0.12]",
    sep = " "
  ))
}

if (!requireNamespace("MVMR", quietly = TRUE)) {
  stop("Package 'MVMR' is required. Install it before running this wrapper.")
}

dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
x <- read_table_auto(args$input)
beta_cols <- strsplit(args$exposure_betas, ",")[[1]]
se_cols <- strsplit(args$exposure_ses, ",")[[1]]
needed <- c(beta_cols, se_cols, args$outcome_beta, args$outcome_se)
missing <- setdiff(needed, colnames(x))
if (length(missing) > 0) stop("Missing columns: ", paste(missing, collapse = ", "))

bx <- as.matrix(x[, beta_cols, drop = FALSE])
bxse <- as.matrix(x[, se_cols, drop = FALSE])
by <- as.numeric(x[[args$outcome_beta]])
byse <- as.numeric(x[[args$outcome_se]])

formatted <- MVMR::format_mvmr(
  BXGs = bx,
  BYG = by,
  seBXGs = bxse,
  seBYG = byse,
  RSID = if ("SNP" %in% colnames(x)) x$SNP else seq_len(nrow(x))
)

strength <- MVMR::strength_mvmr(r_input = formatted, gencov = 0)
pleiotropy <- MVMR::pleiotropy_mvmr(r_input = formatted, gencov = 0)
res <- MVMR::ivw_mvmr(r_input = formatted)

write.csv(res, file.path(args$outdir, "mvmr_ivw_results.csv"), row.names = FALSE)
write.csv(strength, file.path(args$outdir, "mvmr_strength_results.csv"), row.names = FALSE)
write.csv(pleiotropy, file.path(args$outdir, "mvmr_pleiotropy_results.csv"), row.names = FALSE)

if (!is.null(args$total_effect) && !is.null(args$direct_effect)) {
  total <- as.numeric(args$total_effect)
  direct <- as.numeric(args$direct_effect)
  mediation <- data.frame(
    total_effect = total,
    direct_effect = direct,
    indirect_effect = total - direct,
    mediated_fraction = ifelse(total == 0, NA_real_, (total - direct) / total)
  )
  write.csv(mediation, file.path(args$outdir, "two_step_mediation_summary.csv"), row.names = FALSE)
}
