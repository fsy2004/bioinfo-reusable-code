# TwoSampleMR + coloc target evidence-chain wrapper
# Example:
# Rscript 075_twosamplemr_coloc_drug_target.R exposure=exposure.csv outcome=outcome.csv locus=locus.csv output_dir=results/mr_coloc

args <- commandArgs(trailingOnly = TRUE)
kv <- strsplit(args, "=", fixed = TRUE)
opts <- setNames(vapply(kv, function(x) if (length(x) > 1) x[2] else "", ""), vapply(kv, `[`, "", 1))

get_opt <- function(name, default = NULL) {
  value <- opts[[name]]
  if (is.null(value) || identical(value, "")) default else value
}

exposure_file <- get_opt("exposure")
outcome_file <- get_opt("outcome")
locus_file <- get_opt("locus")
output_dir <- get_opt("output_dir", "results/mr_coloc")

if (is.null(exposure_file) || is.null(outcome_file)) {
  stop("Please provide exposure=exposure.csv and outcome=outcome.csv")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required <- c("TwoSampleMR", "coloc", "dplyr", "readr")
missing_pkgs <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "))
}

library(TwoSampleMR)
library(coloc)
library(dplyr)
library(readr)

exp_raw <- read_csv(exposure_file, show_col_types = FALSE)
out_raw <- read_csv(outcome_file, show_col_types = FALSE)

exp_dat <- format_data(
  exp_raw,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval"
)

out_dat <- format_data(
  out_raw,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval"
)

harmonised <- harmonise_data(exp_dat, out_dat)
mr_res <- mr(harmonised)
heterogeneity <- mr_heterogeneity(harmonised)
pleiotropy <- mr_pleiotropy_test(harmonised)

write_csv(harmonised, file.path(output_dir, "harmonised_mr_input.csv"))
write_csv(mr_res, file.path(output_dir, "mr_results.csv"))
write_csv(heterogeneity, file.path(output_dir, "mr_heterogeneity.csv"))
write_csv(pleiotropy, file.path(output_dir, "mr_pleiotropy.csv"))

if (!is.null(locus_file)) {
  locus <- read_csv(locus_file, show_col_types = FALSE)
  needed <- c("SNP", "beta_exposure", "varbeta_exposure", "beta_outcome", "varbeta_outcome", "MAF")
  if (!all(needed %in% colnames(locus))) {
    stop("locus file must contain: ", paste(needed, collapse = ", "))
  }
  coloc_res <- coloc.abf(
    dataset1 = list(
      snp = locus$SNP,
      beta = locus$beta_exposure,
      varbeta = locus$varbeta_exposure,
      MAF = locus$MAF,
      type = "quant"
    ),
    dataset2 = list(
      snp = locus$SNP,
      beta = locus$beta_outcome,
      varbeta = locus$varbeta_outcome,
      MAF = locus$MAF,
      type = "quant"
    )
  )
  write_csv(as.data.frame(t(coloc_res$summary)), file.path(output_dir, "coloc_summary.csv"))
  write_csv(coloc_res$results, file.path(output_dir, "coloc_snp_results.csv"))
}

message("Done. MR/coloc evidence-chain results saved to: ", normalizePath(output_dir))

