# Verify that local PLINK-based LD clumping works without OpenGWAS.
#
# Usage:
#   Rscript verify_local_clumping.R --bfile C:/path/to/1000G/EUR [--n 1000]
# or set MR_LD_BFILE to the PLINK prefix (without .bed/.bim/.fam).

suppressPackageStartupMessages({
  library(dplyr)
  library(ieugwasr)
})

argv <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default = "") {
  key <- paste0("--", name)
  idx <- match(key, argv)
  if (!is.na(idx) && idx < length(argv)) argv[[idx + 1L]] else default
}

bfile <- arg_value("bfile", Sys.getenv("MR_LD_BFILE", unset = ""))
max_variants <- as.integer(arg_value("n", "1000"))
if (is.na(max_variants) || max_variants < 2L) stop("--n must be an integer >= 2.")
if (!nzchar(bfile)) {
  stop("Provide --bfile <PLINK-prefix> or set MR_LD_BFILE.")
}
bfile <- normalizePath(bfile, winslash = "/", mustWork = FALSE)
required <- paste0(bfile, c(".bed", ".bim", ".fam"))
missing <- required[!file.exists(required)]
if (length(missing)) {
  stop("Missing PLINK reference files: ", paste(missing, collapse = ", "))
}

plink <- if (requireNamespace("plinkbinr", quietly = TRUE)) {
  plinkbinr::get_plink_exe()
} else if (requireNamespace("genetics.binaRies", quietly = TRUE)) {
  genetics.binaRies::get_plink_binary()
} else {
  Sys.getenv("MR_PLINK_BIN", unset = "")
}
if (!nzchar(plink) || !file.exists(plink)) {
  stop("No local PLINK 1.9 binary found; install plinkbinr/genetics.binaRies or set MR_PLINK_BIN.")
}

bim <- read.table(paste0(bfile, ".bim"), header = FALSE, stringsAsFactors = FALSE)
chr1 <- bim[bim[[1]] == 1, , drop = FALSE]
if (!nrow(chr1)) stop("The supplied panel contains no chromosome 1 variants.")

set.seed(1)
rsids <- chr1[seq_len(min(max_variants, nrow(chr1))), 2]
dat <- tibble(rsid = rsids, pval = runif(length(rsids), 1e-30, 1e-8), id = "local-clump-check")
clumped <- ld_clump(
  dat,
  plink_bin = plink,
  bfile = bfile,
  clump_r2 = 0.001,
  clump_kb = 10000
)

cat(sprintf("Panel variants: %d\n", nrow(bim)))
cat(sprintf(
  "Local clumping: %d chromosome-1 variants -> %d independent variants (%.1f%% removed)\n",
  nrow(dat), nrow(clumped), 100 * (1 - nrow(clumped) / nrow(dat))
))
cat("PASS: local PLINK clumping completed without an OpenGWAS API call.\n")
