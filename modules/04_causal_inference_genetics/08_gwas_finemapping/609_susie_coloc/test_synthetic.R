# API/guard validation using synthetic summary statistics. No scientific claims.
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
module_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
test_dir <- tempfile("susie_coloc_synthetic_")
dir.create(test_dir)
nvar <- 30L
R <- outer(seq_len(nvar), seq_len(nvar), function(i, j) 0.7^abs(i-j))
ids <- paste0("synthetic_", seq_len(nvar))
ldv <- data.frame(variant_id = ids, ea = "A", nea = "C")
write.table(R, file.path(test_dir, "ld.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(ldv, file.path(test_dir, "ld_variants.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write_summary <- function(path, causal) {
  d <- data.frame(variant_id = ids, ea = "A", nea = "C", beta = as.numeric(R[,causal] * 12 * 0.01), se = 0.01, eaf = 0.2)
  write.table(d, path, sep = "\t", row.names = FALSE, quote = FALSE)
}
write_summary(file.path(test_dir, "trait1.tsv"), 8)
write_summary(file.path(test_dir, "trait2.tsv"), 8)
trait <- function(id) list(trait_id = id, phenotype = "synthetic quantitative trait", genome_build = "GRCh38", ancestry = "SYNTHETIC", beta_unit = "SD", type = "quant", N = 50000, sdY = 1, source = "synthetic contract test", region = "synthetic_region", dense_region_confirmed = TRUE, summary = file.path(test_dir, paste0(id, ".tsv")), ld_matrix = file.path(test_dir, "ld.txt"), ld_variants = file.path(test_dir, "ld_variants.tsv"), ld_genome_build = "GRCh38", ld_ancestry = "SYNTHETIC", ld_metric = "signed_r", ld_coded_allele = "ea", ld_source = "synthetic AR1 matrix")
cfg <- list(seed = 1L, out_dir = file.path(test_dir, "results"), L = 2L, coverage = 0.95, max_iter = 1000L, allow_case_control_rss_approximation = FALSE, prior_grid = data.frame(p1 = c(1e-4,1e-4), p2 = c(1e-4,1e-4), p12 = c(1e-6,1e-5)), trait1 = trait("trait1"), trait2 = trait("trait2"))
config_path <- file.path(test_dir, "config.R")
save_cfg <- function() {
  con <- file(config_path, open = "wt")
  on.exit(close(con))
  writeLines("cfg <-", con)
  dput(cfg, file = con)
}
save_cfg()
rscript <- file.path(R.home("bin"), "Rscript")
run <- function() system2(rscript, c(shQuote(file.path(module_dir, "run_susie_coloc.R")), shQuote(config_path)), stdout = TRUE, stderr = TRUE)
logs <- run()
stopifnot(is.null(attr(logs, "status")))
result <- read.delim(file.path(cfg$out_dir, "susie_coloc_h0_h4_prior_sensitivity.tsv"), check.names = FALSE)
stopifnot(all(c("PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf") %in% names(result)))
stopifnot(all(result$PP.H4.abf > 0.8), length(unique(result$prior_id)) == 2)
# A second run must reject existing results without overwriting any output.
before_md5 <- tools::md5sum(file.path(cfg$out_dir, "variant_pip.tsv"))
logs_existing <- suppressWarnings(run())
stopifnot(!is.null(attr(logs_existing, "status")), attr(logs_existing, "status") != 0)
stopifnot(identical(before_md5, tools::md5sum(file.path(cfg$out_dir, "variant_pip.tsv"))))
# Distinct well-separated causal signals should favor H3 over H4.
write_summary(file.path(test_dir, "trait2.tsv"), 24)
cfg$out_dir <- file.path(test_dir, "distinct_results"); save_cfg()
logs_distinct <- run()
stopifnot(is.null(attr(logs_distinct, "status")))
distinct <- read.delim(file.path(cfg$out_dir, "susie_coloc_h0_h4_prior_sensitivity.tsv"), check.names = FALSE)
stopifnot(all(distinct$PP.H3.abf > distinct$PP.H4.abf))
# Flipped allele coding must stop before fitting; no automatic SNP flipping.
bad <- read.delim(cfg$trait2$summary); bad$ea[[1]] <- "C"; bad$nea[[1]] <- "A"
write.table(bad, cfg$trait2$summary, sep = "\t", row.names = FALSE, quote = FALSE)
cfg$out_dir <- file.path(test_dir, "invalid_results"); save_cfg()
logs_bad <- suppressWarnings(run())
stopifnot(!is.null(attr(logs_bad, "status")), attr(logs_bad, "status") != 0)
cat("PASS synthetic shared/distinct signal API, allele guard and existing-output rejection; not scientific analysis\n")
cat("Temporary synthetic files:", test_dir, "\n")
