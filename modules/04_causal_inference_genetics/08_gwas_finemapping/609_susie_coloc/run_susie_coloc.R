#!/usr/bin/env Rscript
# Dense regional fine-mapping + signal-pair colocalization. No downloads/installations.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript run_susie_coloc.R config.R")
for (pkg in c("susieR", "coloc")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Required package unavailable: ", pkg)
}
config_env <- new.env(parent = baseenv())
sys.source(args[[1]], envir = config_env)
cfg <- config_env$cfg
if (!is.list(cfg)) stop("config.R must define list cfg")
if (is.null(cfg$seed)) stop("Declare a random seed")
if (!all(c("out_dir", "L", "coverage", "max_iter", "trait1", "trait2", "prior_grid") %in% names(cfg))) stop("Missing configuration fields")
if (length(cfg$L) != 1 || !is.finite(cfg$L) || cfg$L < 1 || cfg$L != as.integer(cfg$L)) stop("L must be a positive integer")
if (length(cfg$coverage) != 1 || !is.finite(cfg$coverage) || cfg$coverage <= 0 || cfg$coverage >= 1) stop("coverage must be within (0,1)")
if (length(cfg$max_iter) != 1 || !is.finite(cfg$max_iter) || cfg$max_iter < 1) stop("max_iter must be positive")
set.seed(cfg$seed)
if (length(cfg$out_dir) != 1 || !is.character(cfg$out_dir) || !nzchar(cfg$out_dir)) stop("Declare one output directory")
if (file.exists(cfg$out_dir) && !dir.exists(cfg$out_dir)) stop("Output path exists and is not a directory")
if (dir.exists(cfg$out_dir) && length(list.files(cfg$out_dir, all.files = TRUE, no.. = TRUE)) > 0) stop("Output directory is non-empty; choose a new directory to preserve existing results")
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
assert <- function(condition, message) if (!isTRUE(condition)) stop(message, call. = FALSE)
required_meta <- c("trait_id", "phenotype", "genome_build", "ancestry", "beta_unit", "type", "N", "source", "region")
check_trait <- function(t) {
  assert(all(required_meta %in% names(t)), "Missing required trait metadata")
  assert(t$genome_build %in% c("GRCh37", "GRCh38"), "Declare GRCh37/GRCh38")
  assert(isTRUE(t$dense_region_confirmed), "Provide dense regional summary data, not significant/clumped variants")
  assert(t$type %in% c("quant", "cc"), "Trait type must be quant or cc")
  assert(length(t$N) == 1 && is.finite(t$N) && t$N > 1, "N must be one finite sample size")
  if (t$type == "quant") assert(length(t$sdY) == 1 && is.finite(t$sdY) && t$sdY > 0, "Quantitative trait requires documented sdY")
  if (t$type == "cc") {
    assert(length(t$s) == 1 && is.finite(t$s) && t$s > 0 && t$s < 1, "Case-control trait requires sample case fraction s")
    assert(isTRUE(cfg$allow_case_control_rss_approximation), "Binary trait RSS is an approximation; explicitly acknowledge in config")
  }
  d <- read_tsv(t$summary)
  assert(all(c("variant_id", "ea", "nea", "beta", "se", "eaf") %in% names(d)), "Summary requires variant_id, ea, nea, beta, se, eaf")
  assert(nrow(d) >= 3 && !anyDuplicated(d$variant_id) && !anyNA(d$variant_id), "Need at least 3 unique variant identifiers")
  assert(all(nzchar(d$variant_id)), "Empty variant identifier")
  assert(all(is.finite(d$beta)) && all(is.finite(d$se)) && all(d$se > 0), "Invalid beta/se")
  assert(all(is.finite(d$eaf)) && all(d$eaf > 0 & d$eaf < 1), "Invalid effect allele frequency")
  d$ea <- toupper(d$ea); d$nea <- toupper(d$nea)
  assert(all(grepl("^[ACGT]+$", d$ea)) && all(grepl("^[ACGT]+$", d$nea)) && all(d$ea != d$nea), "Normalize explicit alleles before use")
  if ("n" %in% names(d)) assert(all(is.finite(d$n)) && all(abs(d$n - t$N) < 1e-6), "Per-SNP n differs from N; this entry requires constant-N input")
  ldv <- read_tsv(t$ld_variants)
  assert(all(c("variant_id", "ea", "nea") %in% names(ldv)), "LD variants require variant_id, ea, nea")
  assert(identical(as.character(d$variant_id), as.character(ldv$variant_id)), "Exact summary/LD SNP order is required; no silent intersections")
  assert(identical(d$ea, toupper(ldv$ea)) && identical(d$nea, toupper(ldv$nea)), "LD must be coded to the same effect allele")
  assert(identical(t$ld_genome_build, t$genome_build) && identical(t$ld_ancestry, t$ancestry), "LD build/ancestry must match summary metadata")
  assert(identical(t$ld_metric, "signed_r") && identical(t$ld_coded_allele, "ea") && nzchar(t$ld_source), "Require signed r coded to ea and documented LD source")
  R <- as.matrix(read.table(t$ld_matrix, header = FALSE, check.names = FALSE))
  storage.mode(R) <- "double"
  assert(identical(dim(R), c(nrow(d), nrow(d))) && all(is.finite(R)), "LD dimension/non-finite values")
  assert(max(abs(R - base::t(R))) < 1e-6 && max(abs(diag(R) - 1)) < 1e-6 && max(abs(R)) <= 1 + 1e-6, "LD must be a symmetric correlation matrix")
  assert(min(eigen(R, symmetric = TRUE, only.values = TRUE)$values) >= -1e-6, "LD not positive semidefinite")
  dimnames(R) <- list(d$variant_id, d$variant_id)
  list(table = d, R = R, meta = t)
}
t1 <- check_trait(cfg$trait1)
t2 <- check_trait(cfg$trait2)
assert(identical(t1$meta$genome_build, t2$meta$genome_build), "Two traits use different builds")
assert(identical(t1$meta$ancestry, t2$meta$ancestry), "This single-ancestry entry does not support cross-ancestry coloc; use a dedicated method")
assert(identical(t1$meta$region, t2$meta$region), "Two traits must describe the same regional interval")
assert(identical(t1$table$variant_id, t2$table$variant_id), "Two dense inputs must have exact variant coverage/order; harmonize upstream and audit missing causal mass")
assert(identical(t1$table$ea, t2$table$ea) && identical(t1$table$nea, t2$table$nea), "Two traits need identical effect/non-effect allele orientation")
if (any(abs(t1$table$eaf - t2$table$eaf) > 0.15)) warning("Large allele-frequency differences; inspect ancestry, reference orientation and ascertainment")
assert(is.data.frame(cfg$prior_grid) && all(c("p1", "p2", "p12") %in% names(cfg$prior_grid)), "Supply explicit prior_grid data frame")
assert(nrow(cfg$prior_grid) > 0 && all(is.finite(as.matrix(cfg$prior_grid[c("p1", "p2", "p12")]))), "Invalid priors")
assert(all(cfg$prior_grid$p1 > 0 & cfg$prior_grid$p1 < 1 & cfg$prior_grid$p2 > 0 & cfg$prior_grid$p2 < 1 & cfg$prior_grid$p12 > 0 & cfg$prior_grid$p12 <= pmin(cfg$prior_grid$p1, cfg$prior_grid$p2)), "p12 must be positive and <= min(p1,p2)")

fit_trait <- function(t) {
  z <- setNames(t$table$beta / t$table$se, t$table$variant_id)
  fit <- susieR::susie_rss(z = z, R = t$R, n = t$meta$N,
                         L = cfg$L, coverage = cfg$coverage,
                         max_iter = cfg$max_iter, estimate_residual_variance = FALSE)
  assert(isTRUE(fit$converged), "SuSiE did not converge; inspect LD and model settings before interpretation")
  # Public coloc helper adds identifiers required by coloc.susie.
  coloc::annotate_susie(fit, snp = t$table$variant_id, LD = t$R)
}
fit1 <- fit_trait(t1)
fit2 <- fit_trait(t2)
saveRDS(list(trait1 = fit1, trait2 = fit2), file.path(cfg$out_dir, "susie_fits.rds"))
write.table(data.frame(variant_id = t1$table$variant_id, pip_trait1 = fit1$pip, pip_trait2 = fit2$pip), file.path(cfg$out_dir, "variant_pip.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
to_dataset <- function(t) {
  d <- list(snp = t$table$variant_id, beta = t$table$beta, varbeta = t$table$se^2,
            MAF = pmin(t$table$eaf, 1 - t$table$eaf), N = t$meta$N, type = t$meta$type)
  if (t$meta$type == "quant") d$sdY <- t$meta$sdY else d$s <- t$meta$s
  d
}
all_summaries <- list()
all_abf <- list()
for (i in seq_len(nrow(cfg$prior_grid))) {
  priors <- as.list(cfg$prior_grid[i, c("p1", "p2", "p12")])
  abf <- do.call(coloc::coloc.abf, c(list(dataset1 = to_dataset(t1), dataset2 = to_dataset(t2)), priors))
  all_abf[[i]] <- cbind(prior_id = i, cfg$prior_grid[i, , drop = FALSE], as.data.frame(as.list(abf$summary)))
  if (length(fit1$sets$cs) > 0 && length(fit2$sets$cs) > 0) {
    coloc_result <- do.call(coloc::coloc.susie, c(list(dataset1 = fit1, dataset2 = fit2), priors))
    saveRDS(coloc_result, file.path(cfg$out_dir, paste0("coloc_susie_prior_", i, ".rds")))
    if (is.list(coloc_result) && !is.null(coloc_result$summary)) {
      all_summaries[[length(all_summaries) + 1L]] <- cbind(prior_id = i, cfg$prior_grid[rep(i, nrow(coloc_result$summary)), , drop = FALSE], as.data.frame(coloc_result$summary))
    }
  }
}
write.table(do.call(rbind, all_abf), file.path(cfg$out_dir, "single_causal_abf_baseline.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
status <- if (length(all_summaries)) "SIGNAL_PAIRS_ESTIMATED" else "NO_ESTIMABLE_SIGNAL_PAIRS"
if (length(all_summaries)) write.table(do.call(rbind, all_summaries), file.path(cfg$out_dir, "susie_coloc_h0_h4_prior_sensitivity.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
files <- c(cfg$trait1$summary, cfg$trait2$summary, cfg$trait1$ld_matrix, cfg$trait2$ld_matrix, cfg$trait1$ld_variants, cfg$trait2$ld_variants, args[[1]])
provenance <- list(status = status, config = cfg, input_md5 = tools::md5sum(files),
                   packages = list(susieR = as.character(packageVersion("susieR")), coloc = as.character(packageVersion("coloc"))),
                   session = sessionInfo(),
                   limitations = c("Single-ancestry dense-region analysis; upstream allele/build declarations must be independently verified", "Case-control RSS uses an explicit approximation when requested", "High H4 supports shared signal, not exposure mediation, drug effect or misuse", "No credible set is not proof of no biological effect", "Single-causal ABF is an assumption-limited comparison"))
saveRDS(provenance, file.path(cfg$out_dir, "run_provenance.rds"))
dput(provenance, file.path(cfg$out_dir, "run_provenance.txt"))
writeLines(status, file.path(cfg$out_dir, "status.txt"))
message(status)
