#!/usr/bin/env Rscript
# Native API behavior checks on synthetic summaries only. Not GC research results.
args0 <- commandArgs(FALSE)
this_file <- normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE)[1]), mustWork = TRUE)
repo <- dirname(dirname(this_file))
module <- file.path(repo, "modules/04_causal_inference_genetics/09_method_library/611_mr_coloc_library")
source(file.path(module, "run_native_mr.R")); source(file.path(module, "run_coloc_abf.R"))
native_require("jsonlite"); native_require("TwoSampleMR"); native_require("MendelianRandomization"); native_require("coloc")
scratch <- tempfile("mr-coloc-native-"); dir.create(scratch)
checks <- list(); mr_status <- NULL; mendelian_status <- NULL
check <- function(name, code) {
  failure <- tryCatch({ force(code); NULL }, error = function(e) conditionMessage(e))
  checks[[length(checks) + 1L]] <<- list(name = name, passed = is.null(failure), error = failure)
  cat(if (is.null(failure)) "PASS" else "FAIL", name, if (!is.null(failure)) paste0(": ", failure) else "", "\n")
}
reject <- function(expr, pattern) {
  e <- tryCatch({ force(expr); NULL }, error = function(e) conditionMessage(e))
  if (is.null(e) || !grepl(pattern, e, ignore.case = TRUE)) stop("Expected rejection matching ", pattern, "; got ", if (is.null(e)) "no error" else e)
}
fresh <- function(label) file.path(scratch, paste0(label, "_", length(list.files(scratch))))
write_input <- function(d, label) { p <- paste0(fresh(label), ".csv"); utils::write.csv(d, p, row.names = FALSE); p }
meta <- function(name, unit = "SD") list(construct = name, unit = unit, ancestry = "EUR", build = "GRCh37", source = "SYNTHETIC interface fixture, not participant data", sample_size = 100000)
set.seed(20261021)
n <- 10L
d <- data.frame(snp = paste0("rsSynthetic", seq_len(n)), beta_exposure = seq(.1, .32, length.out = n), se_exposure = rep(.015, n),
                beta_outcome = .3 * seq(.1, .32, length.out = n) + rnorm(n, sd = .012), se_outcome = rep(.025, n),
                effect_allele_exposure = "A", other_allele_exposure = "C", effect_allele_outcome = "A", other_allele_outcome = "C")
mr_input_path <- write_input(d, "mr")
mr_cfg <- function() list(engine = "TwoSampleMR", input = mr_input_path, methods = list("mr_ivw"), outdir = fresh("mr_result"), seed = 20261021,
                         harmonized = TRUE, independent_instruments = TRUE, metadata = list(exposure = meta("syntheticX"), outcome = meta("syntheticY"), sample_overlap = "Synthetic inputs; no real individuals"), parameters = list(nboot = 100))
check("TwoSampleMR installed method registry contains19 native APIs", { stopifnot(nrow(TwoSampleMR::mr_method_list()) == 19L) })
check("TwoSampleMR all18 multiSNP dispatches preserve missingRAPS as failed", {
  c <- mr_cfg(); c$methods <- as.list(setdiff(TwoSampleMR::mr_method_list()$obj, "mr_wald_ratio")); z <- run_native_mr(c)
  mr_status <<- z$results
  expected_fail <- if (requireNamespace("mr.raps", quietly = TRUE)) character() else "mr_raps"
  stopifnot(setequal(z$results$method[z$results$status == "failed"], expected_fail))
  stopifnot(file.exists(file.path(z$outdir, "native_results.rds")), file.exists(file.path(z$outdir, "provenance.json")))
  stopifnot(grepl("not_causal_beta", z$results$estimate_role[z$results$method == "mr_sign"]))
  r <- readRDS(file.path(z$outdir, "native_results.rds")); stopifnot(length(r$native_results) == 18L)
})
check("SingleSNP Wald native roundtrip matches installedAPI", {
  c <- mr_cfg(); c$input <- write_input(d[1, ], "single"); c$methods <- list("mr_wald_ratio"); z <- run_native_mr(c)
  oracle <- TwoSampleMR::mr_wald_ratio(d$beta_exposure[1], d$beta_outcome[1], d$se_exposure[1], d$se_outcome[1], TwoSampleMR::default_parameters())
  stopifnot(isTRUE(all.equal(z$results$estimate, unname(oracle$b))), isTRUE(all.equal(z$results$se, unname(oracle$se))))
})
check("MendelianRandomization11 exported native functions return reloadableS4", {
  c <- mr_cfg(); c$engine <- "MendelianRandomization"; c$methods <- lapply(native_mr_catalog(c$engine), `[[`, "id"); c$parameters <- NULL
  c$method_args <- list(mr_median = list(iterations = 100), mr_mbe = list(iterations = 100),
                       mr_conmix = list(CIMin = -2, CIMax = 2, CIStep = .02), mr_hetpen = list(CIMin = -2, CIMax = 2, CIStep = .02),
                       mr_pivw = list(Boot.Fieller = FALSE), mr_cML = list(DP = FALSE, MA = TRUE, maxit = 100))
  z <- run_native_mr(c); mendelian_status <<- z$results
  if (any(z$results$status != "returned")) stop(paste(z$results$method[z$results$status != "returned"], z$results$error[z$results$status != "returned"], collapse = ";"))
  r <- readRDS(file.path(z$outdir, "native_results.rds")); stopifnot(all(vapply(r$native_results, isS4, logical(1))), nrow(z$results) == 11L)
})
check("SignedLD accepted only for three native correlated summary methods", {
  correlation <- .2^abs(outer(seq_len(n), seq_len(n), "-")); ld <- data.frame(snp = d$snp, correlation, check.names = FALSE); names(ld)[-1] <- d$snp
  ldpath <- write_input(ld, "signed_ld")
  c <- mr_cfg(); c$engine <- "MendelianRandomization"; c$methods <- list("mr_ivw", "mr_egger", "mr_maxlik"); c$parameters <- NULL; c$independent_instruments <- FALSE
  c$ld <- list(file = ldpath, signed_r = TRUE, effect_allele_aligned = TRUE, reference = "syntheticpositive-definiteR", ancestry = "EUR", build = "GRCh37")
  z <- run_native_mr(c); stopifnot(all(z$results$status == "returned"))
  c$outdir <- fresh("bad_ld"); c$ld$signed_r <- FALSE; reject(run_native_mr(c), "signed_r")
  c$ld$signed_r <- TRUE
  shuffled <- ld[n:1, ]; c$ld$file <- write_input(shuffled, "shuffledLD"); reject(run_native_mr(c), "order")
  singular <- ld; singular[-1] <- 1; c$ld$file <- write_input(singular, "singularLD"); reject(run_native_mr(c), "positive definite")
  c$ld$file <- ldpath; c$engine <- "TwoSampleMR"; c$methods <- list("mr_ivw"); reject(run_native_mr(c), "does not support correlated")
})
for (case in c("zeroSE", "duplicate", "mismatchedAlleles", "zeroExposure", "singleEgger", "multipleWald", "unknownMethod", "notHarmonized", "notIndependent", "badBuild", "missingUnit", "ORinsteadOfBeta", "palindromeAmbiguous")) {
  local({ label <- case; check(paste("MRreject", label), {
    c <- mr_cfg(); bad <- d; pattern <- "."
    if (label == "zeroSE") { bad$se_exposure[1] <- 0; pattern <- "SE" }
    if (label == "duplicate") { bad$snp[2] <- bad$snp[1]; pattern <- "unique" }
    if (label == "mismatchedAlleles") { bad$effect_allele_outcome[1] <- "G"; pattern <- "not already aligned" }
    if (label == "zeroExposure") { bad$beta_exposure[1] <- 0; pattern <- "Zero exposure" }
    if (label == "singleEgger") { bad <- bad[1, ]; c$methods <- list("mr_egger_regression"); pattern <- "SNP count" }
    if (label == "multipleWald") { c$methods <- list("mr_wald_ratio"); pattern <- "SNP count" }
    if (label == "unknownMethod") { c$methods <- list("made_up_mr"); pattern <- "Unsupported method" }
    if (label == "notHarmonized") { c$harmonized <- FALSE; pattern <- "harmonized" }
    if (label == "notIndependent") { c$independent_instruments <- FALSE; pattern <- "independent" }
    if (label == "badBuild") { c$metadata$outcome$build <- "GRCh38"; pattern <- "build mismatch" }
    if (label == "missingUnit") { c$metadata$exposure$unit <- NULL; pattern <- "unit" }
    if (label == "ORinsteadOfBeta") { c$metadata$outcome$unit <- "OR"; pattern <- "log-odds beta" }
    if (label == "palindromeAmbiguous") { bad$other_allele_exposure[1] <- bad$other_allele_outcome[1] <- "T"; c$palindromes_resolved <- TRUE; bad$eaf_exposure <- bad$eaf_outcome <- .5; pattern <- "ambiguous" }
    c$input <- write_input(bad, label); reject(run_native_mr(c), pattern)
  }) })
}
check("Existing results are not overwritten", { c <- mr_cfg(); dir.create(c$outdir); writeLines("protected", file.path(c$outdir, "existing")); reject(run_native_mr(c), "must be empty"); stopifnot(readLines(file.path(c$outdir, "existing")) == "protected") })
check("RscriptJSON CLI propagates missingRAPS failure with persisted provenance", {
  if (!requireNamespace("mr.raps", quietly = TRUE)) {
    c <- mr_cfg(); c$methods <- list("mr_raps"); p <- paste0(fresh("config"), ".json"); native_write_json(c, p)
    rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
    z <- suppressWarnings(system2(rscript, c(shQuote(file.path(module, "run_native_mr.R")), "--config", shQuote(p)), stdout = TRUE, stderr = TRUE))
    stopifnot(identical(attr(z, "status"), 2L), file.exists(file.path(c$outdir, "provenance.json")))
    q <- jsonlite::read_json(file.path(c$outdir, "provenance.json")); stopifnot(q$status == "partial_failure", grepl("mr.raps", q$method_status[[1]]$error, fixed = TRUE))
  }
})

set.seed(20261021); k <- 120L; signal <- exp(-((seq_len(k) - 60) / 12)^2)
regional1 <- data.frame(snp = paste0("rsRegion", seq_len(k)), chromosome = "1", position = 1000000 + seq_len(k) * 10,
                        beta = .2 * signal + rnorm(k, sd = .01), se = .03, effect_allele = "A", other_allele = "C", eaf = seq(.05, .45, length.out = k))
regional2 <- regional1; regional2$beta <- .3 * signal + rnorm(k, sd = .015); regional2$se <- .05
r1 <- write_input(regional1, "region1"); r2 <- write_input(regional2, "region2")
co_cfg <- function() list(seed = 20261021, outdir = fresh("coloc"), alleles_aligned = TRUE, min_variants = 100,
                         region = list(chromosome = "1", start = 1000000, end = 1002000), metadata = list(sample_overlap = "synthetic-noindividuals"),
                         trait1 = list(name = "syntheticQuant", input = r1, type = "quant", N = 10000, sdY = 1, sdY_source = "syntheticstandardizedphenotype", full_region = TRUE, association_filtered = FALSE, unclumped = TRUE, metadata = meta("syntheticQuant")),
                         trait2 = list(name = "syntheticCaseControl", input = r2, type = "cc", N = 15000, case_fraction = .2, n_cases = 3000, full_region = TRUE, association_filtered = FALSE, unclumped = TRUE, metadata = meta("syntheticCaseControl", "logodds")),
                         prior_grid = list(list(p1 = 1e-4, p2 = 1e-4, p12 = 1e-5), list(p1 = 1e-4, p2 = 1e-4, p12 = 1e-6)))
check("colocABF densequant+cc nativeRDS and H0-H4 priorgrid", {
  c <- co_cfg(); z <- run_coloc_abf(c); stopifnot(nrow(z$results) == 2L, all(abs(rowSums(z$results[paste0("H", 0:4)]) - 1) < 1e-8))
  saved <- readRDS(file.path(z$outdir, "native_coloc.rds")); oracle <- coloc::coloc.abf(saved$native_datasets$trait1, saved$native_datasets$trait2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
  stopifnot(isTRUE(all.equal(oracle$summary, saved$native_results[[1]]$summary)), saved$native_datasets$trait2$s == .2)
})
check("colocABF explicitly allowed native sdY estimation", { c <- co_cfg(); c$trait1$sdY <- NULL; c$trait1$sdY_source <- NULL; c$trait1$estimate_sdY <- TRUE; z <- run_coloc_abf(c); stopifnot(z$status == "returned") })
check("Coloc intersection requires explicit acceptance and records exclusions", {
  c <- co_cfg(); c$trait2$input <- write_input(regional2[-1, ], "partial"); reject(run_coloc_abf(c), "intersection would drop")
  c$allow_partial_overlap <- TRUE; z <- run_coloc_abf(c); p <- jsonlite::read_json(file.path(z$outdir, "provenance.json")); stopifnot(p$input_rows$shared == 119, p$intersection$trait1_only[[1]] == regional1$snp[1])
})
for (case in c("leadOnly", "filteredTable", "zeroSE", "duplicate", "unaligned", "differentPosition", "wrongChromosome", "missingSdY", "missingCaseFraction", "invalidPrior", "smallOverlap")) {
  local({ label <- case; check(paste("colocReject", label), {
    c <- co_cfg(); bad <- regional2; pattern <- "."
    if (label == "leadOnly") { bad <- bad[1, ]; pattern <- "too few" }
    if (label == "filteredTable") { c$trait1$association_filtered <- TRUE; pattern <- "full_region" }
    if (label == "zeroSE") { bad$se[1] <- 0; pattern <- "SE" }
    if (label == "duplicate") { bad$snp[2] <- bad$snp[1]; pattern <- "unique" }
    if (label == "unaligned") { bad$effect_allele[1] <- "G"; pattern <- "alleles disagree" }
    if (label == "differentPosition") { bad$position[1] <- bad$position[1] + 1; pattern <- "position/alleles disagree" }
    if (label == "wrongChromosome") { bad$chromosome[1] <- "2"; pattern <- "Chromosome" }
    if (label == "missingSdY") { c$trait1$sdY <- NULL; pattern <- "Quant trait requires" }
    if (label == "missingCaseFraction") { c$trait2$case_fraction <- NULL; pattern <- "case_fraction" }
    if (label == "invalidPrior") { c$prior_grid[[1]]$p12 <- .01; pattern <- "prior" }
    if (label == "smallOverlap") { bad$snp[1:30] <- paste0("rsOther", 1:30); c$allow_partial_overlap <- TRUE; pattern <- "common regional" }
    c$trait2$input <- write_input(bad, label); reject(run_coloc_abf(c), pattern)
  }) })
}
receipt <- list(date = "2026-10-21", validation_scope = "SYNTHETIC native-interface/guard behavior only; no glucocorticoid or human scientific validation",
                R_version = R.version.string, packages = as.list(setNames(vapply(c("TwoSampleMR", "MendelianRandomization", "coloc", "jsonlite"), function(p) as.character(utils::packageVersion(p)), character(1)), c("TwoSampleMR", "MendelianRandomization", "coloc", "jsonlite"))),
                mr_raps_installed = requireNamespace("mr.raps", quietly = TRUE), checks = checks, tests_passed = sum(vapply(checks, `[[`, logical(1), "passed")), tests_total = length(checks),
                TwoSampleMR_multi_method_status = mr_status, MendelianRandomization_method_status = mendelian_status,
                code_md5 = native_md5(c(MR = file.path(module, "run_native_mr.R"), coloc = file.path(module, "run_coloc_abf.R"), behavior_test = this_file)))
receipt$status <- if (receipt$tests_passed == receipt$tests_total) "passed" else "failed"
dest <- file.path(repo, "docs/methods/mr_coloc/native_validation_20261021.json"); dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE); native_write_json(receipt, dest)
cat("Validation", receipt$status, receipt$tests_passed, "/", receipt$tests_total, "\n")
if (receipt$status != "passed") quit(status = 1L)
