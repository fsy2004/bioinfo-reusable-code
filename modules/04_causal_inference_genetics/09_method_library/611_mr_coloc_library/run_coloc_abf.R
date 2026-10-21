#!/usr/bin/env Rscript
# 官方 coloc::coloc.abf 单因果变异模型；无LD输入，不替代SuSiE模块609。
COLOC_ABF_SOURCE <- local({
  frames <- sys.frames(); files <- vapply(frames, function(x) if (is.null(x$ofile)) "" else x$ofile, character(1))
  if (any(nzchar(files))) normalizePath(tail(files[nzchar(files)], 1), mustWork = TRUE) else normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = TRUE)
})
source(file.path(dirname(COLOC_ABF_SOURCE), "run_native_mr.R"))

coloc_region_input <- function(spec, base, min_variants) {
  native_metadata(spec$metadata, "trait")
  native_nonempty(spec$name, "trait name")
  if (!isTRUE(spec$full_region) || !identical(spec$association_filtered, FALSE) || !isTRUE(spec$unclumped)) stop("Each trait must declare full_region=true, association_filtered=false, unclumped=true; no lead/significant-only table")
  if (!spec$type %in% c("quant", "cc")) stop("trait type must be quant or cc")
  p <- native_path(spec$input, base); d <- utils::read.csv(p, check.names = FALSE, stringsAsFactors = FALSE)
  if (anyDuplicated(names(d))) stop("Duplicate regional column names")
  if (!"snp" %in% names(d) || anyNA(d$snp) || any(!nzchar(d$snp)) || any(grepl("[[:space:]]", d$snp)) || anyDuplicated(d$snp)) stop("Regional SNP IDs must be nonmissing, whitespace-free and unique")
  d$snp <- as.character(d$snp)
  native_numeric_columns(d, c("position", "beta", "se", "eaf"), spec$name)
  if (any(d$se <= 0) || any(d$eaf <= 0 | d$eaf >= 1)) stop("Regional SE must be >0 and EAF strictly between0 and1")
  if (any(d$position <= 0 | d$position != floor(d$position))) stop("Regional positions must be positive integer coordinates")
  native_alleles(d, "effect_allele", "other_allele", spec$name)
  if (!"chromosome" %in% names(d) || anyNA(d$chromosome) || any(!nzchar(as.character(d$chromosome)))) stop("Regional chromosome column required")
  if (nrow(d) < min_variants) stop("Regional input has too few variants; lead-only input is not supported")
  N <- native_number(spec$N, "trait N")
  x <- list(snp = d$snp, position = d$position, beta = d$beta, varbeta = d$se^2, MAF = pmin(d$eaf, 1 - d$eaf), N = N, type = spec$type)
  if (spec$type == "cc") {
    fraction <- native_number(spec$case_fraction, "case_fraction")
    if (fraction >= 1) stop("case_fraction must be <1")
    if (!is.null(spec$n_cases) && abs(native_number(spec$n_cases, "n_cases") / N - fraction) > 1e-8) stop("case_fraction disagrees with n_cases/N")
    x$s <- fraction
  } else if (!is.null(spec[["sdY"]])) {
    if (isTRUE(spec[["estimate_sdY"]])) stop("Specify sdY or estimate_sdY, not both")
    x$sdY <- native_number(spec[["sdY"]], "sdY"); native_nonempty(spec[["sdY_source"]], "sdY_source")
  } else if (!isTRUE(spec[["estimate_sdY"]])) stop("Quant trait requires positive sdY + sdY_source or explicit estimate_sdY=true with N/MAF/varbeta")
  list(data = x, rows = d, path = p)
}

run_coloc_abf <- function(config) {
  parsed <- native_read_config(config); cfg <- parsed$value; base <- parsed$config_dir
  native_require("coloc")
  if (!isTRUE(cfg$alleles_aligned)) stop("alleles_aligned=true required; no automatic allele flips")
  seed <- native_number(cfg$seed, "seed", FALSE)
  if (seed < 0 || seed != floor(seed) || seed > .Machine$integer.max) stop("seed must be a nonnegative integer")
  min_variants <- if (is.null(cfg$min_variants)) 100 else native_number(cfg$min_variants, "min_variants")
  if (min_variants < 50 || min_variants != floor(min_variants)) stop("min_variants must be an integer >=50; this is an input guard, not a scientific power threshold")
  a <- coloc_region_input(cfg$trait1, base, min_variants); b <- coloc_region_input(cfg$trait2, base, min_variants)
  if (cfg$trait1$metadata$build != cfg$trait2$metadata$build || cfg$trait1$metadata$ancestry != cfg$trait2$metadata$ancestry) stop("Regional trait ancestry/build mismatch")
  native_nonempty(cfg$metadata$sample_overlap, "metadata.sample_overlap")
  native_nonempty(cfg$region$chromosome, "region.chromosome")
  if (any(as.character(a$rows$chromosome) != cfg$region$chromosome) || any(as.character(b$rows$chromosome) != cfg$region$chromosome)) stop("Chromosome disagrees with explicit region; no silent chromosome conversion")
  start <- native_number(cfg$region$start, "region.start"); end <- native_number(cfg$region$end, "region.end")
  if (end <= start || any(a$rows$position < start | a$rows$position > end) || any(b$rows$position < start | b$rows$position > end)) stop("Input coordinates outside explicit region")
  common <- intersect(a$rows$snp, b$rows$snp)
  fraction <- length(common) / max(nrow(a$rows), nrow(b$rows))
  different <- length(common) != nrow(a$rows) || length(common) != nrow(b$rows)
  if (different && !isTRUE(cfg$allow_partial_overlap)) stop("Regional SNP intersection would drop rows; explicit allow_partial_overlap=true required")
  minimum_fraction <- if (is.null(cfg$min_overlap_fraction)) .8 else native_number(cfg$min_overlap_fraction, "min_overlap_fraction")
  if (minimum_fraction > 1 || length(common) < min_variants || fraction < minimum_fraction) stop("Too little common regional coverage")
  ia <- match(common, a$rows$snp); ib <- match(common, b$rows$snp)
  x <- a$rows[ia, ]; y <- b$rows[ib, ]
  if (any(x$position != y$position | x$effect_allele != y$effect_allele | x$other_allele != y$other_allele)) stop("Shared regional SNP position/alleles disagree; harmonize explicitly before analysis")
  subset_dataset <- function(z, index) { for (n in c("snp", "position", "beta", "varbeta", "MAF")) z[[n]] <- z[[n]][index]; z }
  da <- subset_dataset(a$data, ia); db <- subset_dataset(b$data, ib)
  checks <- native_capture({ coloc::check_dataset(da); coloc::check_dataset(db) })
  if (inherits(checks$value, "native_failure")) stop(checks$value$message)
  priors <- cfg$prior_grid
  if (!is.list(priors) || !length(priors)) stop("Explicit prior_grid of p1/p2/p12 objects required")
  for (p in priors) {
    for (k in c("p1", "p2", "p12")) native_number(p[[k]], paste("prior", k))
    if (p$p12 > min(p$p1, p$p2) || any(unlist(p[c("p1", "p2", "p12")]) * length(common) >= 1)) stop("Invalid or too-large prior; prevent native silent prior adjustment")
  }
  outdir <- native_output_dir(cfg$outdir, base); set.seed(as.integer(seed))
  files <- c(trait1 = a$path, trait2 = b$path, if (!is.null(parsed$config_file)) c(config = parsed$config_file))
  prov <- native_provenance(cfg, files, c("coloc", "jsonlite"), c(coloc_adapter = COLOC_ABF_SOURCE, shared_validation = MR_NATIVE_SOURCE))
  prov$input_rows <- c(trait1 = nrow(a$rows), trait2 = nrow(b$rows), shared = length(common))
  prov$intersection <- list(trait1_only = setdiff(a$rows$snp, common), trait2_only = setdiff(b$rows$snp, common), fraction_trait1 = length(common) / nrow(a$rows), fraction_trait2 = length(common) / nrow(b$rows), minimum_coverage_fraction = fraction, accepted_explicitly = !different || isTRUE(cfg$allow_partial_overlap), common_order = "trait1 input order; trait2 explicitly matched by SNP with coordinate+allele checks")
  prov$assumptions <- c("At most one causal variant per trait in the prespecified full region", "Unfiltered regional coverage and appropriate trait variance/sample metadata", "Prespecified prior sensitivity; H4 is not proof of mediation, drug action or causal direction", "No LD supplied/required by coloc.abf; allelic heterogeneity requires a different model such as module609 SuSiE")
  fitted <- list(); rows <- list(); warnings <- list()
  for (i in seq_along(priors)) {
    p <- priors[[i]]
    fit <- native_capture(coloc::coloc.abf(da, db, p1 = p$p1, p2 = p$p2, p12 = p$p12))
    if (inherits(fit$value, "native_failure")) stop("coloc.abf failed: ", fit$value$message)
    fitted[[i]] <- fit$value
    summary <- fit$value$summary
    if (any(!is.finite(summary)) || abs(sum(summary[grep("^PP", names(summary))]) - 1) > 1e-8) stop("Native posterior summary invalid")
    rows[[i]] <- data.frame(prior_index = i, p1 = p$p1, p2 = p$p2, p12 = p$p12,
                            nsnps = unname(summary["nsnps"]), H0 = unname(summary["PP.H0.abf"]), H1 = unname(summary["PP.H1.abf"]), H2 = unname(summary["PP.H2.abf"]), H3 = unname(summary["PP.H3.abf"]), H4 = unname(summary["PP.H4.abf"]))
    warnings[[i]] <- unique(c(checks$warnings, fit$warnings))
  }
  prov$status <- "returned"; prov$warnings <- warnings
  results <- do.call(rbind, rows)
  saveRDS(list(native_datasets = list(trait1 = da, trait2 = db), native_results = fitted, provenance = prov), file.path(outdir, "native_coloc.rds"))
  utils::write.csv(results, file.path(outdir, "posterior_h0_h4.csv"), row.names = FALSE)
  for (i in seq_along(fitted)) utils::write.csv(fitted[[i]]$results, file.path(outdir, paste0("native_variant_results_prior_", i, ".csv")), row.names = FALSE)
  native_write_json(prov, file.path(outdir, "provenance.json")); writeLines(utils::capture.output(utils::sessionInfo()), file.path(outdir, "sessionInfo.txt"))
  list(status = "returned", outdir = outdir, results = results, native = fitted)
}
if (sys.nframe() == 0L) native_cli(run_coloc_abf)
