#!/usr/bin/env Rscript
# 显式配置的原生 MR API 适配器；不下载、不安装、不自动 harmonize。
# Source: installed TwoSampleMR::mr_method_list/default_parameters/native methods,
# MendelianRandomization::mr_input and exported univariable functions.

native_adapter_path <- function() {
  frames <- sys.frames()
  files <- vapply(frames, function(x) if (is.null(x$ofile)) "" else x$ofile, character(1))
  if (any(nzchar(files))) return(normalizePath(tail(files[nzchar(files)], 1), mustWork = TRUE))
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (!length(arg)) stop("Cannot locate adapter source path")
  normalizePath(sub("^--file=", "", arg[[1]]), mustWork = TRUE)
}
MR_NATIVE_SOURCE <- native_adapter_path()

native_require <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) stop("Missing installed package: ", package, "; no automatic installation")
}
native_nonempty <- function(x, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(trimws(x))) stop(name, " must be a nonempty string")
  x
}
native_number <- function(x, name, positive = TRUE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || (positive && x <= 0)) stop(name, " must be finite", if (positive) " and positive" else "")
  x
}
native_path <- function(x, config_dir, must_work = TRUE) {
  native_nonempty(x, "path")
  absolute <- grepl("^([A-Za-z]:[/\\\\]|/|\\\\\\\\)", x)
  p <- if (absolute) x else file.path(config_dir, x)
  normalizePath(p, winslash = "/", mustWork = must_work)
}
native_metadata <- function(x, label) {
  if (!is.list(x)) stop(label, " metadata must be an object")
  for (field in c("construct", "unit", "ancestry", "build", "source")) native_nonempty(x[[field]], paste(label, field))
  if (tolower(trimws(x$unit)) %in% c("or", "odds ratio", "odds_ratio")) stop(label, ": supply log-odds beta, not an odds ratio")
  if (!(x$build %in% c("GRCh37", "GRCh38"))) stop(label, " build must be GRCh37 or GRCh38")
  invisible(x)
}
native_read_config <- function(config) {
  native_require("jsonlite")
  if (is.character(config)) {
    p <- normalizePath(config, mustWork = TRUE)
    x <- jsonlite::fromJSON(p, simplifyVector = FALSE)
    return(list(value = x, config_dir = dirname(p), config_file = p))
  }
  if (!is.list(config)) stop("config must be a JSON filename or list")
  list(value = config, config_dir = getwd(), config_file = NULL)
}
native_output_dir <- function(x, base) {
  p <- native_path(x, base, FALSE)
  if (dir.exists(p) && length(list.files(p, all.files = TRUE, no.. = TRUE))) stop("Output directory must be empty: ", p)
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(p)) stop("Cannot create output directory")
  p
}
native_write_json <- function(x, p) jsonlite::write_json(x, p, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null", digits = 16)
native_md5 <- function(paths) as.list(stats::setNames(unname(tools::md5sum(paths)), names(paths)))
native_provenance <- function(cfg, files, packages, code_files) {
  list(schema_version = "1.0", timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
       configuration = cfg, input_md5 = native_md5(files), code_md5 = native_md5(code_files),
       R_version = R.version.string, platform = R.version$platform,
       packages = as.list(stats::setNames(vapply(packages, function(x) if (requireNamespace(x, quietly = TRUE)) as.character(utils::packageVersion(x)) else "NOT_INSTALLED", character(1)), packages)),
       scientific_validity = "Native execution and declared-input checks do not establish IV validity, drug action, causal identification or scientific validation")
}
native_capture <- function(expr) {
  warnings <- character()
  value <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  }), error = function(e) structure(list(message = conditionMessage(e)), class = "native_failure"))
  list(value = value, warnings = unique(warnings))
}
native_flatten <- function(x, prefix = "") {
  if (isS4(x)) x <- stats::setNames(lapply(methods::slotNames(x), function(n) methods::slot(x, n)), methods::slotNames(x))
  if (is.list(x)) {
    out <- lapply(seq_along(x), function(i) native_flatten(x[[i]], paste0(prefix, if (nzchar(prefix)) "." else "", if (!is.null(names(x)) && nzchar(names(x)[i])) names(x)[i] else i)))
    return(if (length(out)) do.call(rbind, out) else data.frame(field = character(), value = character()))
  }
  if (!is.atomic(x) || !length(x)) return(data.frame(field = character(), value = character()))
  data.frame(field = paste0(prefix, if (length(x) > 1) paste0("[", seq_along(x), "]") else ""), value = as.character(x), stringsAsFactors = FALSE)
}
native_scalar <- function(x, fields) {
  for (f in fields) {
    value <- if (isS4(x) && f %in% methods::slotNames(x)) methods::slot(x, f) else if (is.list(x)) x[[f]] else NULL
    if (is.numeric(value) && length(value) == 1L) return(as.numeric(value))
  }
  NA_real_
}
native_alleles <- function(d, ea, oa, label) {
  for (x in c(ea, oa)) if (!x %in% names(d) || anyNA(d[[x]]) || any(!grepl("^[ACGT]$", d[[x]]))) stop(label, ": require uppercase single-base A/C/G/T alleles in ", x)
  if (any(d[[ea]] == d[[oa]])) stop(label, ": effect and other allele are identical")
}
native_numeric_columns <- function(d, fields, label) {
  for (f in fields) if (!f %in% names(d) || !is.numeric(d[[f]]) || any(!is.finite(d[[f]]))) stop(label, ": missing/nonfinite numeric column ", f)
}
native_ld <- function(spec, ids, base, metadata) {
  if (!is.list(spec) || !isTRUE(spec$signed_r) || !isTRUE(spec$effect_allele_aligned)) stop("LD requires signed_r=true and effect_allele_aligned=true; r2 is not a signed correlation matrix")
  for (f in c("reference", "ancestry", "build")) native_nonempty(spec[[f]], paste("LD", f))
  if (spec$build != metadata$exposure$build || spec$ancestry != metadata$exposure$ancestry) stop("LD ancestry/build mismatch")
  p <- native_path(spec$file, base)
  d <- utils::read.csv(p, check.names = FALSE, stringsAsFactors = FALSE)
  if (anyDuplicated(names(d))) stop("Duplicate LD column names")
  if (!identical(names(d)[1], "snp") || !identical(d[[1]], ids) || !identical(names(d)[-1], ids)) stop("Signed LD SNP row/column order must exactly match input; no silent reordering")
  if (any(!vapply(d[-1], is.numeric, logical(1)))) stop("LD entries must be numeric")
  a <- as.matrix(d[-1]); rownames(a) <- ids
  if (nrow(a) != ncol(a) || any(!is.finite(a)) || max(abs(a - t(a))) > 1e-8 || max(abs(diag(a) - 1)) > 1e-8 || any(abs(a) > 1 + 1e-8)) stop("LD must be finite square symmetric correlation with unit diagonal")
  if (min(eigen(a, symmetric = TRUE, only.values = TRUE)$values) <= 1e-8) stop("LD must be positive definite; resolve singular SNP sets explicitly")
  list(matrix = a, path = p)
}

native_mr_catalog <- function(engine) {
  if (engine == "TwoSampleMR") {
    native_require("TwoSampleMR")
    ids <- TwoSampleMR::mr_method_list()$obj
    mins <- stats::setNames(rep(3L, length(ids)), ids)
    mins[c("mr_wald_ratio", "mr_two_sample_ml", "mr_ivw", "mr_ivw_radial", "mr_ivw_mre", "mr_ivw_fe", "mr_uwr")] <- c(1L, 2L, 2L, 2L, 2L, 2L, 2L)
    mins["mr_sign"] <- 6L
    return(lapply(ids, function(n) list(id = n, fn = n, min_snps = unname(mins[n]), signed_ld = FALSE)))
  }
  if (engine == "MendelianRandomization") {
    ids <- c("mr_ivw", "mr_egger", "mr_median", "mr_maxlik", "mr_mbe", "mr_conmix", "mr_lasso", "mr_divw", "mr_pivw", "mr_hetpen", "mr_cML")
    return(lapply(ids, function(n) list(id = n, fn = n, min_snps = if (n %in% c("mr_ivw", "mr_maxlik")) 2L else 3L,
                                     signed_ld = n %in% c("mr_ivw", "mr_egger", "mr_maxlik"))))
  }
  stop("engine must be TwoSampleMR or MendelianRandomization")
}

run_native_mr <- function(config) {
  parsed <- native_read_config(config); cfg <- parsed$value; base <- parsed$config_dir
  native_nonempty(cfg$engine, "engine"); catalog <- native_mr_catalog(cfg$engine)
  native_require(cfg$engine)
  native_metadata(cfg$metadata$exposure, "exposure"); native_metadata(cfg$metadata$outcome, "outcome")
  native_nonempty(cfg$metadata$sample_overlap, "metadata.sample_overlap")
  if (cfg$metadata$exposure$build != cfg$metadata$outcome$build) stop("Exposure/outcome build mismatch")
  if (cfg$metadata$exposure$ancestry != cfg$metadata$outcome$ancestry) stop("Exposure/outcome ancestry mismatch; this adapter requires matched declared ancestry")
  if (!isTRUE(cfg$harmonized)) stop("harmonized=true required; adapter does not harmonize alleles")
  seed <- native_number(cfg$seed, "seed", FALSE)
  if (seed != floor(seed) || seed < 0 || seed > .Machine$integer.max) stop("seed must be a nonnegative R integer")
  ids <- unlist(cfg$methods, use.names = FALSE)
  if (!is.character(ids) || !length(ids) || anyNA(ids) || anyDuplicated(ids)) stop("methods must be nonempty unique native method IDs")
  choices <- stats::setNames(catalog, vapply(catalog, `[[`, character(1), "id"))
  if (any(!ids %in% names(choices))) stop("Unsupported method: ", paste(setdiff(ids, names(choices)), collapse = ", "))
  path <- native_path(cfg$input, base); d <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (anyDuplicated(names(d))) stop("Duplicate summary column names")
  if (!"snp" %in% names(d) || !nrow(d) || anyNA(d$snp) || any(!nzchar(d$snp)) || any(grepl("[[:space:]]", d$snp)) || anyDuplicated(d$snp)) stop("Require nonempty unique snp identifiers without whitespace; duplicates cannot be silently removed")
  d$snp <- as.character(d$snp)
  native_numeric_columns(d, c("beta_exposure", "se_exposure", "beta_outcome", "se_outcome"), "MR")
  if (any(d$se_exposure <= 0 | d$se_outcome <= 0)) stop("SE must be strictly positive")
  if (any(d$beta_exposure == 0)) stop("Zero exposure beta is not an informative instrument")
  native_alleles(d, "effect_allele_exposure", "other_allele_exposure", "exposure")
  native_alleles(d, "effect_allele_outcome", "other_allele_outcome", "outcome")
  if (any(d$effect_allele_exposure != d$effect_allele_outcome | d$other_allele_exposure != d$other_allele_outcome)) stop("Alleles are not already aligned; no automatic strand/allele flips")
  for (f in intersect(c("eaf_exposure", "eaf_outcome"), names(d))) {
    native_numeric_columns(d, f, "EAF")
    if (any(d[[f]] <= 0 | d[[f]] >= 1)) stop("EAF must be strictly between0 and1")
  }
  pal <- paste0(d$effect_allele_exposure, d$other_allele_exposure) %in% c("AT", "TA", "CG", "GC")
  if (any(pal)) {
    if (!isTRUE(cfg$palindromes_resolved)) stop("Palindromic variants require explicit palindromes_resolved=true and both EAF columns")
    native_numeric_columns(d, c("eaf_exposure", "eaf_outcome"), "palindromes")
    if (any(d$eaf_exposure <= 0 | d$eaf_exposure >= 1 | d$eaf_outcome <= 0 | d$eaf_outcome >= 1) || any(abs(d$eaf_exposure[pal] - .5) < .08 | abs(d$eaf_outcome[pal] - .5) < .08) || any(abs(d$eaf_exposure[pal] - d$eaf_outcome[pal]) > .1)) stop("Palindromic allele frequency is ambiguous or inconsistent")
  }
  for (n in ids) {
    spec <- choices[[n]]
    if (nrow(d) < spec$min_snps || (n == "mr_wald_ratio" && nrow(d) != 1L)) stop(n, " incompatible SNP count: ", nrow(d))
    if (n == "mr_hetpen" && nrow(d) > 15L) stop("mr_hetpen exhaustive subsets capped at15 SNPs in adapter; use an explicitly reviewed alternative")
  }
  ld <- NULL
  if (!is.null(cfg$ld)) {
    if (isTRUE(cfg$independent_instruments)) stop("Do not declare independent_instruments=true together with signed LD")
    if (cfg$engine != "MendelianRandomization" || any(!vapply(choices[ids], `[[`, logical(1), "signed_ld"))) stop("Selected engine/method does not support correlated instruments")
    ld <- native_ld(cfg$ld, d$snp, base, cfg$metadata)
  } else if (!isTRUE(cfg$independent_instruments)) stop("independent_instruments=true required when no signed LD supplied")
  outdir <- native_output_dir(cfg$outdir, base)
  packages <- c(cfg$engine, "jsonlite", if ("mr_raps" %in% ids) "mr.raps")
  files <- c(summary = path, if (!is.null(ld)) c(signed_ld = ld$path), if (!is.null(parsed$config_file)) c(config = parsed$config_file))
  prov <- native_provenance(cfg, files, packages, c(adapter = MR_NATIVE_SOURCE))
  prov$input_rows <- nrow(d); prov$methods_requested <- ids; prov$LD <- if (is.null(ld)) "independence declared by user, not established by adapter" else "signed correlation validated; source/allele alignment are user declarations"
  set.seed(as.integer(seed)); fitted <- list(); rows <- list(); flat <- list()
  object <- NULL
  if (cfg$engine == "MendelianRandomization") {
    args <- list(bx = d$beta_exposure, bxse = d$se_exposure, by = d$beta_outcome, byse = d$se_outcome,
                 snps = d$snp, effect_allele = d$effect_allele_exposure, other_allele = d$other_allele_exposure,
                 exposure = cfg$metadata$exposure$construct, outcome = cfg$metadata$outcome$construct)
    if (!is.null(ld)) args$correlation <- ld$matrix
    object <- do.call(MendelianRandomization::mr_input, args)
  }
  for (id in ids) {
    captured <- native_capture({
      if (cfg$engine == "TwoSampleMR") {
        if (id == "mr_raps") native_require("mr.raps")
        parameters <- TwoSampleMR::default_parameters()
        user <- cfg$parameters
        if (!is.null(user)) {
          if (!is.list(user) || any(!names(user) %in% names(parameters))) stop("Unknown TwoSampleMR parameter")
          for (n in names(user)) parameters[[n]] <- user[[n]]
        }
        do.call(getExportedValue("TwoSampleMR", id), list(b_exp = d$beta_exposure, b_out = d$beta_outcome, se_exp = d$se_exposure, se_out = d$se_outcome, parameters = parameters))
      } else {
        f <- getExportedValue("MendelianRandomization", id); formal <- names(formals(f)); user <- cfg$method_args[[id]]
        if (is.null(user)) user <- list()
        if (!is.list(user) || any(!names(user) %in% formal) || any(names(user) %in% c("object", "correl", "seed", "random_seed", "n"))) stop("Unknown or reserved method_args for ", id)
        if (!is.null(ld) && (isTRUE(user$robust) || isTRUE(user$penalized))) stop("Correlated native IVW/Egger adapters do not combine robust/penalized options with LD")
        args <- c(list(object = object), user)
        if ("correl" %in% formal) args$correl <- !is.null(ld)
        if ("seed" %in% formal) args$seed <- as.integer(seed)
        if ("random_seed" %in% formal) args$random_seed <- as.integer(seed)
        if (id == "mr_cML") args$n <- native_number(cfg$metadata$exposure$sample_size, "cML exposure sample_size")
        do.call(f, args)
      }
    })
    value <- captured$value; failed <- inherits(value, "native_failure")
    estimate <- if (failed) NA_real_ else native_scalar(value, c("b", "Estimate"))
    if (!failed && !is.finite(estimate)) { failed <- TRUE; value <- structure(list(message = "Native method returned a nonfinite estimate"), class = "native_failure") }
    fitted[[id]] <- captured$value
    row <- data.frame(engine = cfg$engine, method = id, status = if (failed) "failed" else "returned",
                      estimate_role = if (id == "mr_sign") "sign_concordance_score_not_causal_beta" else "effect_on_declared_native_scale",
                      estimate = estimate, se = if (failed) NA_real_ else native_scalar(value, c("se", "StdError", "StdError.Est")),
                      pvalue = if (failed) NA_real_ else native_scalar(value, c("pval", "Pvalue", "Pvalue.Est")),
                      ci_lower = if (failed) NA_real_ else native_scalar(value, c("CILower", "CILower.Est")),
                      ci_upper = if (failed) NA_real_ else native_scalar(value, c("CIUpper", "CIUpper.Est")),
                      input_snps = nrow(d), nsnp = if (failed) NA_real_ else native_scalar(value, c("nsnp", "SNPs")),
                      error = if (failed) value$message else "", warnings = paste(captured$warnings, collapse = " | "), stringsAsFactors = FALSE)
    rows[[id]] <- row
    if (!failed) { detail <- native_flatten(value); detail$method <- id; flat[[id]] <- detail }
  }
  results <- do.call(rbind, rows); rownames(results) <- NULL
  prov$method_status <- lapply(rows, function(x) list(method = x$method, status = x$status, error = x$error, warnings = x$warnings))
  prov$status <- if (any(results$status == "failed")) "partial_failure" else "returned"
  saveRDS(list(native_input = object, native_results = fitted, provenance = prov), file.path(outdir, "native_results.rds"))
  utils::write.csv(results, file.path(outdir, "results.csv"), row.names = FALSE, na = "")
  if (length(flat)) utils::write.csv(do.call(rbind, flat), file.path(outdir, "native_fields.csv"), row.names = FALSE, na = "")
  native_write_json(prov, file.path(outdir, "provenance.json"))
  writeLines(utils::capture.output(utils::sessionInfo()), file.path(outdir, "sessionInfo.txt"))
  list(status = prov$status, results = results, outdir = outdir, native = fitted)
}

native_cli <- function(runner) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2L || args[[1]] != "--config") { cat("Usage: Rscript <adapter.R> --config <explicit.json>\n", file = stderr()); quit(status = 2L) }
  code <- tryCatch({ x <- runner(args[[2]]); cat("Native adapter status:", x$status, "\nOutputs:", x$outdir, "\n"); if (x$status == "returned") 0L else 2L }, error = function(e) { cat("Adapter rejected/failed:", conditionMessage(e), "\n", file = stderr()); 2L })
  quit(status = code)
}
if (sys.nframe() == 0L) native_cli(run_native_mr)
