# Original missingness, multiple-imputation and propensity-weight diagnostics.
run_missingness_confounding <- function(config) {
  required <- c("input", "output_dir", "mode", "id")
  if (!all(required %in% names(config))) stop("Config requires: ", paste(required, collapse = ", "))
  if (!config$mode %in% c("audit", "mi_glm", "ps_diagnostics")) stop("mode must be audit, mi_glm or ps_diagnostics")
  dat <- read.csv(config$input, stringsAsFactors = FALSE, check.names = FALSE)
  if (anyDuplicated(names(dat)) || !config$id %in% names(dat) || anyNA(dat[[config$id]]) || anyDuplicated(dat[[config$id]])) stop("One row per nonmissing ID and unique columns required")
  factors <- unlist(config$factor_variables, use.names = FALSE)
  if (!all(factors %in% names(dat))) stop("Unknown factor_variables")
  dat[factors] <- lapply(dat[factors], factor)
  missingness <- data.frame(variable = names(dat), missing_n = vapply(dat, function(x) sum(is.na(x)), integer(1)),
                            observed_n = vapply(dat, function(x) sum(!is.na(x)), integer(1)),
                            missing_fraction = vapply(dat, function(x) mean(is.na(x)), numeric(1)), row.names = NULL)
  results <- list()
  if (config$mode != "audit") {
    if (!all(c("exposure", "confounders", "dag_assumptions") %in% names(config)) || !nzchar(config$dag_assumptions)) stop("Specify exposure, DAG-justified confounders and dag_assumptions")
    if (!is.null(config$survey_weight) || isTRUE(config$survey_design)) stop("These MI/PS adapters do not combine survey design/weights; choose a separately justified analysis")
    confounders <- unlist(config$confounders, use.names = FALSE)
    if (!length(confounders) || config$exposure %in% confounders || config$id %in% confounders) stop("Require an explicit nonempty pre-exposure adjustment set excluding ID and exposure")
    vars <- c(config$exposure, confounders)
    if (!all(vars %in% names(dat)) || any(!grepl("^[A-Za-z.][A-Za-z0-9._]*$", vars))) stop("Model variables missing or nonsyntactic")
  }
  if (config$mode == "mi_glm") {
    for (pkg in c("mice", "broom")) if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing MI dependency: ", pkg)
    mi_required <- c("outcome", "imputation_variables", "imputation_methods", "m", "maxit", "seed", "family")
    if (!all(mi_required %in% names(config))) stop("MI config requires: ", paste(mi_required, collapse = ", "))
    ivars <- unlist(config$imputation_variables, use.names = FALSE)
    if (anyDuplicated(ivars) || config$id %in% ivars || !all(c(vars, config$outcome) %in% ivars) || !all(ivars %in% names(dat))) stop("imputation_variables must include model variables and exclude participant ID")
    if (config$outcome %in% confounders || config$outcome == config$exposure) stop("Outcome cannot also be exposure or a confounder")
    if (any(!grepl("^[A-Za-z.][A-Za-z0-9._]*$", ivars))) stop("Use syntactic imputation column names")
    x <- dat[ivars]
    if (any(vapply(x, is.character, logical(1)))) stop("Declare character columns in factor_variables or supply numeric encoding")
    m <- as.integer(config$m); maxit <- as.integer(config$maxit)
    if (length(m) != 1L || is.na(m) || m < 2L || length(maxit) != 1L || is.na(maxit) || maxit < 1L) stop("m >= 2 and maxit >= 1 required")
    methods <- setNames(rep("", length(ivars)), ivars)
    supplied <- unlist(config$imputation_methods, use.names = TRUE)
    if (!length(supplied) || is.null(names(supplied)) || !all(names(supplied) %in% ivars)) stop("imputation_methods must be a named variable-to-method map")
    methods[names(supplied)] <- supplied
    incomplete <- names(x)[vapply(x, anyNA, logical(1))]
    if (any(!nzchar(methods[incomplete]))) stop("Every incomplete variable requires an explicitly selected method")
    predictors <- matrix(1L, length(ivars), length(ivars), dimnames = list(ivars, ivars)); diag(predictors) <- 0L
    if (!is.null(config$predictor_matrix)) {
      predictors <- as.matrix(read.csv(config$predictor_matrix, row.names = 1, check.names = FALSE))
      if (!identical(rownames(predictors), ivars) || !identical(colnames(predictors), ivars) || anyNA(predictors) || !all(predictors %in% c(0, 1)) || any(diag(predictors) != 0)) stop("Predictor matrix must match ordered imputation_variables, use 0/1 and zero diagonal")
    }
    if (!config$family %in% c("binomial", "gaussian")) stop("family must be binomial or gaussian")
    if (config$family == "binomial" && !all(is.na(x[[config$outcome]]) | x[[config$outcome]] %in% c(0, 1))) stop("Binary outcome must be 0/1/NA")
    imputed <- mice::mice(x, m = m, maxit = maxit, method = methods, predictorMatrix = predictors, seed = as.integer(config$seed), printFlag = FALSE)
    formula <- stats::reformulate(c(config$exposure, confounders), response = config$outcome)
    fam <- if (config$family == "binomial") stats::binomial() else stats::gaussian()
    # Fit each completed data set, then pool parameter estimates using Rubin's rules.
    fits <- lapply(seq_len(m), function(i) stats::glm(formula, data = mice::complete(imputed, i), family = fam, na.action = na.fail))
    pooled <- mice::pool(mice::as.mira(fits))
    results <- list(imputations = imputed, fits = fits, pooled = pooled)
    table <- summary(pooled, conf.int = TRUE)
    explanation <- "MI GLM under specified imputation/analysis and MAR assumptions; pool estimates, never average completed data sets. No complex-survey design correction. Inspect chain traces, loggedEvents, compatibility and MNAR sensitivity separately."
  }
  if (config$mode == "ps_diagnostics") {
    if (anyNA(dat[vars])) stop("PS inputs contain missingness; explicitly resolve it before estimating a propensity score")
    z <- dat[[config$exposure]]
    if (!is.numeric(z) || !all(z %in% c(0, 1)) || length(unique(z)) != 2L) stop("PS adapter requires a numeric binary 0/1 exposure with both groups")
    form <- stats::reformulate(confounders, response = config$exposure)
    model <- stats::glm(form, data = dat, family = binomial())
    if (!isTRUE(model$converged) || anyNA(stats::coef(model))) stop("PS model not converged or rank deficient")
    ps <- stats::predict(model, type = "response")
    if (any(!is.finite(ps) | ps <= 0 | ps >= 1)) stop("Propensity scores at numerical 0/1: inspect separation and positivity")
    if (is.null(config$estimand) || !config$estimand %in% c("ATE", "ATT", "ATO")) stop("Prespecify estimand: ATE, ATT or ATO")
    if (is.null(config$positivity_bounds) || length(config$positivity_bounds) != 2L) stop("Specify positivity_bounds as two diagnostic thresholds; no automatic trimming is performed")
    bounds <- as.numeric(config$positivity_bounds)
    if (any(!is.finite(bounds)) || bounds[1] <= 0 || bounds[2] >= 1 || bounds[1] >= bounds[2]) stop("Require 0 < lower < upper < 1")
    prevalence <- mean(z)
    w <- switch(config$estimand, ATE = ifelse(z == 1, prevalence / ps, (1 - prevalence) / (1 - ps)),
                ATT = ifelse(z == 1, 1, ps / (1 - ps)), ATO = ifelse(z == 1, 1 - ps, ps))
    # Fixed unweighted pooled SD makes before/after SMD directly comparable.
    mm <- stats::model.matrix(stats::reformulate(confounders), data = dat)
    mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
    weighted_mean <- function(x, weights) sum(x * weights) / sum(weights)
    balance <- lapply(seq_len(ncol(mm)), function(j) {
      x <- mm[, j]
      scale <- sqrt((stats::var(x[z == 1]) + stats::var(x[z == 0])) / 2)
      before <- mean(x[z == 1]) - mean(x[z == 0])
      after <- weighted_mean(x[z == 1], w[z == 1]) - weighted_mean(x[z == 0], w[z == 0])
      data.frame(covariate = colnames(mm)[j], smd_unweighted = if (is.finite(scale) && scale > 0) before / scale else NA_real_,
                 smd_weighted = if (is.finite(scale) && scale > 0) after / scale else NA_real_, denominator = scale)
    })
    table <- do.call(rbind, balance)
    diagnostics <- do.call(rbind, lapply(0:1, function(g) {
      p <- ps[z == g]; wg <- w[z == g]
      data.frame(exposure = g, n = length(p), ps_min = min(p), ps_max = max(p),
                 outside_bounds_n = sum(p < bounds[1] | p > bounds[2]), max_weight = max(wg),
                 effective_sample_size = sum(wg)^2 / sum(wg^2))
    }))
    overlap <- c(lower = max(min(ps[z == 0]), min(ps[z == 1])), upper = min(max(ps[z == 0]), max(ps[z == 1])))
    if (overlap[1] >= overlap[2]) warning("No empirical score-range overlap; causal comparison is unsupported")
    results <- list(model = model, balance = table, diagnostics = diagnostics, empirical_overlap = overlap,
                    record_weights = data.frame(person_id = dat[[config$id]], exposure = z, propensity_score = ps, analysis_weight = w))
    explanation <- "Diagnostic propensity weights only; no outcome treatment effect or survey-weight multiplication is performed. Positivity is not guaranteed by a nonzero fitted score. Balance uses a fixed unweighted pooled SD; factor levels are separate indicator columns. ATE uses stabilized weights; ATT and ATO are their explicit target-population weights."
  }
  if (dir.exists(config$output_dir) && length(list.files(config$output_dir, all.files = TRUE, no.. = TRUE))) stop("Output directory must be empty")
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(missingness, file.path(config$output_dir, "missingness.csv"), row.names = FALSE)
  if (config$mode != "audit") {
    write.csv(table, file.path(config$output_dir, if (config$mode == "mi_glm") "mi_pooled_coefficients.csv" else "covariate_balance.csv"), row.names = FALSE)
    saveRDS(results, file.path(config$output_dir, "results.rds"))
    if (config$mode == "ps_diagnostics") {
      write.csv(results$diagnostics, file.path(config$output_dir, "positivity_ess.csv"), row.names = FALSE)
      write.csv(results$record_weights, file.path(config$output_dir, "record_weights.csv"), row.names = FALSE)
      write.csv(data.frame(lower = results$empirical_overlap[1], upper = results$empirical_overlap[2]), file.path(config$output_dir, "empirical_score_overlap.csv"), row.names = FALSE)
    } else {
      write.csv(imputed$loggedEvents, file.path(config$output_dir, "mice_logged_events.csv"), row.names = FALSE)
      saveRDS(imputed$chainMean, file.path(config$output_dir, "mice_chain_means.rds"))
    }
    writeLines(explanation, file.path(config$output_dir, "interpretation.txt"))
  }
  dput(config, file = file.path(config$output_dir, "config_used.R"))
  writeLines(c(paste("Input MD5:", unname(tools::md5sum(config$input))), "Missingness audit cannot determine MCAR/MAR/MNAR or the correct adjustment set."), file.path(config$output_dir, "receipt.txt"))
  writeLines(capture.output(sessionInfo()), file.path(config$output_dir, "sessionInfo.txt"))
  invisible(list(missingness = missingness, results = results))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2L || args[1] != "--config") stop("Usage: Rscript 607_missingness_confounding.R --config config.json")
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Missing CLI dependency: jsonlite")
  run_missingness_confounding(jsonlite::fromJSON(args[2], simplifyVector = TRUE))
}
