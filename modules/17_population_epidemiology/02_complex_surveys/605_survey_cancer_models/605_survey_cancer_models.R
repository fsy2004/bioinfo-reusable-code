# Original adapter to survey/survival. No data download or package installation.
run_survey_cancer <- function(config) {
  for (pkg in c("survey", "survival")) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing dependency: ", pkg)
  }
  required <- c("input", "output_dir", "id", "psu", "strata", "weight", "outcome", "exposure", "weight_provenance")
  if (!all(required %in% names(config))) stop("Config requires: ", paste(required, collapse = ", "))
  if (!nzchar(config$weight_provenance)) stop("Document weight variable, cycles and pooling rule in weight_provenance")
  vars <- unique(unlist(config[c("id", "psu", "strata", "weight", "outcome", "exposure", "covariates", "domain", "modifier", "descriptives", "group")], use.names = FALSE))
  if (!is.null(config$mortality)) vars <- unique(c(vars, config$mortality$time, config$mortality$status))
  if (any(!grepl("^[A-Za-z.][A-Za-z0-9._]*$", vars))) stop("Use syntactic column names (letters, digits, underscore, dot)")
  dat <- read.csv(config$input, check.names = FALSE, stringsAsFactors = FALSE)
  if (anyDuplicated(names(dat)) || !all(vars %in% names(dat))) stop("Missing or duplicate input columns")
  if (anyNA(dat[[config$id]]) || anyDuplicated(dat[[config$id]])) stop("Input must contain one row per nonmissing participant ID")
  core <- c(config$psu, config$strata, config$weight)
  if (anyNA(dat[core])) stop("Missing survey-design variables; resolve before analysis")
  if (!is.numeric(dat[[config$weight]]) || any(!is.finite(dat[[config$weight]]) | dat[[config$weight]] <= 0)) stop("Weights must be finite and positive")
  binary <- function(x, label) {
    if (!all(is.na(x) | x %in% c(0, 1))) stop(label, " must be explicitly coded 0/1 with NA permitted")
  }
  binary(dat[[config$outcome]], config$outcome)
  # Explicitly chosen policy; do not silently select a lonely-PSU correction.
  lonely <- config$lonely_psu
  if (is.null(lonely) || !lonely %in% c("fail", "remove", "certainty", "adjust", "average")) stop("Specify lonely_psu policy")
  old <- options(survey.lonely.psu = lonely)
  on.exit(options(old), add = TRUE)
  one <- function(name) stats::reformulate(name)
  design <- survey::svydesign(ids = one(config$psu), strata = one(config$strata),
                             weights = one(config$weight), data = dat, nest = TRUE)
  full_n <- nrow(dat)
  if (!is.null(config$domain)) {
    binary(dat[[config$domain]], config$domain)
    keep <- !is.na(dat[[config$domain]]) & dat[[config$domain]] == 1
    design <- design[keep, ] # Domain subset is applied after constructing the full design.
  }
  if (!nrow(design$variables)) stop("Empty domain")
  covars <- unlist(config$covariates, use.names = FALSE)
  rhs <- unique(c(config$exposure, covars))
  form <- stats::reformulate(rhs, response = config$outcome)
  models <- list(logistic = survey::svyglm(form, design = design, family = quasibinomial(), na.action = na.omit))
  tests <- list()
  if (!is.null(config$modifier)) {
    term <- paste(config$exposure, config$modifier, sep = "*")
    f <- stats::reformulate(c(term, setdiff(covars, config$modifier)), response = config$outcome)
    models$interaction <- survey::svyglm(f, design = design, family = quasibinomial(), na.action = na.omit)
    tests$interaction <- survey::regTermTest(models$interaction, stats::reformulate(paste(config$exposure, config$modifier, sep = ":")))
  }
  if (!is.null(config$spline)) {
    knots <- as.numeric(unlist(config$spline$knots))
    boundary <- as.numeric(unlist(config$spline$boundary_knots))
    if (!is.numeric(dat[[config$exposure]]) || length(boundary) != 2L || length(knots) < 1L ||
        any(!is.finite(c(knots, boundary))) || anyDuplicated(knots) || is.unsorted(knots, strictly = TRUE) ||
        boundary[1] >= boundary[2] || any(knots <= boundary[1] | knots >= boundary[2])) stop("Supply numeric exposure and prespecified interior/boundary knots")
    term <- sprintf("splines::ns(%s, knots=c(%s), Boundary.knots=c(%s))", config$exposure, paste(knots, collapse = ","), paste(boundary, collapse = ","))
    f <- stats::reformulate(c(term, covars), response = config$outcome)
    models$natural_spline <- survey::svyglm(f, design = design, family = quasibinomial(), na.action = na.omit)
  }
  if (!is.null(config$mortality)) {
    m <- config$mortality
    if (is.null(m$endpoint) || !m$endpoint %in% c("all_cause_mortality", "cancer_mortality")) stop("mortality.endpoint must state all_cause_mortality or cancer_mortality")
    binary(dat[[m$status]], m$status)
    if (!is.numeric(dat[[m$time]]) || any(dat[[m$time]] < 0 | !is.finite(dat[[m$time]]), na.rm = TRUE)) stop("Mortality follow-up must be finite, nonnegative and have documented units")
    response <- sprintf("survival::Surv(%s,%s)", m$time, m$status)
    models$mortality_cox <- survey::svycoxph(stats::reformulate(rhs, response = response), design = design, na.action = na.omit)
  }
  descriptive_vars <- unique(c(config$outcome, unlist(config$descriptives, use.names = FALSE)))
  descriptives <- lapply(descriptive_vars, function(v) survey::svymean(one(v), design = design, na.rm = TRUE))
  names(descriptives) <- descriptive_vars
  grouped <- if (!is.null(config$group)) survey::svyby(one(config$outcome), one(config$group), design, survey::svymean, na.rm = TRUE, vartype = c("se", "ci")) else NULL
  coef_table <- function(model, name) {
    beta <- stats::coef(model)
    se <- sqrt(diag(stats::vcov(model)))
    ci <- suppressMessages(stats::confint(model))
    data.frame(model = name, term = names(beta), estimate = as.numeric(beta), std_error = as.numeric(se),
               lower = ci[, 1], upper = ci[, 2], exp_estimate = exp(beta), exp_lower = exp(ci[, 1]), exp_upper = exp(ci[, 2]),
               used_n = stats::nobs(model), row.names = NULL)
  }
  tables <- do.call(rbind, Map(coef_table, models, names(models)))
  if (dir.exists(config$output_dir) && length(list.files(config$output_dir, all.files = TRUE, no.. = TRUE))) stop("Output directory must be empty")
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(tables, file.path(config$output_dir, "model_coefficients.csv"), row.names = FALSE)
  if (!is.null(grouped)) write.csv(grouped, file.path(config$output_dir, "weighted_prevalence_by_group.csv"), row.names = FALSE)
  saveRDS(list(models = models, descriptives = descriptives, tests = tests), file.path(config$output_dir, "survey_models.rds"))
  writeLines(capture.output(print(descriptives)), file.path(config$output_dir, "weighted_descriptives.txt"))
  writeLines(capture.output(print(tests)), file.path(config$output_dir, "interaction_tests.txt"))
  dput(config, file = file.path(config$output_dir, "config_used.R"))
  writeLines(capture.output(sessionInfo()), file.path(config$output_dir, "sessionInfo.txt"))
  receipt <- c(paste("Full design rows:", full_n), paste("Domain rows:", nrow(design$variables)),
               paste("Input MD5:", unname(tools::md5sum(config$input))), paste("Weight provenance:", config$weight_provenance),
               "Model-specific missing rows are omitted explicitly; used_n is exported for each model.",
               "Logistic: prevalent cancer odds association. Mortality Cox: death hazard association, never incident-cancer risk.",
               "Natural-spline basis coefficients are not exposure-specific odds ratios; exp_estimate is for each basis coefficient only.")
  writeLines(receipt, file.path(config$output_dir, "receipt.txt"))
  invisible(list(models = models, coefficients = tables, design_n = full_n, domain_n = nrow(design$variables)))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2L || args[1] != "--config") stop("Usage: Rscript 605_survey_cancer_models.R --config config.json")
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Missing CLI dependency: jsonlite")
  run_survey_cancer(jsonlite::fromJSON(args[2], simplifyVector = TRUE))
}
