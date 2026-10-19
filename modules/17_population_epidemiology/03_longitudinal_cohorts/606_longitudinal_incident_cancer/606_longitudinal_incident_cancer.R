# Original adapters; inputs must already encode a defensible cohort and time origin.
run_incident_cancer <- function(config) {
  if (!requireNamespace("survival", quietly = TRUE)) stop("Missing dependency: survival")
  required <- c("input", "output_dir", "mode", "id", "exposure", "baseline_cancer", "time_origin", "time_unit")
  if (!all(required %in% names(config))) stop("Config requires: ", paste(required, collapse = ", "))
  modes <- c("cause_specific_cox", "interval_aft", "discrete_time")
  if (!config$mode %in% modes) stop("mode must be one of: ", paste(modes, collapse = ", "))
  if (!nzchar(config$time_origin) || !nzchar(config$time_unit)) stop("Document time origin and unit")
  covars <- unlist(config$covariates, use.names = FALSE)
  vars <- unique(c(config$id, config$exposure, config$baseline_cancer, covars))
  if (any(!grepl("^[A-Za-z.][A-Za-z0-9._]*$", vars))) stop("Use syntactic column names")
  dat <- read.csv(config$input, stringsAsFactors = FALSE, check.names = FALSE)
  timevars <- switch(config$mode, cause_specific_cox = c("entry_time", "exit_time", "status"),
                     interval_aft = c("left_time", "right_time", "status"),
                     discrete_time = c("period_start", "period_end", "period", "status"))
  if (anyDuplicated(names(dat)) || !all(c(vars, timevars) %in% names(dat))) stop("Missing or duplicate input columns")
  if (anyNA(dat[[config$id]])) stop("Missing participant ID")
  if (!all(is.na(dat[[config$baseline_cancer]]) | dat[[config$baseline_cancer]] %in% c(0, 1))) stop("baseline_cancer must be coded 0/1/NA")
  baseline_by_person <- split(dat[[config$baseline_cancer]], dat[[config$id]])
  if (any(vapply(baseline_by_person, function(x) length(unique(x)) != 1L, logical(1)))) stop("Baseline cancer status inconsistent within participant")
  flow <- c(input_rows = nrow(dat), input_people = length(unique(dat[[config$id]])),
            excluded_baseline_cancer_people = sum(vapply(baseline_by_person, function(x) !is.na(x[1]) && x[1] == 1, logical(1))),
            excluded_unknown_baseline_people = sum(vapply(baseline_by_person, function(x) is.na(x[1]), logical(1))))
  dat <- dat[!is.na(dat[[config$baseline_cancer]]) & dat[[config$baseline_cancer]] == 0, , drop = FALSE]
  if (!nrow(dat)) stop("No participants with documented cancer-free baseline")
  if (anyNA(dat[c(config$exposure, covars, timevars)])) {
    # right_time=NA is the explicit right-censoring representation in interval mode.
    check <- c(config$exposure, covars, setdiff(timevars, if (config$mode == "interval_aft") "right_time" else character()))
    if (anyNA(dat[check])) stop("Missing model/time inputs; resolve before analysis rather than silently dropping rows")
  }
  if (!all(dat$status %in% c(0, 1, 2))) stop("status must be 0=censored, 1=incident cancer, 2=competing death")
  if (!any(dat$status == 1)) stop("No incident cancer events")
  rhs <- unique(c(config$exposure, covars))
  ph_test <- NULL
  numeric_finite <- function(columns) {
    if (any(!vapply(dat[columns], is.numeric, logical(1))) || any(!is.finite(as.matrix(dat[columns])))) stop("Time inputs must be finite numeric")
  }
  if (config$mode == "cause_specific_cox") {
    if (anyDuplicated(dat[[config$id]])) stop("Cause-specific mode requires one first-event record per person")
    numeric_finite(c("entry_time", "exit_time"))
    if (any(dat$entry_time < 0 | dat$exit_time <= dat$entry_time)) stop("Require 0 <= entry_time < exit_time")
    f <- stats::reformulate(rhs, response = "survival::Surv(entry_time, exit_time, status == 1)")
    model <- survival::coxph(f, data = dat, ties = "efron", x = TRUE, model = TRUE)
    ph_test <- survival::cox.zph(model)
    interpretation <- "Cause-specific cancer HR among people still alive and cancer free. Competing deaths exit the risk set at their observed times. Does not estimate cumulative incidence or a Fine-Gray HR."
  } else if (config$mode == "interval_aft") {
    if (anyDuplicated(dat[[config$id]])) stop("Interval mode requires one cancer interval/censoring record per person")
    numeric_finite("left_time")
    if (!is.numeric(dat$right_time) || any(!is.finite(dat$right_time[!is.na(dat$right_time)]))) stop("right_time must be numeric or NA for right censoring")
    if (any(dat$status == 2)) stop("Interval AFT adapter does not handle competing deaths. Preserve these records and use an explicitly justified competing-risk analysis.")
    if (any(dat$left_time < 0) || any(dat$status == 0 & !is.na(dat$right_time)) ||
        any(dat$status == 1 & (is.na(dat$right_time) | dat$right_time < dat$left_time | dat$right_time <= 0))) stop("Inconsistent interval/status: censored right_time=NA; event 0<=left_time<=right_time with positive right_time")
    dist <- config$aft_distribution
    if (is.null(dist) || !dist %in% c("weibull", "lognormal", "loglogistic")) stop("Prespecify aft_distribution: weibull/lognormal/loglogistic")
    f <- stats::reformulate(rhs, response = "survival::Surv(left_time, right_time, type='interval2')")
    model <- survival::survreg(f, data = dat, dist = dist)
    interpretation <- "Parametric interval-censored AFT time ratio, conditional on model and censoring assumptions. No diagnosis midpoint imputation; no left truncation or competing-death extension in this adapter."
  } else {
    numeric_finite(c("period_start", "period_end"))
    if (any(dat$period_start < 0 | dat$period_end <= dat$period_start)) stop("Require valid positive-width person-period intervals")
    period_groups <- split(seq_len(nrow(dat)), dat$period)
    if (any(vapply(period_groups, function(indexes) length(unique(dat$period_start[indexes])) != 1L ||
                   length(unique(dat$period_end[indexes])) != 1L, logical(1)))) {
      stop("Each period must have identical time boundaries across participants. Irregular interview intervals require a different time model.")
    }
    groups <- split(seq_len(nrow(dat)), dat[[config$id]])
    for (indexes in groups) {
      z <- dat[indexes, , drop = FALSE]
      z <- z[order(z$period_start), , drop = FALSE]
      if (nrow(z) > 1L && any(z$period_start[-1] != head(z$period_end, -1))) stop("Intervals must be contiguous and nonoverlapping within person")
      if (anyDuplicated(z$period) || sum(z$status != 0) > 1L || (any(z$status != 0) && tail(z$status, 1) == 0)) stop("At most one cancer/death event, in final period; period labels unique within person")
      if (any(vapply(z[c(config$exposure, covars)], function(x) length(unique(x)) != 1L, logical(1)))) stop("This adapter uses baseline exposure/covariates only; time-varying effects require a separate design")
    }
    # A saturated shared period factor accommodates survey-wave baseline hazards.
    dat$period <- factor(dat$period)
    f <- stats::reformulate(c(rhs, "period"), response = "I(status == 1)")
    model <- stats::glm(f, data = dat, family = binomial(link = "cloglog"))
    # Cluster sandwich for repeated person-period rows, not independent row-level SEs.
    score <- stats::model.matrix(model) * as.numeric((model$y - stats::fitted(model)) * model$family$mu.eta(stats::predict(model, type = "link")) / model$family$variance(stats::fitted(model)))
    if (anyNA(stats::coef(model))) stop("Rank-deficient discrete-time design; revise prespecified covariates/periods")
    clusters <- rowsum(score, group = dat[[config$id]], reorder = FALSE)
    bread <- summary(model)$cov.unscaled
    robust_vcov <- bread %*% crossprod(clusters) %*% bread
    interpretation <- "Cause-specific discrete-time complementary-log-log association with participant-cluster sandwich SE. Each period has identical boundaries across participants; period factors absorb duration differences between periods. Final competing-death period has no cancer event. Period risk does not recover exact event dates or cumulative incidence."
  }
  beta <- stats::coef(model)
  vc <- if (config$mode == "discrete_time") robust_vcov else stats::vcov(model)
  # survreg includes an additional log(scale) covariance parameter.
  se <- sqrt(diag(vc))[seq_along(beta)]
  coefficients <- data.frame(term = names(beta), estimate = as.numeric(beta), std_error = as.numeric(se),
                             exp_estimate = exp(beta), exp_lower = exp(beta - 1.96 * se), exp_upper = exp(beta + 1.96 * se), row.names = NULL)
  flow <- c(flow, included_rows = nrow(dat), included_people = length(unique(dat[[config$id]])),
            incident_events = sum(dat$status == 1), competing_deaths = sum(dat$status == 2))
  if (dir.exists(config$output_dir) && length(list.files(config$output_dir, all.files = TRUE, no.. = TRUE))) stop("Output directory must be empty")
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(coefficients, file.path(config$output_dir, "coefficients.csv"), row.names = FALSE)
  write.csv(data.frame(stage = names(flow), n = as.integer(flow)), file.path(config$output_dir, "cohort_flow.csv"), row.names = FALSE)
  saveRDS(list(model = model, ph_test = ph_test, vcov = vc), file.path(config$output_dir, "model.rds"))
  if (!is.null(ph_test)) writeLines(capture.output(ph_test), file.path(config$output_dir, "proportional_hazards_test.txt"))
  dput(config, file = file.path(config$output_dir, "config_used.R"))
  writeLines(c(interpretation, paste("Input MD5:", unname(tools::md5sum(config$input))),
               "No survey weights, attrition weights, MI, target-trial emulation or causal identification are implied by this adapter."), file.path(config$output_dir, "receipt.txt"))
  writeLines(capture.output(sessionInfo()), file.path(config$output_dir, "sessionInfo.txt"))
  invisible(list(model = model, coefficients = coefficients, flow = flow))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2L || args[1] != "--config") stop("Usage: Rscript 606_longitudinal_incident_cancer.R --config config.json")
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Missing CLI dependency: jsonlite")
  run_incident_cancer(jsonlite::fromJSON(args[2], simplifyVector = TRUE))
}
