# Artificial interface fixtures only. No scientific estimates are interpreted.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1])
repo <- normalizePath(file.path(dirname(script), ".."))
root <- file.path(repo, "modules", "17_population_epidemiology")
base <- tempfile("observational_interface_")
if (dir.exists(base)) stop("Validation output already exists; choose a new validation root")
dir.create(base, recursive = TRUE)
source(file.path(root, "02_complex_surveys", "605_survey_cancer_models", "605_survey_cancer_models.R"))
source(file.path(root, "03_longitudinal_cohorts", "606_longitudinal_incident_cancer", "606_longitudinal_incident_cancer.R"))
source(file.path(root, "06_confounding_missingness", "607_missingness_confounding", "607_missingness_confounding.R"))
passed <- character()
record <- function(label) { passed <<- c(passed, label); cat("PASS", label, "\n") }
rejects <- function(code) inherits(tryCatch({force(code); NULL}, error = identity), "error")
set.seed(604607)
n <- 480L
d <- data.frame(person_id = seq_len(n), strata = rep(1:12, each = 40), psu = rep(rep(1:2, each = 20), 12),
                analysis_weight = runif(n, 0.5, 3), eligible_domain = rep(c(1,1,1,0), length.out = n),
                age = runif(n, 20, 80), sex = sample(c("F", "M"), n, TRUE), smoking = rbinom(n, 1, 0.25))
d$exposure <- rbinom(n, 1, plogis(-1 + 0.01 * d$age + 0.3 * d$smoking))
d$prevalent_cancer <- rbinom(n, 1, plogis(-3 + 0.02 * d$age + 0.2 * d$exposure))
d$followup <- rexp(n, 0.05) + 0.1; d$death <- rbinom(n, 1, 0.25)
input <- file.path(base, "artificial_people.csv"); write.csv(d, input, row.names = FALSE)
survey_config <- list(input = input, output_dir = file.path(base, "survey"), id = "person_id", psu = "psu", strata = "strata",
                      weight = "analysis_weight", weight_provenance = "Artificial weights for software interface check only", lonely_psu = "fail",
                      outcome = "prevalent_cancer", exposure = "age", covariates = c("sex", "smoking"), domain = "eligible_domain",
                      modifier = "sex", descriptives = "age", group = "sex", spline = list(knots = c(40,60), boundary_knots = c(20,80)),
                      mortality = list(time = "followup", status = "death", endpoint = "all_cause_mortality"))
s <- run_survey_cancer(survey_config)
stopifnot(length(s$models) == 4L, s$design_n == n, s$domain_n == sum(d$eligible_domain == 1), all(s$coefficients$used_n <= s$domain_n))
record("605 full design/domain, descriptives, logistic, interaction, spline and mortality Cox")
bad <- survey_config; bad$weight_provenance <- ""
stopifnot(rejects(run_survey_cancer(bad)))
record("605 missing weight provenance is rejected")
cohort <- d
cohort$baseline_cancer <- rep(c(0, 0, 0, 0, 1, NA), length.out = n)
cohort$entry_time <- 0; cohort$exit_time <- rexp(n, .15) + .1
cohort$status <- sample(c(0,1,2), n, TRUE, c(.5,.3,.2))
cohort_file <- file.path(base,"artificial_cohort.csv"); write.csv(cohort, cohort_file, row.names = FALSE)
cconfig <- list(input = cohort_file, output_dir = file.path(base,"cause_cox"), mode = "cause_specific_cox", id = "person_id",
                exposure = "exposure", covariates = c("age","sex"), baseline_cancer = "baseline_cancer", time_origin = "Artificial baseline", time_unit = "years")
cfit <- run_incident_cancer(cconfig)
stopifnot(cfit$flow["included_people"] == sum(cohort$baseline_cancer == 0, na.rm = TRUE), cfit$flow["competing_deaths"] > 0)
record("606 cause-specific Cox with baseline exclusions and retained competing deaths")
interval <- cohort; interval$status[interval$status == 2] <- 0
interval$left_time <- floor(interval$exit_time)
interval$right_time <- ifelse(interval$status == 1, ceiling(interval$exit_time), NA_real_)
interval$left_time[interval$status == 0] <- interval$exit_time[interval$status == 0]
interval_file <- file.path(base,"artificial_intervals.csv"); write.csv(interval, interval_file, row.names = FALSE)
aconfig <- cconfig; aconfig$input <- interval_file; aconfig$output_dir <- file.path(base,"aft"); aconfig$mode <- "interval_aft"; aconfig$aft_distribution <- "weibull"
afit <- run_incident_cancer(aconfig); stopifnot(inherits(afit$model, "survreg"))
record("606 interval-censored AFT with explicit distribution")
periods <- do.call(rbind, lapply(seq_len(n), function(i) {
  last <- sample(2:4, 1); status <- sample(c(0,1,2), 1, prob = c(.5,.3,.2))
  z <- cohort[rep(i,last), , drop = FALSE]
  z$period <- seq_len(last); z$period_start <- z$period - 1; z$period_end <- z$period
  z$status <- c(rep(0,last-1),status); z
}))
period_file <- file.path(base,"artificial_periods.csv"); write.csv(periods, period_file, row.names = FALSE)
pconfig <- cconfig; pconfig$input <- period_file; pconfig$output_dir <- file.path(base,"discrete"); pconfig$mode <- "discrete_time"
pfit <- run_incident_cancer(pconfig); stopifnot(inherits(pfit$model,"glm"), all(is.finite(pfit$coefficients$std_error)))
record("606 discrete-time model and participant-cluster sandwich")
irregular <- periods; k <- which(!is.na(irregular$baseline_cancer) & irregular$baseline_cancer == 0 & irregular$period == 1)[1]
irregular$period_end[k] <- irregular$period_end[k] + 0.1
irregular_file <- file.path(base,"artificial_irregular_periods.csv"); write.csv(irregular, irregular_file, row.names = FALSE)
bad <- pconfig; bad$input <- irregular_file
stopifnot(rejects(run_incident_cancer(bad)))
record("606 unequal person-specific period boundaries are rejected")
bad <- aconfig; bad$input <- cohort_file
stopifnot(rejects(run_incident_cancer(bad)))
record("606 missing interval schema is rejected")
audit_config <- list(input = input, output_dir = file.path(base,"missing_audit"), mode = "audit", id = "person_id")
audit <- run_missingness_confounding(audit_config); stopifnot(all(audit$missingness$missing_n == 0))
record("607 missingness audit")
psconfig <- list(input = input, output_dir = file.path(base,"ps"), mode = "ps_diagnostics", id = "person_id", exposure = "exposure",
                 confounders = c("age","sex","smoking"), factor_variables = "sex", dag_assumptions = "Artificial pre-exposure covariates; no scientific claim",
                 estimand = "ATE", positivity_bounds = c(.01,.99))
psfit <- run_missingness_confounding(psconfig)
stopifnot(nrow(psfit$results$record_weights) == n, all(psfit$results$diagnostics$effective_sample_size <= psfit$results$diagnostics$n + 1e-8))
record("607 stabilized ATE propensity weights, fixed-denominator balance, overlap and ESS")
bad <- psconfig; bad$survey_weight <- "analysis_weight"
stopifnot(rejects(run_missingness_confounding(bad)))
record("607 survey/PS weight mixing is rejected")
mi <- d; mi$age[sample.int(n,40)] <- NA; mi$smoking[sample.int(n,30)] <- NA
mi_file <- file.path(base,"artificial_missing.csv"); write.csv(mi, mi_file, row.names = FALSE)
mconfig <- psconfig; mconfig$input <- mi_file; mconfig$output_dir <- file.path(base,"mi"); mconfig$mode <- "mi_glm"
mconfig$outcome <- "prevalent_cancer"; mconfig$family <- "binomial"; mconfig$factor_variables <- c("sex","smoking")
mconfig$imputation_variables <- c("exposure","age","sex","smoking","prevalent_cancer")
mconfig$imputation_methods <- list(age = "pmm", smoking = "logreg"); mconfig$m <- 3; mconfig$maxit <- 3; mconfig$seed <- 604607
mfit <- run_missingness_confounding(mconfig)
stopifnot(inherits(mfit$results$pooled,"mipo"), length(mfit$results$fits) == 3L)
record("607 explicit MICE methods, per-imputation GLM and Rubin pooling")
writeLines(c("Artificial software-interface checks only; not scientific validation.", passed), file.path(base,"validation_receipt.txt"))
