#!/usr/bin/env Rscript

# FAERS disproportionality wrapper.
# Input either raw drug-event-case rows or precomputed n11/n10/n01/n00 counts.

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list()
  i <- 1
  while (i <= length(args)) {
    key <- sub("^--", "", args[[i]])
    val <- if (i + 1 <= length(args) && !grepl("^--", args[[i + 1]])) args[[i + 1]] else TRUE
    out[[key]] <- val
    i <- i + if (isTRUE(val)) 1 else 2
  }
  out
}

read_table_auto <- function(path) {
  sep <- if (grepl("\\.csv$", path, ignore.case = TRUE)) "," else "\t"
  read.table(path, header = TRUE, sep = sep, quote = "\"", stringsAsFactors = FALSE, check.names = FALSE)
}

build_counts <- function(x, drug_col, event_col, case_col) {
  x[[drug_col]] <- trimws(as.character(x[[drug_col]]))
  x[[event_col]] <- trimws(as.character(x[[event_col]]))
  x[[case_col]] <- as.character(x[[case_col]])
  x <- unique(x[, c(case_col, drug_col, event_col)])
  all_cases <- unique(x[[case_col]])
  drugs <- sort(unique(x[[drug_col]]))
  events <- sort(unique(x[[event_col]]))
  rows <- vector("list", length(drugs) * length(events))
  k <- 1
  for (d in drugs) {
    cases_d <- unique(x[[case_col]][x[[drug_col]] == d])
    for (e in events) {
      cases_e <- unique(x[[case_col]][x[[event_col]] == e])
      n11 <- length(intersect(cases_d, cases_e))
      n10 <- length(setdiff(cases_d, cases_e))
      n01 <- length(setdiff(cases_e, cases_d))
      n00 <- length(setdiff(all_cases, union(cases_d, cases_e)))
      rows[[k]] <- data.frame(drug = d, event = e, n11 = n11, n10 = n10, n01 = n01, n00 = n00)
      k <- k + 1
    }
  }
  do.call(rbind, rows)
}

calc_signals <- function(x) {
  for (nm in c("n11", "n10", "n01", "n00")) x[[nm]] <- as.numeric(x[[nm]])
  a <- x$n11 + 0.5
  b <- x$n10 + 0.5
  c <- x$n01 + 0.5
  d <- x$n00 + 0.5
  n <- a + b + c + d
  x$ROR <- (a * d) / (b * c)
  x$ROR_log <- log(x$ROR)
  x$ROR_se <- sqrt(1 / a + 1 / b + 1 / c + 1 / d)
  x$ROR025 <- exp(x$ROR_log - 1.96 * x$ROR_se)
  x$ROR975 <- exp(x$ROR_log + 1.96 * x$ROR_se)
  x$PRR <- (a / (a + b)) / (c / (c + d))
  expected <- ((a + b) * (a + c)) / n
  x$IC <- log2(a / expected)
  x$IC025 <- x$IC - 1.96 * sqrt(1 / pmax(a, 1))
  x$chi_square <- n * (a * d - b * c)^2 / ((a + b) * (c + d) * (a + c) * (b + d))
  x$signal_ROR <- x$n11 >= 3 & x$ROR025 > 1
  x$signal_PRR <- x$n11 >= 3 & x$PRR >= 2 & x$chi_square >= 4
  x$signal_BCPNN <- x$n11 >= 3 & x$IC025 > 0
  x$EBGM_proxy <- exp(x$IC * log(2))
  x$signal_EBGM_proxy <- x$n11 >= 3 & x$EBGM_proxy >= 2
  x
}

args <- parse_args()
if (is.null(args$input) || is.null(args$outdir)) {
  stop(paste(
    "Usage:",
    "Rscript 078_FAERS_ROR_PRR_BCPNN_EBGM药物警戒.R",
    "--input faers_drug_event.tsv --outdir results/faers",
    "[--drug_col drugname --event_col pt --case_col primaryid]",
    sep = " "
  ))
}

dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
x <- read_table_auto(args$input)

if (all(c("n11", "n10", "n01", "n00") %in% colnames(x))) {
  counts <- x
  if (!"drug" %in% colnames(counts)) counts$drug <- NA_character_
  if (!"event" %in% colnames(counts)) counts$event <- NA_character_
} else {
  drug_col <- if (!is.null(args$drug_col)) args$drug_col else "drug"
  event_col <- if (!is.null(args$event_col)) args$event_col else "event"
  case_col <- if (!is.null(args$case_col)) args$case_col else "case_id"
  missing <- setdiff(c(drug_col, event_col, case_col), colnames(x))
  if (length(missing) > 0) stop("Missing input columns: ", paste(missing, collapse = ", "))
  counts <- build_counts(x, drug_col, event_col, case_col)
}

signals <- calc_signals(counts)
signals <- signals[order(signals$signal_ROR | signals$signal_PRR | signals$signal_BCPNN, signals$ROR, decreasing = TRUE), ]
write.csv(signals, file.path(args$outdir, "faers_disproportionality_signals.csv"), row.names = FALSE)
write.csv(subset(signals, signal_ROR | signal_PRR | signal_BCPNN | signal_EBGM_proxy),
          file.path(args$outdir, "faers_positive_signals.csv"), row.names = FALSE)

summary <- data.frame(
  pairs = nrow(signals),
  ror_signals = sum(signals$signal_ROR, na.rm = TRUE),
  prr_signals = sum(signals$signal_PRR, na.rm = TRUE),
  bcpnn_signals = sum(signals$signal_BCPNN, na.rm = TRUE),
  ebgm_proxy_signals = sum(signals$signal_EBGM_proxy, na.rm = TRUE)
)
write.csv(summary, file.path(args$outdir, "faers_signal_summary.csv"), row.names = FALSE)
