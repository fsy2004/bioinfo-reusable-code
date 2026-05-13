# scPerturb / scperturbR perturbation distance and E-test wrapper
# Example:
# Rscript 067_scPerturb_扰动数据Etest.R input_rds=seurat.rds group_col=perturbation control=control output_dir=results/scperturb

args <- commandArgs(trailingOnly = TRUE)
kv <- strsplit(args, "=", fixed = TRUE)
opts <- setNames(vapply(kv, function(x) if (length(x) > 1) x[2] else "", ""), vapply(kv, `[`, "", 1))

get_opt <- function(name, default = NULL) {
  value <- opts[[name]]
  if (is.null(value) || identical(value, "")) default else value
}

input_rds <- get_opt("input_rds")
group_col <- get_opt("group_col", "perturbation")
control <- get_opt("control", "control")
output_dir <- get_opt("output_dir", "results/scperturb")
assay <- get_opt("assay", "RNA")

if (is.null(input_rds)) stop("Please provide input_rds=path/to/seurat.rds")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required <- c("Seurat", "scperturbR")
missing_pkgs <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "),
       "\nInstall scperturbR from CRAN and Seurat before running this script.")
}

obj <- readRDS(input_rds)
if (!group_col %in% colnames(obj@meta.data)) {
  stop("group_col not found in Seurat metadata: ", group_col)
}

DefaultAssay(obj) <- assay
expr <- as.matrix(Seurat::GetAssayData(obj, assay = assay, slot = "data"))
meta <- obj@meta.data

saveRDS(list(expression = expr, metadata = meta), file.path(output_dir, "scperturb_input_expression_metadata.rds"))

ns <- asNamespace("scperturbR")
has_edist <- exists("edist", envir = ns, inherits = FALSE)
has_etest <- exists("etest", envir = ns, inherits = FALSE)

if (!has_edist || !has_etest) {
  writeLines(c(
    "scperturbR is installed, but edist/etest functions were not found.",
    "The Seurat expression and metadata bundle has been saved for downstream scPerturb/pertpy use.",
    "Output: scperturb_input_expression_metadata.rds"
  ), file.path(output_dir, "README_next_step.txt"))
  quit(save = "no", status = 0)
}

edata <- data.frame(t(expr), check.names = FALSE)
edata[[group_col]] <- meta[[group_col]]

edist_fun <- get("edist", envir = ns)
etest_fun <- get("etest", envir = ns)

edist_res <- edist_fun(edata, obs_key = group_col)
etest_res <- etest_fun(edata, obs_key = group_col, control = control)

write.csv(edist_res, file.path(output_dir, "scperturb_edistance.csv"), row.names = TRUE)
write.csv(etest_res, file.path(output_dir, "scperturb_etest_vs_control.csv"), row.names = FALSE)

message("Done. Results saved to: ", normalizePath(output_dir))

