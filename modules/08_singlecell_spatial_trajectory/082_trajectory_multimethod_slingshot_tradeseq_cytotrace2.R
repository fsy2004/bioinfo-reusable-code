#!/usr/bin/env Rscript

# Trajectory consensus wrapper for Slingshot, tradeSeq and CytoTRACE2.

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

args <- parse_args()
if (is.null(args$seurat_rds) || is.null(args$outdir)) {
  stop(paste(
    "Usage:",
    "Rscript 082_trajectory_multimethod_slingshot_tradeseq_cytotrace2.R",
    "--seurat_rds object.rds --outdir results/trajectory_consensus",
    "[--cluster_col seurat_clusters] [--reduction umap] [--start_cluster 0]",
    sep = " "
  ))
}

needed <- c("Seurat", "SingleCellExperiment", "slingshot")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) stop("Missing packages: ", paste(missing, collapse = ", "))

dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
cluster_col <- if (!is.null(args$cluster_col)) args$cluster_col else "seurat_clusters"
reduction <- if (!is.null(args$reduction)) args$reduction else "umap"

obj <- readRDS(args$seurat_rds)
if (!cluster_col %in% colnames(obj@meta.data)) stop("cluster_col not found in Seurat metadata.")
if (!reduction %in% names(obj@reductions)) stop("reduction not found in Seurat reductions.")

sce <- Seurat::as.SingleCellExperiment(obj)
SingleCellExperiment::reducedDim(sce, toupper(reduction)) <- Seurat::Embeddings(obj, reduction = reduction)
sce$cluster_for_trajectory <- as.character(obj@meta.data[[cluster_col]])

start.clus <- if (!is.null(args$start_cluster)) args$start_cluster else NULL
sce <- slingshot::slingshot(
  sce,
  clusterLabels = "cluster_for_trajectory",
  reducedDim = toupper(reduction),
  start.clus = start.clus
)

pt <- slingshot::slingPseudotime(sce)
write.csv(pt, file.path(args$outdir, "slingshot_pseudotime.csv"))
saveRDS(sce, file.path(args$outdir, "slingshot_sce.rds"))

if (requireNamespace("tradeSeq", quietly = TRUE)) {
  counts <- as.matrix(SingleCellExperiment::counts(sce))
  set.seed(1)
  fit <- tradeSeq::fitGAM(counts = counts, sds = slingshot::SlingshotDataSet(sce), nknots = 6)
  assoc <- tradeSeq::associationTest(fit)
  write.csv(assoc, file.path(args$outdir, "tradeseq_association_test.csv"))
  saveRDS(fit, file.path(args$outdir, "tradeseq_fit.rds"))
}

if (requireNamespace("CytoTRACE2", quietly = TRUE)) {
  expr <- as.matrix(SingleCellExperiment::counts(sce))
  cyt <- CytoTRACE2::cytotrace2(expr)
  saveRDS(cyt, file.path(args$outdir, "cytotrace2_result.rds"))
}

summary <- data.frame(
  cells = ncol(sce),
  genes = nrow(sce),
  clusters = length(unique(sce$cluster_for_trajectory)),
  reduction = reduction
)
write.csv(summary, file.path(args$outdir, "trajectory_consensus_summary.csv"), row.names = FALSE)
