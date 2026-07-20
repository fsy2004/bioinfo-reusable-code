#!/usr/bin/env Rscript

# NicheNet ligand-target wrapper.
# Typical use: connect CellChat/COMMOT sender-receiver pairs to receiver DE genes.

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

read_gene_vector <- function(path) {
  x <- read.table(path, header = TRUE, sep = "", stringsAsFactors = FALSE, check.names = FALSE)
  gene_col <- intersect(c("gene", "genes", "symbol", "Gene", "GENE"), colnames(x))[1]
  if (is.na(gene_col)) gene_col <- colnames(x)[1]
  unique(na.omit(as.character(x[[gene_col]])))
}

read_any <- function(path) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) return(readRDS(path))
  read.table(path, header = TRUE, sep = "", stringsAsFactors = FALSE, check.names = FALSE)
}

args <- parse_args()
if (is.null(args$receiver_genes) || is.null(args$background_genes) ||
    is.null(args$ligand_target_matrix) || is.null(args$lr_network) ||
    is.null(args$outdir)) {
  stop(paste(
    "Usage:",
    "Rscript 077_nichenet_ligand_target.R",
    "--receiver_genes receiver_de_genes.tsv",
    "--background_genes expressed_genes.tsv",
    "--ligand_target_matrix ligand_target_matrix.rds",
    "--lr_network lr_network.rds",
    "--outdir results/nichenet",
    "[--expressed_ligands expressed_ligands.tsv]",
    "[--top_n 30]",
    sep = " "
  ))
}

if (!requireNamespace("nichenetr", quietly = TRUE)) {
  stop("Package 'nichenetr' is required. Install it before running this wrapper.")
}

dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
receiver_genes <- read_gene_vector(args$receiver_genes)
background_genes <- read_gene_vector(args$background_genes)
ligand_target_matrix <- readRDS(args$ligand_target_matrix)
lr_network <- read_any(args$lr_network)
expressed_ligands <- if (!is.null(args$expressed_ligands)) read_gene_vector(args$expressed_ligands) else unique(lr_network$from)
top_n <- if (!is.null(args$top_n)) as.integer(args$top_n) else 30

candidate_ligands <- intersect(expressed_ligands, rownames(ligand_target_matrix))
if (length(candidate_ligands) == 0) {
  stop("No candidate ligands overlap with ligand_target_matrix rownames.")
}

ligand_activities <- nichenetr::predict_ligand_activities(
  geneset = receiver_genes,
  background_expressed_genes = background_genes,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = candidate_ligands
)
ligand_activities <- ligand_activities[order(ligand_activities$pearson, decreasing = TRUE), , drop = FALSE]

top_ligands <- head(ligand_activities$test_ligand, top_n)
ligand_target_links <- nichenetr::get_weighted_ligand_target_links(
  ligands = top_ligands,
  geneset = receiver_genes,
  ligand_target_matrix = ligand_target_matrix,
  n = 250
)

write.csv(ligand_activities, file.path(args$outdir, "nichenet_ligand_activities.csv"), row.names = FALSE)
write.csv(ligand_target_links, file.path(args$outdir, "nichenet_ligand_target_links.csv"), row.names = FALSE)
writeLines(top_ligands, file.path(args$outdir, "nichenet_top_ligands.txt"))

summary <- data.frame(
  receiver_genes = length(receiver_genes),
  background_genes = length(background_genes),
  candidate_ligands = length(candidate_ligands),
  top_ligands = length(top_ligands)
)
write.csv(summary, file.path(args$outdir, "nichenet_run_summary.csv"), row.names = FALSE)
