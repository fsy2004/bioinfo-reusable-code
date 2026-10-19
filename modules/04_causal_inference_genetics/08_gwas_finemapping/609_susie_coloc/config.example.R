# Illustrative paths and metadata only; replace from actual source documentation.
cfg <- list(
  seed = 20261019L, out_dir = "results", L = 10L, coverage = 0.95,
  max_iter = 1000L, allow_case_control_rss_approximation = FALSE,
  prior_grid = data.frame(p1 = rep(1e-4, 3), p2 = rep(1e-4, 3), p12 = c(1e-6, 5e-6, 1e-5)),
  trait1 = list(trait_id = "REPLACE_EQTL_ID", phenotype = "REPLACE gene expression",
    genome_build = "GRCh38", ancestry = "REPLACE", beta_unit = "REPLACE documented expression unit",
    type = "quant", N = 10000, sdY = 1, source = "REPLACE accession + release",
    region = "REPLACE chr:start-end", dense_region_confirmed = TRUE,
    summary = "trait1.tsv", ld_matrix = "ld1.txt", ld_variants = "ld1_variants.tsv",
    ld_genome_build = "GRCh38", ld_ancestry = "REPLACE", ld_metric = "signed_r",
    ld_coded_allele = "ea", ld_source = "REPLACE reference panel + release"),
  trait2 = list(trait_id = "REPLACE_OUTCOME_ID", phenotype = "REPLACE quantitative phenotype",
    genome_build = "GRCh38", ancestry = "REPLACE", beta_unit = "REPLACE documented unit",
    type = "quant", N = 50000, sdY = 1, source = "REPLACE accession + release",
    region = "REPLACE chr:start-end", dense_region_confirmed = TRUE,
    summary = "trait2.tsv", ld_matrix = "ld2.txt", ld_variants = "ld2_variants.tsv",
    ld_genome_build = "GRCh38", ld_ancestry = "REPLACE", ld_metric = "signed_r",
    ld_coded_allele = "ea", ld_source = "REPLACE reference panel + release")
)
# For a binary cancer outcome use type="cc", beta_unit="log odds", N=cases+controls,
# s=cases/N (instead of sdY), and explicitly set allow_case_control_rss_approximation=TRUE.
# The numeric N/sdY placeholders above are examples, never defaults for real datasets.
