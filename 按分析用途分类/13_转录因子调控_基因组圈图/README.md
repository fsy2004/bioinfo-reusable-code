# 13 · Transcription factor regulation / genomic circos plots

TF/motif regulatory interpretation of candidate gene sets, regulon activity inference, and chromosomal circos display.

| Module | Purpose | Language | Output figures | Status |
|------|------|------|--------|:---:|
| [053 Chromosome circos](053_circlize染色体圈图/) | Circos of gene chromosomal positions | R | Circos plot | Ready |
| 047 RcisTarget motif-TF network | motif/TF enrichment + regulatory network | R | Network, Sankey | Requires cisTarget DB |
| 081 pySCENIC regulon | GRN + ctx + AUCell, TF activity | Python | UMAP, heatmap | Heavy environment |

> **053**: a coordinate table produces the circos plot (follows the [unified framework conventions](../_framework/CONVENTIONS.md)).
> **047**: requires the RcisTarget motif ranking database (GB scale); **081**: pySCENIC (Python + GRNBoost, heavy). Both retain the original scripts for reference.

## Recommended workflow
WGCNA/DEG/marker genes, RcisTarget motif, pySCENIC regulon, AUCell/decoupler activity, trajectory/spatial/perturbation interpretation.
