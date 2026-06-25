# 16 · Spatial communication and cell fate

Methods for spatial cell-cell communication, single-cell to spatial mapping, cell fate inference, spatial niches, and TF/pathway activity scoring.

## Status

These methods rely on the Python spatial-omics stack (CellRank/COMMOT/Tangram/cell2location/Squidpy/decoupler) or the large NicheNet prior network. Example figures are not rendered locally; the original scripts are kept for reference, with dependencies and reproducible commands in each script header. For figure conventions see [unified framework](../_framework/CONVENTIONS.md).

## Scripts

| Script | Purpose |
|---|---|
| `072_cellrank_fate_drivers.py` | Runs CellRank on a scVelo velocity graph, outputting terminal states, fate probabilities, and driver genes. |
| `073_commot_spatial_communication.py` | Computes spatial sender/receiver communication scores from spatial coordinates and a ligand-receptor database. |
| `074_tangram_sc_to_spatial.py` | Maps scRNA-seq cells onto spatial transcriptomics spots. |
| `076_decoupler_tf_pathway_activity.py` | Infers TF/pathway activity using CollecTRI or PROGENy. |
| `077_nichenet_ligand_target.R` | Infers ligand activity and ligand-target chains from receiver DE genes and an LR/target prior. |
| `080_cell2location_squidpy_niche.py` | Integrates a cell abundance table and runs Squidpy spatial neighborhood enrichment. |

## Input

- AnnData h5ad or Seurat conversion output.
- Spatial coordinates and tissue section information.
- Ligand-receptor database, receiver differential genes, background expressed genes.
- Cell type annotations or cell abundance matrix.

## Outputs

- CellRank fate probabilities and driver genes.
- COMMOT spatial communication network.
- Tangram single-cell to spatial mapping results.
- NicheNet ligand activity and ligand-target links.
- Squidpy neighborhood enrichment.
