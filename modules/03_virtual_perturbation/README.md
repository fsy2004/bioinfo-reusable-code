# 03 · 虚拟扰动技术 — Virtual perturbation

本域共 15 个条目。完整字段见 [`../CATALOG.md`](../CATALOG.md)。


## 虚拟敲除与扰动模拟 — In-silico knockout

- [026_scrna_seurat_sctenifoldknk_ko.R](01_insilico_knockout/026_scrna_seurat_sctenifoldknk_ko.R) — Full QC/cluster/annotate + virtual KO pipeline
- [067_scperturb_etest.R](01_insilico_knockout/067_scperturb_etest.R) — scPerturb perturbation distance + E-test (light)
- [068_gears_combo_perturbation.py](01_insilico_knockout/068_gears_combo_perturbation.py) — GEARS single/combo perturbation prediction (GPU)
- [069_celloracle_grn_perturbation.py](01_insilico_knockout/069_celloracle_grn_perturbation.py) — CellOracle GRN virtual knockout (heavy)
- [085_squidiff_diffusion_perturbation.py](01_insilico_knockout/085_squidiff_diffusion_perturbation.py) — Squidiff/PerturbDiff diffusion perturbation (GPU)
- [494_genki_vgae_ko.py](01_insilico_knockout/494_genki_vgae_ko.py) — GenKI graph-VGAE virtual KO (KL ranking)
- [495_bulkvgk_sctenifoldknk_di.R](01_insilico_knockout/495_bulkvgk_sctenifoldknk_di.R) — Bulk co-expression virtual KO + differential influence
- [507_geneformer_insilico](01_insilico_knockout/507_geneformer_insilico) — Geneformer zero-shot embedding + in-silico deletion (baseline local)
- [561_regvelo_grn_velocity](01_insilico_knockout/561_regvelo_grn_velocity) — (用途待补)

## 基因调控网络推断 — GRN inference

- [047_rcistarget_tf_motif_network.R](02_grn_inference/047_rcistarget_tf_motif_network.R) — RcisTarget motif/TF enrichment + regulatory network
- [081_pyscenic_regulon_tf_activity.py](02_grn_inference/081_pyscenic_regulon_tf_activity.py) — pySCENIC GRN + ctx + AUCell wrapper
- [511_tf_convergence_depmap_jaspar](02_grn_inference/511_tf_convergence_depmap_jaspar) — Three-evidence convergence to core TFs

## 药物扰动与响应 — Drug perturbation

- [070_chemcpa_drug_perturbation.py](04_drug_perturbation/070_chemcpa_drug_perturbation.py) — chemCPA drug-perturbation expression prediction (GPU)
- [071_scdrug_response_prediction.py](04_drug_perturbation/071_scdrug_response_prediction.py) — scDrug single-cell drug response (heavy)
- [518_beyondcell_drug_response](04_drug_perturbation/518_beyondcell_drug_response) — beyondcell core re-impl: BCS + therapeutic clusters
