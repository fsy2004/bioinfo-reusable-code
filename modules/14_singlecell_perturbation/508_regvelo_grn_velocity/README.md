# 508 · RegVelo — gene-regulatory-informed RNA velocity & regulon perturbation

Uses **RegVelo** (Wang et al., *Cell* 2026) to couple **splicing dynamics with a gene-regulatory
network**, so velocity is driven by regulation rather than by expression kinetics alone. It also
supports **in-silico regulon perturbation** through the CellRank framework. Ships with a
**runnable baseline** (plain scVelo velocity + CellRank fate probabilities, no GRN) so the
GRN-informed model is always measured against a simple comparator, never reported alone.

| | |
|---|---|
| Language / deps | Python ≥3.10 · baseline: `scanpy` `scvelo` `cellrank`; RegVelo: `regvelo` (`scvi-tools<1.2.1` `torchode`) · **GPU recommended** |
| Purpose | GRN-informed RNA velocity, latent time, velocity uncertainty, regulon perturbation |
| Input | spliced/unspliced AnnData (`adata.layers['spliced'|'unspliced']`) + a TF→target prior (skeleton) |
| Output | `results/` velocity + fate + perturbation ranking; preview in `assets/` |
| Runtime | baseline CPU minutes · RegVelo training needs GPU |

## Method

**Runnable baseline (always, CPU):**
1. Standard preprocessing → `scvelo` moments → dynamical/stochastic velocity.
2. CellRank kernel from the velocity field → terminal states → fate probabilities.
3. This is the "no-GRN" floor: any RegVelo claim must beat or add to it.

**RegVelo path (`--run-regvelo`, GPU):**
Real exported API (verified from the package): `regvelo.REGVELOVI`, `VELOVAE`,
`ModelComparison`, plus `pp` / `tl` / `pl` / `mt` submodules. The scvi-tools pattern applies
(`setup_anndata` → construct → `train`). **Check the official tutorial for exact signatures
before a production run** — they are not pinned here.

- Docs / tutorials: https://regvelo.readthedocs.io
- Repo: https://github.com/theislab/regvelo

## ⏭️ Needs install + GPU (RegVelo path)

```bash
pip install regvelo          # Python >=3.10; pulls scvi-tools<1.2.1, torchode, cellrank
```

The script guards on import and GPU availability and otherwise runs the baseline only,
printing why. Nothing silently degrades.

## When to use this vs 069 / 507

| Module | Engine | Perturbation logic |
|---|---|---|
| 069 CellOracle | GRN + signal propagation on a vector field | knock a TF to zero, propagate, score the shift |
| 507 Geneformer | foundation-model embedding | delete a gene, measure embedding shift |
| **508 RegVelo** | **GRN coupled to splicing dynamics** | **regulon perturbation via CellRank fates** |

They are not interchangeable. RegVelo needs spliced/unspliced counts; CellOracle and
Geneformer do not. Using two engines that share assumptions does not buy independence —
pick engines that differ in what they assume.

## Citation

Wang W, Hu Z, Weiler P, Mayes S, Lange M, Wang J, Xue Z, Sauka-Spengler T, Theis FJ.
RegVelo: gene-regulatory-informed dynamics of single cells. *Cell* 2026.
doi:10.1016/j.cell.2026.04.022 · PMID 42119563
