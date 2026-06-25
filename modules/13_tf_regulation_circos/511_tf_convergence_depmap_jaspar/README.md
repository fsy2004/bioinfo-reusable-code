# 511 · TF convergence — regulon × JASPAR motif × DepMap essentiality

Triangulates three **orthogonal** lines of evidence into a single rankable score for
"is this transcription factor a genuine core regulator?", guarding against the
single-evidence false positives that a lone SCENIC (or lone motif, or lone essentiality)
hit tends to produce.

| | |
|---|---|
| Language / deps | R · `ggplot2` `ggrepel` (+ shared `theme_pub.R`) |
| Purpose | Robust core-TF shortlist from convergent multi-evidence |
| Input | `tf_evidence.csv` (regulon / motif / DepMap); synthetic by default |
| Output | `results/tf_convergence.csv`; preview in `assets/` |

## Method

For each TF, three evidences are rank-normalised to 0–1 and averaged into a convergence
score; TFs scoring high on **all three** are flagged convergent:

1. **Regulon activity** — SCENIC/pySCENIC regulon AUCell mean in the cells of interest.
2. **Motif support** — fraction of the regulon's target promoters carrying the TF's
   JASPAR binding motif.
3. **DepMap essentiality** — CRISPR gene-effect (more negative = more essential),
   sign-flipped so larger = more essential.

## Input

`tf_evidence.csv` columns: `TF`, `regulon_activity`, `motif_support`,
`depmap_gene_effect`. Demo data is synthetic (12 TFs; TF01–03 designed as true
convergent, TF04/05/06 as single-evidence traps), generated on first run.

**On real data**: regulon activity from pySCENIC (module 081), motif support from
JASPAR2024 matching against regulon targets, essentiality from the DepMap CRISPR
(Chronos) gene-effect table.

## Use

A defensible final filter after regulon inference: instead of trusting one SCENIC run,
require that a candidate master TF is also motif-supported and functionally essential —
the kind of multi-evidence convergence reviewers expect for a "key regulator" claim.

## Outputs

| File | Type | Description |
|------|------|------|
| `results/tf_convergence.csv` | table | per-TF evidences, rank scores, convergence, flag |
| `assets/tf_convergence_scatter.png` | scatter | motif × regulon, colour = essentiality, convergent ringed |
| `assets/evidence_heatmap.png` | heatmap | TF × three rank-normalised evidences |
| `assets/convergence_lollipop.png` | lollipop | convergence ranking (convergent TFs highlighted) |

![TF convergence](assets/tf_convergence_scatter.png)

## Run

```bash
Rscript 511_tf_convergence_depmap_jaspar.R
```

## Dependencies

```r
install.packages(c("ggplot2","ggrepel"))
```
