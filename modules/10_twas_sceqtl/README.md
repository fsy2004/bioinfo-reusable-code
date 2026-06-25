# 10 · TWAS / single-cell eQTL weights

Weight training and model fitting for transcriptome-wide association studies (TWAS).

| Module | Purpose | Language | Status |
|------|------|------|:---:|
| 036-039 OneK1K TWAS weights | sc-eQTL homogeneous/heterogeneous component fitting, weight preprocessing and generation | R | Heavy |
| 040-042 FUSION TWAS | FUSION targetC / S+targetC / S+allC models | R | Heavy |

The TWAS workflow depends on FUSION, plink, large LD reference panels, and sc-eQTL weight files. It is computationally intensive and relies on an external toolchain, so it is not rendered locally; the original scripts are kept for reference. Upstream input comes from the 09 category (GWAS processing); after weight generation, association testing uses the official FUSION workflow. For figure conventions, see [unified framework](../_framework/CONVENTIONS.md).
