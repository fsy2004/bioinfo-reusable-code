# 04_comorbidity_network

Disease-pair comorbidity network (co-occurrence association, igraph network,
Louvain communities, centrality).

| Module | Purpose | Language |
|--------|---------|----------|
| [530 comorbidity network](530_comorbidity_network/) | 2×2 phi/OR/Jaccard → igraph → Louvain modules + hubs | R |

Turnkey on a synthetic per-patient disease table; igraph construction grounded in
`../99_external_sources/comorbidity_networks/` and `CSB-IG_Comorbidity_Networks/`.
See the module README for the directed-network / Louvain / zero-cell caveats.
