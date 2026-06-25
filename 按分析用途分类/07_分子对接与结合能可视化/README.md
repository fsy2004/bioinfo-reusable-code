# 07 · Molecular Docking and Binding Energy Visualization

Modules for visualizing docking binding energies and running docking plus molecular dynamics pipelines.

| Module | Purpose | Language | Output figures | Status |
|------|------|------|--------|:---:|
| [022 Docking binding energy visualization](022_分子对接结合能可视化/) | Binding energy matrix to heatmap/bubble | R | Binding energy heatmap, bubble | Available |
| 086 Vina+GROMACS+MMPBSA | Automated docking and molecular dynamics | Python | RMSD/RMSF/Rg/MM-PBSA | Heavy environment |

Notes:
- **022**: produces a heatmap from a binding energy matrix.
- **086**: depends on external molecular dynamics toolchains including AutoDock Vina, GROMACS, gmx_MMPBSA, and MDAnalysis. It cannot be rendered locally; the original scripts are kept for reference.
- Follows the [shared framework conventions](../_framework/CONVENTIONS.md).
