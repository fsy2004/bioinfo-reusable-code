# 基准与评述 —— 用新方法前先过这一关

这里收的不是工具，是**尺子和警告**：顶刊上专门评测/批评某类方法的论文。
库里加了一批 2025–2026 的新方法（基础模型、虚拟扰动、空间深度学习），
它们几乎都自称超过前人。这些论文的作用，是让我们在写进论文之前先问一句
「跟朴素基线比过了吗、用的指标站得住吗」。

对应本库两条既有纪律：
- `CONVENTIONS.md` §0 —— 复用真实工具，不从头写；
- 每个新方法模块都必须带**可跑的朴素基线**（见各模块 README 的「特点」段）。

全部条目的 PMID 均经 NCBI E-utilities 核实（2026-07-20）。

---

## 1. 单细胞基础模型：先证明它打得过 PCA

| 论文 | 出处 | PMID | 结论 |
|---|---|---|---|
| Zero-shot evaluation reveals limitations of single-cell foundation models | Genome Biol 2025 Apr 18 | 40251685 | 零样本设置下，单细胞基础模型未能稳定超过 PCA 等简单基线 |
| Biology-driven insights into the power of single-cell foundation models | Genome Biol 2025 Oct 3 | 41044630 | 从生物学角度重估基础模型能力，提出 scGraph-OntoRWR 等以本体为基准的评估指标 |

**怎么用**：`01_single_cell/10_foundation_models/` 与
`02_spatial_transcriptomics/07_foundation_models/` 下的模块（scPRINT、Nicheformer、
EpiAgent、CAPTAIN、Novae 等），任何"基础模型嵌入更好"的说法，都要在同一份数据上
跟 PCA / Harmony / scVI 比过再写。零样本尤其要小心。

## 2. 整合基准：silhouette 不是可靠指标

| 论文 | 出处 | PMID | 结论 |
|---|---|---|---|
| Shortcomings of silhouette in single-cell integration benchmarking | Nat Biotechnol 2026 Jun（在线 2025-07-30） | 40739072 | 指出 silhouette 类指标在单细胞整合基准中的缺陷 |

**怎么用**：评价整合效果（`01_single_cell/02_integration_batch/`）不要只报 silhouette
或 ASW。配合 `565_scmultibench_integration_benchmark` 的多指标评估，并明确写出用了哪些指标。

## 3. 扰动预测：泛化性才是关键

| 资源 | 出处 | 位置 |
|---|---|---|
| scPerturBench（27 种扰动响应预测方法的泛化性基准） | Nature Methods 2026 | `03_virtual_perturbation/05_benchmark/590_scperturbench_generalization` |
| scArchon | Genome Biology 2026 | `03_virtual_perturbation/05_benchmark/591_scarchon_perturbation_benchmark` |

**怎么用**：虚拟扰动（CellOracle、GEARS、RegVelo、scCausalVI 等）给出的排序，
在写进论文前应当有一个"换个数据/换个扰动还成不成立"的检验。已有的做法是跨先验、
跨队列复现；scPerturBench 提供了同口径的评估指标。

## 4. Geneformer：官方推荐流程

| 论文 | 出处 | PMID |
|---|---|---|
| Discovery of candidate therapeutic targets with Geneformer | Nat Protoc 2026 Apr 23 | 42026145 |

**怎么用**：`03_virtual_perturbation/01_insilico_knockout/507_geneformer_insilico`
的参数与流程，以这份 Nature Protocols 为准，不要照抄网上零散示例。

---

## 一句话原则

> 新方法进论文的门槛不是"发在顶刊上"，而是"在**我们自己的数据**上，
> 跟一个**朴素基线**比过，用**站得住的指标**"。
> 上面这几篇就是干这个用的。
