# 497 · scSurvival — 单细胞队列生存(虚拟生存模块)

- **来源**：cliffren/scSurvival（v1.3.0，2026 活跃；GitHub 接入包装），下载 2026-06-06。repo: https://github.com/cliffren/scSurvival
- **语言/环境**：Python ≥3.11，PyTorch/GPU（torch+cu，scanpy/lifelines/sklearn）。安装 `pip install -e .`（GPU 版 torch 需手动装）。

## 用途（它干嘛）
从**单细胞队列**做**生存分析**、细胞级风险画像：VAE 学 batch-invariant 单细胞表征 → 聚合到病人级 → **multi-head attention 多示例 Cox 回归**；不仅出病人级生存风险预测,还**识别与生存风险最相关的细胞亚群**及其风险倾向。

## 硬性输入（用它的前提）
1. **单细胞队列**：多病人、每人贡献细胞（`adata.obs['sample']`，例 30 病人）。
2. **病人级生存数据**：`surv_time` + `surv_status`（**时间-事件 + 删失**）。
3. 可选：病人协变量、batch。

## ★ 适用题型
带 **TCGA/队列生存随访(OS/PFS/DFS)** + **单细胞队列** 的**癌症预后题**。= "单细胞 → 预后"的虚拟生存模块；和 Scissor（058，单细胞↔二元/连续表型）互补——**Scissor 配二元表型、scSurvival 配生存随访**。
> 例:已弃的 disulfidptosis×ccRCC（TCGA-KIRC 生存 + ccRCC scRNA 队列)本可用它;留给下一个同型癌症预后题。

## ✗ 不适配（别硬上=跑不起来的装饰）
- **二元风险事件 / 无随访**的题（如 **UC→VTE**：结局是 VTE 二元风险,无 time-to-event,且无单细胞队列）→ 缺生存标签 + 缺队列,**两样都没有,加进来跑不起来**。
- 无单细胞队列(只有 1-2 个 scRNA 数据集)的题。
- **原则**:生存工具(scSurvival / Mime-101 496)需 time-to-event;配错题=负分。先看题有没有"单细胞队列 + 生存随访",有才上。

## 入口
`README.md`（框架图 pics/）、`examples/`（quick start）、`scSurvival/`（package）。
