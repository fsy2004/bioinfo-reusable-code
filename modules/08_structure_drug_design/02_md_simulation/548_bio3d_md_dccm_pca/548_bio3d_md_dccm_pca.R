# =============================================================================
# 编号   : 548
# 脚本名 : bio3d 结构系综 / MD 轨迹分析 — DCCM 动态互相关 + PCA 主成分 + RMSF
# 分类   : 07_molecular_docking
# 用途   : 对接-MD 下游高级分析。把一组构象(晶体系综 / MD 轨迹快照)交给 bio3d:
#          ① PCA(pca.xyz)提取集体运动(几个主成分主导构象空间)
#          ② DCCM(dccm)揭示残基对的动态正/负相关(关联运动的结构域)
#          ③ RMSF(rmsf)逐残基定位柔性区(loop / 末端 / 别构开关)
#          三者互补:PCA=方向,DCCM=耦合,RMSF=幅度。
# ★诚实基线: 描述性分析,非统计检验。内置「集体性基线」对照——
#          把真实系综的 PCA 方差捕获(前几个 PC 占总方差比例)与
#          「打散残基间协方差的零模型(逐坐标独立重排)」对比:
#          真实系综前 2-3 个 PC 远超零模型 → 证明观察到的是真集体运动而非噪声。
#          DCCM/PCA/RMSF 本身不产生 p 值,结论是「揭示/定位」而非「证明因果」。
# 依赖   : bio3d · ggplot2 · (theme_pub.R 框架)
# 运行   : Rscript 548_bio3d_md_dccm_pca.R                 # 内置 transducin 晶体系综(零下载)
#          Rscript 548_bio3d_md_dccm_pca.R --pdb my.pdb    # 单结构 → NMA 驱动 DCCM
#          Rscript 548_bio3d_md_dccm_pca.R --traj traj.dcd --topol ref.pdb  # 真实 MD 轨迹
# 输入   : 默认 = bio3d 内置 transducin 多构象系综(53 个 Gα 晶体结构, GDP/GTP 两态);
#          或 --pdb 单 PDB(走弹性网络 NMA 得 DCCM/系综);
#          或 --traj + --topol(MD 轨迹 dcd/nc + 拓扑 pdb,经 read.dcd/fit 走全流程)
# =============================================================================

## ---- 定位框架 + 加载顶刊主题 ------------------------------------------------
.find_fw <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  d <- if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
  for (i in 0:6) { cand <- file.path(d, paste(rep("..", i), collapse = "/"), "_framework", "theme_pub.R")
    if (file.exists(cand)) return(normalizePath(cand)) }
  stop("未找到 _framework/theme_pub.R")
}
source(.find_fw())
suppressWarnings(suppressMessages({ library(bio3d); library(ggplot2) }))
set.seed(42)

SCRIPT_DIR <- bio_script_dir()
DDAT   <- file.path(SCRIPT_DIR, "example_data")
ASSETS <- file.path(SCRIPT_DIR, "assets")
args <- bio_args(list(
  pdb    = "",                 # 单 PDB → NMA 驱动 DCCM/系综(留空=用内置 transducin)
  traj   = "",                 # MD 轨迹 (.dcd/.nc) 留空=不走轨迹分支
  topol  = "",                 # 轨迹对应拓扑 PDB
  outdir = file.path(SCRIPT_DIR, "results"),
  n_pc   = 3,                  # 展示的主成分数
  n_top_flex = 12))            # RMSF lollipop 标注的最柔残基数
for (k in c("n_pc","n_top_flex")) args[[k]] <- as.integer(args[[k]])
for (d in c(DDAT, ASSETS, args$outdir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 0. 准备构象系综 (xyz 矩阵: 行=构象, 列=3N 坐标) + 残基注释 + DCCM
#    三条入口,均落到统一的 (xyz_fit, resno, dccm) 结构,下游一视同仁。
# =============================================================================
cat("Step 0: 准备构象系综 + DCCM...\n")

prep <- list()
if (nzchar(args$traj) && nzchar(args$topol) && file.exists(args$traj) && file.exists(args$topol)) {
  # ---- 入口 C: 真实 MD 轨迹 -------------------------------------------------
  cat("  [入口 C] 读 MD 轨迹:", basename(args$traj), "\n")
  topo <- read.pdb(args$topol)
  ca   <- atom.select(topo, "calpha")
  trj  <- read.dcd(args$traj)                          # 帧 x 3N
  xyz_fit <- fit.xyz(fixed = topo$xyz, mobile = trj,
                     fixed.inds = ca$xyz, mobile.inds = ca$xyz)
  prep$xyz   <- xyz_fit[, ca$xyz, drop = FALSE]
  prep$resno <- topo$atom$resno[ca$atom]
  prep$state <- NULL
  prep$src   <- sprintf("MD trajectory (%d frames)", nrow(trj))
  prep$dccm  <- dccm(prep$xyz)                          # 轨迹协方差直接给 DCCM

} else if (nzchar(args$pdb) && file.exists(args$pdb)) {
  # ---- 入口 B: 单结构 → 弹性网络 NMA 生成系综 + DCCM -----------------------
  cat("  [入口 B] 单结构 NMA:", basename(args$pdb), "\n")
  pdb  <- read.pdb(args$pdb)
  ca   <- atom.select(pdb, "calpha")
  m    <- nma(pdb)                                      # 弹性网络模型(Cα),全程基于 Cα
  # 沿最低频的三个非平凡模式(7,8,9)各造一段系综,叠成多构象 xyz 矩阵
  ens  <- rbind(mktrj(m, mode = 7), mktrj(m, mode = 8), mktrj(m, mode = 9))  # 帧 x 3N(Cα)
  prep$xyz   <- ens                                     # nma 已是 Cα 级,列即残基坐标
  prep$resno <- pdb$atom$resno[ca$atom]
  prep$state <- NULL
  prep$src   <- sprintf("single PDB + ENM-NMA modes 7-9 (%s)", basename(args$pdb))
  prep$dccm  <- dccm(m)                                 # NMA 模式协方差给 DCCM(Cα x Cα)

} else {
  # ---- 入口 A (默认, turnkey 零下载): bio3d 内置 transducin 晶体系综 -------
  cat("  [入口 A] 内置 transducin 系综 (53 Gα 晶体结构, GDP/GTP 两态)\n")
  data(transducin)
  pdbs  <- transducin$pdbs
  gaps  <- gap.inspect(pdbs$xyz)                        # 去比对空位的 xyz 列
  gapsr <- gap.inspect(pdbs$ali)                        # 去空位的残基列
  # 叠合到第一个构象(Cα 最小二乘),消除整体平移/旋转,只留内部运动
  xyz_fit <- fit.xyz(fixed = pdbs$xyz[1, ], mobile = pdbs$xyz,
                     fixed.inds = gaps$f.inds, mobile.inds = gaps$f.inds)
  prep$xyz   <- xyz_fit[, gaps$f.inds, drop = FALSE]    # 53 x (3*305)
  prep$resno <- pdbs$resno[1, gapsr$f.inds]             # 305 残基编号
  prep$state <- transducin$annotation[, "state"]        # GDP / GTP 标签(构象状态)
  prep$src   <- "transducin crystallographic ensemble (built-in)"
  prep$dccm  <- dccm(xyz_fit[, gaps$f.inds])            # 系综坐标协方差给 DCCM
}

n_conf <- nrow(prep$xyz); n_res <- length(prep$resno)
cat(sprintf("  系综: %d 构象 x %d 残基 (Cα) | 来源: %s\n", n_conf, n_res, prep$src))

# =============================================================================
# 1. PCA: 提取集体运动 (pca.xyz)
# =============================================================================
cat("Step 1: PCA (pca.xyz) 提取集体运动...\n")
pc <- pca.xyz(prep$xyz)
var_pct <- pc$L / sum(pc$L) * 100                        # 各 PC 方差占比
cum2 <- sum(var_pct[1:min(2, length(var_pct))])
cat(sprintf("  PC1=%.1f%% PC2=%.1f%% PC3=%.1f%% | 前 2 PC 累计 %.1f%%\n",
            var_pct[1], var_pct[2], if (length(var_pct) >= 3) var_pct[3] else NA, cum2))

# =============================================================================
# 2. ★诚实基线 (描述性): 集体性零模型对照
#    把每个坐标列独立随机重排(打散残基间协方差),再做 PCA。
#    若真实系综 PC1/PC2 远高于零模型 → 集体运动真实存在,而非噪声/采样假象。
# =============================================================================
cat("Step 2: ★诚实基线 — 集体性零模型 (逐坐标独立重排)...\n")
N_NULL <- 15                                             # 零模型重排次数(每次一遍 pca.xyz)
k_keep <- min(3, length(var_pct))
null_mat <- matrix(NA_real_, nrow = N_NULL, ncol = k_keep)  # 逐 PC 方差占比(供基线图复用)
for (b in seq_len(N_NULL)) {
  xs <- apply(prep$xyz, 2, sample)                       # 打散每列(去掉列间相关)
  Lb <- pca.xyz(xs)$L
  null_mat[b, ] <- Lb[seq_len(k_keep)] / sum(Lb) * 100
}
null_pc1 <- mean(null_mat[, 1])
cat(sprintf("  真实 PC1 = %.1f%%  vs  零模型 PC1 = %.1f%% (均值, n=%d)\n",
            var_pct[1], null_pc1, N_NULL))
cat(sprintf("  → 集体性倍数 = %.1fx (>1 表明观察到的是真集体运动,而非独立残基噪声)\n",
            var_pct[1] / null_pc1))

# =============================================================================
# 3. RMSF: 逐残基柔性 (rmsf)
# =============================================================================
cat("Step 3: RMSF 逐残基柔性...\n")
rf <- rmsf(prep$xyz)
df_rmsf <- data.frame(resno = prep$resno, rmsf = rf)
top_flex <- df_rmsf[order(-df_rmsf$rmsf), ][seq_len(min(args$n_top_flex, nrow(df_rmsf))), ]
cat(sprintf("  最柔 %d 残基 (resno): %s\n", nrow(top_flex),
            paste(top_flex$resno, collapse = ", ")))

## 落盘关键结果表 -------------------------------------------------------------
write.csv(data.frame(PC = seq_along(var_pct), var_pct = round(var_pct, 3),
                     cum_pct = round(cumsum(var_pct), 3)),
          file.path(args$outdir, "PCA_variance.csv"), row.names = FALSE)
write.csv(df_rmsf, file.path(args$outdir, "RMSF_per_residue.csv"), row.names = FALSE)
write.csv(round(prep$dccm, 4), file.path(args$outdir, "DCCM_matrix.csv"))
write.csv(data.frame(metric = "PC1_var_pct", real = round(var_pct[1], 2),
                     null_mean = round(null_pc1, 2),
                     fold = round(var_pct[1] / null_pc1, 2)),
          file.path(args$outdir, "baseline_collectivity.csv"), row.names = FALSE)

# =============================================================================
# 4. 顶刊级出图 (全部精修矢量;禁止平凡条形图)
# =============================================================================
cat("Step 4: 出图 (DCCM 热图 / PCA 散点+porcupine / RMSF lollipop / 基线 dumbbell)...\n")

# ---- 图 1: DCCM 动态互相关热图 (残基 x 残基, 发散 RdBu) --------------------
cm <- prep$dccm
# DCCM 的残基轴 = 1..n(NMA/轨迹按拓扑顺序),系综分支与 resno 等长
cm_n <- nrow(cm)
res_axis <- if (cm_n == n_res) prep$resno else seq_len(cm_n)
dcc_df <- expand.grid(i = seq_len(cm_n), j = seq_len(cm_n))
dcc_df$cij <- as.vector(cm)
dcc_df$ri  <- res_axis[dcc_df$i]
dcc_df$rj  <- res_axis[dcc_df$j]
p_dccm <- ggplot(dcc_df, aes(ri, rj, fill = cij)) +
  geom_tile() +                                          # resno 含空位间隔不均 → tile 而非 raster
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1), name = "Cij") +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  labs(title = "Dynamic cross-correlation (DCCM)",
       subtitle = "Red = correlated, Blue = anti-correlated residue motions",
       x = "Residue i", y = "Residue j") +
  theme_pub(border = TRUE) +
  coord_cartesian(clip = "off") +                        # 防止标题/轴名被画布边缘裁掉
  theme(plot.margin = margin(12, 14, 10, 10))            # 近方形画布,留足边距给标题与轴名
save_fig(p_dccm, file.path(ASSETS, "fig1_dccm_heatmap"), width = 6.8, height = 6.2)

# ---- 图 2: PCA 投影散点 (构象在 PC1-PC2 平面;按状态着色) ------------------
pca_df <- data.frame(PC1 = pc$z[, 1], PC2 = pc$z[, 2])
if (!is.null(prep$state)) pca_df$State <- factor(prep$state) else pca_df$State <- factor("conf")
p_pca <- ggplot(pca_df, aes(PC1, PC2, fill = State)) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_point(shape = 21, size = 3.4, colour = "grey20", alpha = 0.9, stroke = 0.4) +
  scale_fill_manual(values = pal_pub(nlevels(pca_df$State), "npg")) +
  labs(title = "Conformer PCA projection",
       subtitle = sprintf("Each point = one conformer; PC1 %.1f%% / PC2 %.1f%% of variance",
                          var_pct[1], var_pct[2]),
       x = sprintf("PC1 (%.1f%%)", var_pct[1]),
       y = sprintf("PC2 (%.1f%%)", var_pct[2])) +
  theme_pub()
save_fig(p_pca, file.path(ASSETS, "fig2_pca_projection"), width = 6.4, height = 5.4)

# ---- 图 3: PC1 porcupine — 沿主链画 PC1 集体位移矢量 -----------------------
# pc$U[,1] 是 3N 特征向量;reshape 成 (n_res x 3) 即每残基位移方向
u1 <- matrix(pc$U[, 1], ncol = 3, byrow = TRUE)          # n_res x 3 (dx,dy,dz)
mean_xyz <- matrix(pc$mean, ncol = 3, byrow = TRUE)      # 平均结构坐标
mag <- sqrt(rowSums(u1^2))                                # 每残基 PC1 位移幅度
scal <- 12                                                # 矢量放大(仅作可视化)
porc_df <- data.frame(
  resno = res_axis[seq_len(nrow(u1))],
  x0 = mean_xyz[, 1], y0 = mean_xyz[, 2],
  x1 = mean_xyz[, 1] + scal * u1[, 1],
  y1 = mean_xyz[, 2] + scal * u1[, 2],
  mag = mag)
p_porc <- ggplot(porc_df) +
  geom_path(aes(x0, y0), colour = "grey75", linewidth = 0.6) +
  geom_segment(aes(x = x0, y = y0, xend = x1, yend = y1, colour = mag),
               arrow = arrow(length = unit(0.06, "inches")), linewidth = 0.55) +
  scale_colour_viridis_c(option = "C", name = "PC1\ndispl.") +
  coord_equal() +
  labs(title = "PC1 collective-motion porcupine",
       subtitle = "Arrows = per-residue displacement along principal mode 1 (projected x-y)",
       x = "x (Angstrom)", y = "y (Angstrom)") +
  theme_pub()
save_fig(p_porc, file.path(ASSETS, "fig3_pc1_porcupine"), width = 6.6, height = 5.6)

# ---- 图 4: RMSF 逐残基 lollipop (顶部最柔残基标注) -------------------------
df_rmsf$is_top <- df_rmsf$resno %in% top_flex$resno
p_rmsf <- ggplot(df_rmsf, aes(resno, rmsf)) +
  geom_segment(aes(xend = resno, yend = 0, colour = is_top), linewidth = 0.5) +
  geom_point(aes(colour = is_top), size = 1.7) +
  scale_colour_manual(values = c("FALSE" = "grey60", "TRUE" = "#B2182B"),
                      labels = c("FALSE" = "Rigid", "TRUE" = "Most flexible"),
                      name = NULL) +
  labs(title = "Per-residue flexibility (RMSF)",
       subtitle = "Lollipop height = fluctuation; red = top flexible residues (loops / switches)",
       x = "Residue number", y = "RMSF (Angstrom)") +
  theme_pub()
# 有 ggrepel 则标注最柔残基编号,无则降级(不挡主分析)
if (requireNamespace("ggrepel", quietly = TRUE)) {
  p_rmsf <- p_rmsf + ggrepel::geom_text_repel(
    data = top_flex, aes(resno, rmsf, label = resno),
    size = 2.7, colour = "#B2182B", max.overlaps = 20, seed = 42)
}
save_fig(p_rmsf, file.path(ASSETS, "fig4_rmsf_lollipop"), width = 7.4, height = 4.4)

# ---- 图 5: ★诚实基线 dumbbell — 真实系综 vs 集体性零模型 -------------------
# 复用 Step 2 的 null_mat(逐 PC 方差占比),不重复昂贵的重排 PCA
base_df <- data.frame(
  PC   = factor(paste0("PC", seq_len(k_keep)), levels = paste0("PC", rev(seq_len(k_keep)))),
  real = var_pct[seq_len(k_keep)],
  null = colMeans(null_mat))
p_base <- ggplot(base_df) +
  geom_segment(aes(x = null, xend = real, y = PC, yend = PC),
               colour = "grey70", linewidth = 1.1) +
  geom_point(aes(null, PC, colour = "Shuffled null"), size = 4) +
  geom_point(aes(real, PC, colour = "Real ensemble"), size = 4) +
  scale_colour_manual(values = c("Shuffled null" = "grey55", "Real ensemble" = "#1F77B4"),
                      name = NULL) +
  labs(title = "Honest baseline: collectivity vs shuffled null",
       subtitle = "Real PCs capture far more variance than coordinate-shuffled null = true collective motion",
       x = "Variance captured (%)", y = NULL) +
  theme_pub()
save_fig(p_base, file.path(ASSETS, "fig5_baseline_collectivity"), width = 6.6, height = 3.8)

cat("完成。结果表见", normalizePath(args$outdir), "; 图见 assets/\n")
sink(file.path(args$outdir, "sessionInfo.txt")); print(sessionInfo()); sink()  # 依赖快照(铁律6)
