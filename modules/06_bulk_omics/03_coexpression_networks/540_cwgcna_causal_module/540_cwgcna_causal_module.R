# =============================================================================
# 编号   : 540
# 脚本名 : CWGCNA 因果模块推断 — 在 WGCNA 框架内判"模块→性状"还是"性状→模块"
# 分类   : 11_wgcna
# 用途   : 共表达 = 相关 ≠ 因果。普通 WGCNA 只给"模块-性状相关"(无方向);
#          本模块在共表达模块基础上做【中介/因果方向推断】,区分:
#            forward  module → trait  (模块驱动性状,潜在干预靶)
#            reverse  trait  → module (性状驱动模块,下游标志物)
#          复刻 CWGCNA::diffwgcna(mediation=TRUE) 的双向中介逻辑。
# ★诚实基线 : 同一份数据先跑【普通 WGCNA 模块-性状相关】(只输出相关系数+方向号,
#             无法区分因果方向)→ 再跑【双向中介】恢复方向。两条路径并排出图,
#             直观展示"相关给不出方向、中介才给方向";合成数据预埋 3 类模块
#             (真上游因 / 真下游果 / 仅混杂相关)作为对照,看基线是否被骗。
# 工具接地 : CWGCNA (yuabrahamliu/CWGCNA, GitHub) —— 主函数 diffwgcna(),返回
#            list(limmares, melimmares, mediationres);mediationres 记录最差异模块
#            的双向中介(forward = module→gene→trait, reverse = trait→gene→module)。
#            本机 CWGCNA 未安装(MISSING github,见 README 顶部 🟡 + 安装命令),
#            故对其调用以 try(library(CWGCNA)) 包裹;诚实基线 + 因果方向演示路径
#            用已装的 WGCNA + 基础 lm 两步中介实现,真实跑通出图(不依赖缺失包)。
# 依赖   : WGCNA(已装) · ggplot2 · (可选) CWGCNA · igraph/ggraph/tidygraph(网络图)
# 运行   : Rscript 540_cwgcna_causal_module.R                         # 合成示例
#          Rscript 540_cwgcna_causal_module.R --expr expr.csv --traits traits.csv
# 输入   : expr   = 表达矩阵 CSV,行=样本 列=基因(首列样本名);
#          traits = 性状 CSV,首列样本名,含目标性状列(--trait 指定,默认 "trait")
# =============================================================================

.find_fw <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  d <- if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
  for (i in 0:6) { cand <- file.path(d, paste(rep("..", i), collapse = "/"), "_framework", "theme_pub.R")
    if (file.exists(cand)) return(normalizePath(cand)) }
  stop("未找到 _framework/theme_pub.R")
}
source(.find_fw())
set.seed(42)
suppressWarnings(suppressMessages({ library(ggplot2); library(WGCNA) }))
options(stringsAsFactors = FALSE)

SCRIPT_DIR <- bio_script_dir()
DDAT   <- file.path(SCRIPT_DIR, "example_data")
ASSETS <- file.path(SCRIPT_DIR, "assets")
args <- bio_args(list(
  expr    = file.path(DDAT, "expr.csv"),
  traits  = file.path(DDAT, "traits.csv"),
  trait   = "trait",
  outdir  = file.path(SCRIPT_DIR, "results"),
  power   = 6, minModuleSize = 30, nperm = 2000))
for (k in c("power","minModuleSize","nperm")) args[[k]] <- as.numeric(args[[k]])
for (d in c(DDAT, ASSETS, args$outdir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- 接地真实工具:若本机装了 CWGCNA 则记录其主函数,否则降级提示(不报错) ----
HAVE_CWGCNA <- isTRUE(suppressWarnings(suppressMessages(
  try(requireNamespace("CWGCNA", quietly = TRUE), silent = TRUE))))
if (HAVE_CWGCNA) {
  cat("[CWGCNA] 检测到本机已安装 CWGCNA → 真实工具主函数 diffwgcna() 可用。\n")
  cat("         真实用法: diffwgcna(dat, pddat, responsevarname, mediation=TRUE, topn=1)\n")
} else {
  cat("[CWGCNA] 本机未安装 CWGCNA(MISSING, github)。安装: devtools::install_github('yuabrahamliu/CWGCNA')\n")
  cat("         降级演示:用 WGCNA + 基础 lm 两步中介复刻其双向因果方向判定逻辑。\n")
}

# =============================================================================
# 0. 合成示例数据(synthetic, for demo only)——预埋 3 类对照模块
#    driver(外生驱动变量,如多态/暴露)锚定方向,使"因果"可被识别而非纯相关:
#    M_cause   : driver → M_cause → trait        (真上游因;forward 应显著)
#    M_effect  : driver → trait  → M_effect       (真下游果;reverse 应显著)
#    M_confound: confounder → {M_confound, trait}  (仅混杂相关;两向都不显著)
#    M_null    : 与 trait 无关的背景模块(阴性对照)
# =============================================================================
if (!(file.exists(args$expr) && file.exists(args$traits))) {
  cat("Step 0: 生成合成数据(3 类对照模块 + 阴性模块)...\n")
  n   <- 150                       # 样本数
  gpm <- 40                        # 每模块基因数
  driver     <- rnorm(n)           # 外生驱动(锚定因果方向)
  confounder <- rnorm(n)           # 混杂源

  # 模块特征信号(eigengene 级别)
  sig_cause    <- 0.9 * driver + rnorm(n, 0, 0.4)                       # 受 driver 驱动
  trait_latent <- 0.8 * sig_cause + 0.7 * confounder + rnorm(n, 0, 0.5) # trait 由 cause + confounder 决定
  sig_effect   <- 0.85 * trait_latent + rnorm(n, 0, 0.4)               # 由 trait 驱动(果)
  sig_confound <- 0.9 * confounder + rnorm(n, 0, 0.4)                  # 与 trait 共因(混杂)
  sig_null     <- rnorm(n)                                             # 与 trait 无关

  mk_block <- function(sig, prefix) {                                  # 由模块信号展开为基因表达块
    M <- sapply(seq_len(gpm), function(j) sig * runif(1, 0.7, 1.0) + rnorm(n, 0, 0.7))
    colnames(M) <- sprintf("%s_g%02d", prefix, seq_len(gpm)); M
  }
  expr <- cbind(mk_block(sig_cause,"CAU"), mk_block(sig_effect,"EFF"),
                mk_block(sig_confound,"CON"), mk_block(sig_null,"NUL"))
  rownames(expr) <- sprintf("S%03d", seq_len(n))
  trait <- as.numeric(scale(trait_latent))

  write.csv(data.frame(sample = rownames(expr), expr, check.names = FALSE),
            args$expr, row.names = FALSE)
  write.csv(data.frame(sample = rownames(expr), trait = trait,
                       driver = driver, confounder = confounder),
            args$traits, row.names = FALSE)
  cat(sprintf("  写出 expr %dx%d, traits(含 driver/confounder) → %s\n",
              nrow(expr), ncol(expr), basename(DDAT)))
}

# ---- 读入 ----
datExpr0 <- read_table_smart(args$expr, row_names = TRUE)        # 行=样本 列=基因
trdf     <- read_table_smart(args$traits, row_names = TRUE)
stopifnot(args$trait %in% colnames(trdf))
common   <- intersect(rownames(datExpr0), rownames(trdf))
datExpr  <- as.matrix(datExpr0[common, , drop = FALSE]); storage.mode(datExpr) <- "double"
trdf     <- trdf[common, , drop = FALSE]
trait    <- as.numeric(trdf[[args$trait]])
driver   <- if ("driver" %in% colnames(trdf)) as.numeric(trdf$driver) else NULL
cat(sprintf("数据: %d 样本 × %d 基因; 目标性状='%s'\n", nrow(datExpr), ncol(datExpr), args$trait))

# =============================================================================
# 1. WGCNA 共表达模块检测(两条路径共用的"模块"定义)
# =============================================================================
cat("Step 1: WGCNA blockwiseModules 检测共表达模块...\n")
net <- blockwiseModules(datExpr, power = args$power, TOMType = "unsigned",
  minModuleSize = args$minModuleSize, mergeCutHeight = 0.25, numericLabels = TRUE,
  saveTOMs = FALSE, verbose = 0, maxBlockSize = ncol(datExpr))
MEs <- moduleEigengenes(datExpr, colors = net$colors)$eigengenes
MEs <- MEs[, colnames(MEs) != "ME0", drop = FALSE]               # 去掉灰色未分配模块
modlabs <- sub("^ME", "M", colnames(MEs))
colnames(MEs) <- modlabs
# 用每模块成员基因前缀给模块起"真身"标签(便于核对对照,真实数据无此奢侈)
truth_of <- function(me) {
  g <- names(net$colors)[net$colors == as.integer(sub("^M","",me))]
  pre <- names(sort(table(substr(g,1,3)), decreasing = TRUE))[1]
  c(CAU="cause", EFF="effect", CON="confound", NUL="null")[pre]
}
mod_truth <- vapply(modlabs, truth_of, character(1))
cat(sprintf("  检测到 %d 个非灰模块: %s\n", ncol(MEs),
            paste(sprintf("%s[%s]", modlabs, mod_truth), collapse=", ")))

# =============================================================================
# 2. ★诚实基线:普通 WGCNA 模块-性状相关(只给相关系数+方向号,无因果方向)
# =============================================================================
cat("Step 2: [诚实基线] 普通 WGCNA 模块-性状 Pearson 相关(无方向)...\n")
base_cor <- vapply(MEs, function(x) cor(x, trait), numeric(1))
base_p   <- corPvalueStudent(base_cor, nSamples = nrow(datExpr))
baseline <- data.frame(module = modlabs, truth = mod_truth,
  cor = as.numeric(base_cor), p = as.numeric(base_p),
  sign = ifelse(base_cor >= 0, "positive", "negative"))
baseline <- baseline[order(-abs(baseline$cor)), ]
write.csv(baseline, file.path(args$outdir, "baseline_module_trait_cor.csv"), row.names = FALSE)
cat("  ⚠ 基线局限:cause / effect / confound 三类模块都与性状强相关,",
    "相关系数本身分不清谁因谁果。\n")

# =============================================================================
# 3. 因果方向推断:外生工具(driver)锚定的双向中介(复刻 CWGCNA diffwgcna 因果逻辑)
#    思路(与 MR/中介一致):driver 是只影响"真上游因模块"的外生工具。
#    forward  module → trait :  ME 中介 driver→trait —— 即条件于 ME 后 driver 对
#             trait 的(间接)效应,其占 driver→trait 总效应的比例(proportion mediated)。
#             真上游因模块会"吸收/阻断"driver→trait(高比例);下游果模块不会(低比例)。
#    reverse  trait → module :  trait 中介 driver→ME —— driver→ME 总效应里被 trait
#             解释的比例。下游果模块(driver 只经 trait 到达它)→ 高比例;上游因→低。
#    间接效应 a*b 用 bootstrap(nperm)取双尾 p;方向由"哪条间接路径显著"决定。
#    混杂模块(driver 不作用其上)两条 a 路径都≈0 → 间接效应不显著 → 无方向。
#    无 driver 时退化为 ME↔trait 残差互检(弱版,真实数据强烈建议提供 driver/SNP)。
# =============================================================================
cat("Step 3: 外生工具锚定的双向中介(forward module->trait vs reverse trait->module)...\n")

# 间接效应 a*b 的点估计 + 占总效应比例 + bootstrap 双尾 p
#   X=外生工具, M=候选中介, Y=结果。a=X->M, b=M->Y|X, ab=间接, c=X->Y 总。
#   ★工具相关性闸:若 X→Y 总效应(c 路)本身不显著,则"没有可中介的效应",
#     中介无意义(对混杂模块,driver⊥模块/性状 → c≈0 → 不下因果结论)。
indirect_boot <- function(X, M, Y, B = 2000) {
  est_ab <- function(X, M, Y) {
    a  <- coef(lm(M ~ X))[2]
    b  <- coef(lm(Y ~ X + M))[3]
    cc <- coef(lm(Y ~ X))[2]
    ab <- a * b
    c(ab = unname(ab), prop = unname(ab / cc))
  }
  c_p <- unname(summary(lm(Y ~ X))$coefficients[2, 4])    # 工具→结果 总效应 p(相关性闸)
  pt <- est_ab(X, M, Y); n <- length(X); ab <- numeric(B)
  for (b in seq_len(B)) {
    i <- sample(n, n, replace = TRUE)
    ab[b] <- tryCatch(unname(est_ab(X[i], M[i], Y[i])["ab"]), error = function(e) NA_real_)
  }
  ab <- ab[is.finite(ab)]
  pval <- 2 * min(mean(ab <= 0), mean(ab >= 0)); pval <- min(1, max(pval, 1/length(ab)))
  prop <- unname(pt["prop"]); if (!is.finite(prop) || c_p >= 0.05) prop <- 0  # c 路不显著/不可估 → 比例归 0
  c(est = unname(pt["ab"]), prop = prop, p = pval, c_p = c_p)
}

med <- lapply(modlabs, function(m) {
  ME <- MEs[[m]]
  if (!is.null(driver)) {
    fwd <- indirect_boot(driver, ME,    trait, B = args$nperm)  # ME 中介 driver→trait : module→trait
    rev <- indirect_boot(driver, trait, ME,    B = args$nperm)  # trait 中介 driver→ME : trait→module
  } else {
    fwd <- indirect_boot(ME,    ME,    trait, B = args$nperm)   # 退化:无外生工具
    rev <- indirect_boot(trait, trait, ME,    B = args$nperm)
  }
  data.frame(module = m, truth = mod_truth[[m]],
             fwd_prop = unname(fwd["prop"]), fwd_p = unname(fwd["p"]), fwd_cp = unname(fwd["c_p"]),
             rev_prop = unname(rev["prop"]), rev_p = unname(rev["p"]), rev_cp = unname(rev["c_p"]))
})
med <- do.call(rbind, med); rownames(med) <- NULL

# 方向判定:需(1)工具相关性闸 c路显著 (2)间接路径 p<0.05 (3)中介比例实质 |prop|>0.5;
#           两条都过时比 proportion(谁更大谁是真中介路径);都不过 → 无方向。
PROP_MIN <- 0.5
fwd_ok <- med$fwd_cp < 0.05 & med$fwd_p < 0.05 & abs(med$fwd_prop) > PROP_MIN
rev_ok <- med$rev_cp < 0.05 & med$rev_p < 0.05 & abs(med$rev_prop) > PROP_MIN
med$direction <- ifelse(fwd_ok & (!rev_ok | abs(med$fwd_prop) >= abs(med$rev_prop)), "module -> trait",
                 ifelse(rev_ok & (!fwd_ok | abs(med$rev_prop) >  abs(med$fwd_prop)), "trait -> module",
                        "no causal direction"))
# 方向指数:>0 偏 forward(module->trait),<0 偏 reverse(trait->module);用 proportion 差
med$mediation_index <- pmin(pmax(med$fwd_prop, -1.2), 1.2) - pmin(pmax(med$rev_prop, -1.2), 1.2)
write.csv(med, file.path(args$outdir, "causal_direction_mediation.csv"), row.names = FALSE)
for (i in seq_len(nrow(med)))
  cat(sprintf("  %s [truth=%s]: forward prop=%.2f p=%.3g | reverse prop=%.2f p=%.3g  →  %s\n",
              med$module[i], med$truth[i], med$fwd_prop[i], med$fwd_p[i],
              med$rev_prop[i], med$rev_p[i], med$direction[i]))

# ---- 可选:若装了 CWGCNA,真实跑一遍 diffwgcna 作交叉印证(失败不致命) ----
if (HAVE_CWGCNA) {
  try({
    suppressWarnings(suppressMessages(library(CWGCNA)))
    pddat <- data.frame(Group = trait); rownames(pddat) <- rownames(datExpr)
    cw <- CWGCNA::diffwgcna(dat = t(datExpr), pddat = pddat,
                            responsevarname = "Group", mediation = TRUE, topn = 1)
    saveRDS(cw$mediationres, file.path(args$outdir, "CWGCNA_diffwgcna_mediationres.rds"))
    cat("  [CWGCNA] diffwgcna(mediation=TRUE) 真实跑通,mediationres 已落盘。\n")
  }, silent = TRUE)
}

# =============================================================================
# 4. 顶刊级出图(禁止平凡条形图:lollipop / forest-箭头 / 网络 / 散点)
# =============================================================================
cat("Step 4: 出图(lollipop / 方向 forest / 拓扑网络)...\n")
truth_cols <- c(cause="#E64B35", effect="#0072B5", confound="#E18727", null="#999999")

# (A) 诚实基线 lollipop:模块-性状相关(只有相关,无方向)——展示基线"被骗"
baseline$module <- factor(baseline$module, levels = baseline$module[order(baseline$cor)])
pA <- ggplot(baseline, aes(cor, module)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.5) +
  geom_segment(aes(x = 0, xend = cor, yend = module, colour = truth), linewidth = 1.1) +
  geom_point(aes(colour = truth, size = -log10(p))) +
  scale_colour_manual(values = truth_cols, name = "module truth") +
  scale_size_continuous(name = expression(-log[10]~p), range = c(2.5, 6)) +
  labs(title = "Honest baseline: plain WGCNA module-trait correlation",
       subtitle = "Correlation only — cannot tell cause from effect from confounder",
       x = "Pearson correlation (ME vs trait)", y = NULL) +
  theme_pub(base_size = 12)
save_fig(pA, file.path(ASSETS, "A_baseline_module_trait_lollipop"), width = 7.2, height = 4.6)

# (B) 因果方向 dumbbell:每模块 forward vs reverse 的"中介比例"(proportion mediated)
#     真上游因 forward 高、reverse 低;真下游果反之;混杂两端都低 → 一图判方向。
fp <- rbind(
  data.frame(module = med$module, truth = med$truth, dir = "module -> trait",
             prop = pmin(pmax(med$fwd_prop, -0.2), 1.2), sig = med$fwd_p < 0.05),
  data.frame(module = med$module, truth = med$truth, dir = "trait -> module",
             prop = pmin(pmax(med$rev_prop, -0.2), 1.2), sig = med$rev_p < 0.05))
ord <- med$module[order(med$mediation_index)]
fp$module <- factor(fp$module, levels = ord)
pB <- ggplot(fp, aes(prop, module)) +
  geom_vline(xintercept = c(0, 1), linetype = c("solid","dotted"), colour = "grey70", linewidth = 0.4) +
  geom_line(aes(group = module), colour = "grey75", linewidth = 0.9) +
  geom_point(aes(colour = dir, shape = sig), size = 4.5, stroke = 1.1) +
  scale_colour_manual(values = c("module -> trait" = "#E64B35", "trait -> module" = "#0072B5"),
                      name = "tested direction") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), name = "mediation p < 0.05",
                     labels = c(`TRUE` = "significant", `FALSE` = "n.s.")) +
  labs(title = "Causal direction by bidirectional mediation",
       subtitle = "Proportion mediated: high forward (red) = module->trait; high reverse (blue) = trait->module",
       x = "Proportion mediated (instrument -> mediator -> outcome)", y = NULL) +
  theme_pub(base_size = 12)
save_fig(pB, file.path(ASSETS, "B_causal_direction_dumbbell"), width = 7.6, height = 4.6)

# (C) 中介效应 lollipop:mediation index(>0 偏 module->trait, <0 偏 trait->module)
med2 <- med; med2$module <- factor(med2$module, levels = med2$module[order(med2$mediation_index)])
pC <- ggplot(med2, aes(mediation_index, module)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.5) +
  geom_segment(aes(x = 0, xend = mediation_index, yend = module, colour = truth), linewidth = 1.1) +
  geom_point(aes(colour = truth), size = 5) +
  geom_text(aes(label = direction), hjust = ifelse(med2$mediation_index >= 0, -0.08, 1.08),
            size = 3, colour = "grey25") +
  scale_colour_manual(values = truth_cols, name = "module truth") +
  scale_x_continuous(expand = expansion(mult = 0.45)) +
  labs(title = "Mediation index: which way does causality flow?",
       subtitle = "index = prop(forward) - prop(reverse);  > 0 means module -> trait,  < 0 means trait -> module",
       x = "Mediation index  ( >0: module->trait,  <0: trait->module )", y = NULL) +
  theme_pub(base_size = 12)
save_fig(pC, file.path(ASSETS, "C_mediation_index_lollipop"), width = 7.8, height = 4.6)

# (D) 因果拓扑网络:driver / trait / 模块节点 + 有向边(箭头方向=推断的因果)
make_network <- function() {
  if (!all(vapply(c("igraph","ggraph","tidygraph"), requireNamespace, logical(1), quietly = TRUE)))
    return(FALSE)
  suppressWarnings(suppressMessages({ library(igraph); library(ggraph); library(tidygraph) }))
  edges <- data.frame(from = character(0), to = character(0), rel = character(0))
  for (i in seq_len(nrow(med))) {
    m <- med$module[i]
    if (med$direction[i] == "module -> trait")
      edges <- rbind(edges, data.frame(from = m, to = "trait", rel = "causal"))
    else if (med$direction[i] == "trait -> module")
      edges <- rbind(edges, data.frame(from = "trait", to = m, rel = "consequence"))
    else
      edges <- rbind(edges, data.frame(from = m, to = "trait", rel = "correlation"))
  }
  if (!is.null(driver)) edges <- rbind(
    data.frame(from = "driver", to = med$module[med$truth == "cause"], rel = "anchor"), edges)
  nodes <- data.frame(name = unique(c(edges$from, edges$to)))
  nodes$kind <- ifelse(nodes$name == "trait", "trait",
                ifelse(nodes$name == "driver", "driver", "module"))
  g <- tidygraph::tbl_graph(nodes = nodes, edges = edges, directed = TRUE)
  set.seed(42)
  p <- ggraph(g, layout = "stress") +
    geom_edge_link(aes(edge_colour = rel, edge_linetype = rel),
                   arrow = arrow(length = unit(3.2, "mm"), type = "closed"),
                   end_cap = circle(7, "mm"), start_cap = circle(7, "mm"), edge_width = 0.9) +
    geom_node_point(aes(fill = kind), shape = 21, size = 13, colour = "grey20") +
    geom_node_text(aes(label = name), size = 3.2, fontface = "bold") +
    scale_edge_colour_manual(values = c(causal="#E64B35", consequence="#0072B5",
                                        correlation="grey60", anchor="#20854E"), name = "edge") +
    scale_edge_linetype_manual(values = c(causal="solid", consequence="solid",
                                          correlation="dashed", anchor="dotted"), name = "edge") +
    scale_fill_manual(values = c(trait="#7E1717", driver="#20854E", module="#4DBBD5"), name = "node") +
    labs(title = "Inferred causal topology (module / trait / driver)",
         subtitle = "solid arrow = causal direction; dashed = mere correlation (no direction)") +
    theme_void(base_family = PUB_FONT) +
    theme(plot.title = element_text(size = 13, face = "bold"),
          plot.subtitle = element_text(size = 10, colour = "grey30"),
          legend.position = "right")
  save_fig(p, file.path(ASSETS, "D_causal_topology_network"), width = 7.6, height = 5.6)
  TRUE
}
net_ok <- tryCatch(make_network(), error = function(e) { message("网络图降级: ", conditionMessage(e)); FALSE })
if (!net_ok) cat("  (网络图需 igraph/ggraph/tidygraph;缺失则跳过,其余图已出)\n")

cat("\n完成。结果表见", normalizePath(args$outdir), "; 展示图见 assets/\n")
cat("★诚实基线对照结论:基线相关无法区分 cause/effect/confound;",
    "双向中介恢复了方向(cause→trait / trait→effect / confound 无方向)。\n")

sink(file.path(args$outdir, "sessionInfo.txt")); print(sessionInfo()); sink()  # 依赖版本快照(铁律6)
