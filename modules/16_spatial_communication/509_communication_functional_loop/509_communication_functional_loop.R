# =============================================================================
# 509 · 细胞通讯功能复现闭环 (ligand activity -> UCell -> enrichment -> Venn)
# -----------------------------------------------------------------------------
# 把"配体活性预测 → 受体细胞功能打分 → 富集 → 共识靶点"串成一个可验证的闭环,
# 解决通讯分析常被诟病的"只给 LR 对、不证下游功能"问题:
#   ① NicheNet 式【配体活性】= 配体先验靶权重 与 受体细胞 DEG logFC 的相关(排序);
#   ② UCell 式【单细胞签名打分】= 对 top 配体的预测靶基因做秩和打分,按细胞群比较;
#   ③ 富集 = 预测靶 在受体 DEG 中的超几何检验(GSEA 的离散版);
#   ④ Venn 收口 = 预测靶 ∩ 受体 DEG ∩ 通路基因集 → 高可信"功能靶"。
#
# ★为 turnkey 免大模型下载,先验矩阵为合成;接真数据时把先验换成 NicheNet
#   ligand_target_matrix、DEG 换成受体细胞 FindMarkers 结果即可(见 README)。
# Turnkey: Rscript 509_communication_functional_loop.R   (合成→results/+assets/)
# 复用 _framework/theme_pub.R(含 venn_pub);无条形图(lollipop+violin+Venn)。
# =============================================================================
suppressWarnings(suppressMessages({ library(ggplot2) }))

## ---- 定位脚本 + 载框架 -----------------------------------------------------
.this <- tryCatch({ a <- commandArgs(FALSE); m <- grep("^--file=", a)
  if(length(m)) dirname(normalizePath(sub("^--file=","",a[m[1]]))) else getwd() },
  error=function(e) getwd())
.fw <- NULL; .p <- .this
for(i in 1:6){ cand <- file.path(.p,"_framework","theme_pub.R"); if(file.exists(cand)){ .fw <- cand; break }; .p <- dirname(.p) }
if(!is.null(.fw)) source(.fw) else stop("需要 _framework/theme_pub.R")

set.seed(42)
DIR <- .this
DDAT <- file.path(DIR,"example_data"); DRES <- file.path(DIR,"results"); DAST <- file.path(DIR,"assets")
for(d in c(DDAT,DRES,DAST)) dir.create(d, showWarnings=FALSE, recursive=TRUE)

## ---- UCell 式单细胞签名打分(秩和 / Mann-Whitney U 归一)-------------------
ucell_score <- function(expr, signature, maxRank=1500){
  # expr: 基因 x 细胞;对每个细胞按表达排名,签名基因的 U 统计归一为 0-1 分
  sig <- intersect(signature, rownames(expr)); ns <- length(sig)
  apply(expr, 2, function(x){
    r <- rank(-x, ties.method="average")           # 高表达=小秩
    r[r>maxRank] <- maxRank+1
    u <- sum(r[match(sig, rownames(expr))])
    1 - (u - ns*(ns+1)/2) / (ns*maxRank)            # UCell 归一(越大=签名越活跃)
  })
}

## ---- 1. 合成:受体 scRNA(近/远 sender 两群)+ 配体-靶先验 + sender 配体 ----
ng<-150; nc<-300; genes<-sprintf("G%03d",1:ng); cells<-sprintf("C%03d",1:nc)
grp <- rep(c("near_sender","distal"), each=nc/2)             # 近 sender 群被通讯激活
target_prog <- genes[1:25]                                   # 真·下游靶程序(25 基因)
E <- matrix(rpois(ng*nc, 5), ng, nc, dimnames=list(genes,cells))
E[target_prog, grp=="near_sender"] <- E[target_prog, grp=="near_sender"] +
  rpois(length(target_prog)*sum(grp=="near_sender"), 8)      # 近 sender 群靶程序上调
# 配体-靶先验矩阵(8 配体 x 150 基因);L1/L2 强权重落在 target_prog(真活性配体)
ligs <- sprintf("LIG%d",1:8); P <- matrix(rexp(8*ng,30), 8, ng, dimnames=list(ligs,genes))
P["LIG1",target_prog] <- P["LIG1",target_prog] + runif(25,0.5,0.9)
P["LIG2",target_prog] <- P["LIG2",target_prog] + runif(25,0.3,0.6)
write.csv(round(E,2), file.path(DDAT,"receiver_expr.csv"))
write.csv(round(P,4), file.path(DDAT,"ligand_target_prior.csv"))
write.csv(data.frame(cell=cells, group=grp), file.path(DDAT,"receiver_groups.csv"), row.names=FALSE)
cat(sprintf("[gen] synthetic: receiver %dx%d, %d ligands; true active=LIG1/LIG2 (demo only)\n", ng,nc,length(ligs)))

## ---- 2. 受体 DEG(near vs distal 的 logFC)---------------------------------
lfc <- log2((rowMeans(E[,grp=="near_sender"])+1)/(rowMeans(E[,grp=="distal"])+1))
deg <- names(sort(lfc, decreasing=TRUE))[lfc[order(lfc,decreasing=TRUE)] > 0.5]

## ---- 3. NicheNet 式配体活性 = 先验靶权重 vs DEG logFC 的相关,排序 ---------
activity <- apply(P, 1, function(w) cor(w, lfc, method="pearson"))
act_df <- data.frame(ligand=names(activity), activity=as.numeric(activity))
act_df <- act_df[order(-act_df$activity),]
top_lig <- act_df$ligand[1]
write.csv(act_df, file.path(DRES,"ligand_activity.csv"), row.names=FALSE)

## ---- 4. top 配体的预测靶 → UCell 打分(按细胞群)---------------------------
pred_targets <- names(sort(P[top_lig,], decreasing=TRUE))[1:30]
uc <- ucell_score(E, pred_targets)
uc_df <- data.frame(cell=cells, group=grp, UCell=uc)
write.csv(uc_df, file.path(DRES,"ucell_scores.csv"), row.names=FALSE)
wt <- wilcox.test(UCell~group, data=uc_df)

## ---- 5. 富集(预测靶 在 DEG 中的超几何检验)+ Venn 收口 -------------------
pathway <- genes[1:40]                                        # 假定的相关通路基因集(含真靶程序)
ov <- length(intersect(pred_targets, deg))
enr_p <- phyper(ov-1, length(deg), ng-length(deg), length(pred_targets), lower.tail=FALSE)
sets <- list(`Predicted targets`=pred_targets, `Receiver DEGs`=deg, `Pathway set`=pathway)
consensus <- Reduce(intersect, sets)
write.csv(data.frame(consensus_functional_target=consensus), file.path(DRES,"consensus_targets.csv"), row.names=FALSE)
cat(sprintf("[loop] top ligand=%s (activity r=%.2f); UCell near>distal p=%.2g; enrich p=%.2g; consensus=%d genes\n",
            top_lig, act_df$activity[1], wt$p.value, enr_p, length(consensus)))

## ---- 6. 出图(lollipop + violin + Venn;无条形图)-------------------------
col <- pal_pub(3,"npg")
# Fig1: 配体活性 lollipop(top 高亮)
act_df$ligand <- factor(act_df$ligand, levels=rev(act_df$ligand))
act_df$hl <- act_df$ligand==top_lig
p1 <- ggplot(act_df, aes(activity, ligand)) +
  geom_segment(aes(x=0, xend=activity, yend=ligand, color=hl), linewidth=1.2, alpha=0.7) +
  geom_point(aes(color=hl), size=4) +
  scale_color_manual(values=c(`FALSE`="grey65",`TRUE`=col[1]), guide="none") +
  labs(x="Ligand activity (prior-target vs DEG correlation)", y=NULL,
       title="NicheNet-style ligand activity") + theme_pub(base_size=11)
save_fig(p1, file.path(DAST,"ligand_activity"), width=5.6, height=3.4)

# Fig2: UCell 预测靶签名 violin+box+jitter(近 vs 远)
p2 <- ggplot(uc_df, aes(group, UCell, fill=group)) +
  geom_violin(alpha=0.55, color=NA, trim=FALSE) +
  geom_boxplot(width=0.18, outlier.shape=NA, alpha=0.9) +
  geom_jitter(width=0.08, size=0.5, alpha=0.35) +
  scale_fill_manual(values=col[2:3], guide="none") +
  labs(x=NULL, y="UCell score (predicted-target signature)",
       title=sprintf("Receiver activity by group (p=%.2g)", wt$p.value)) +
  theme_pub(base_size=11)
save_fig(p2, file.path(DAST,"ucell_violin"), width=4.4, height=3.6)

# Fig3: 闭环 Venn(预测靶 ∩ DEG ∩ 通路)
p3 <- venn_pub(sets, title=sprintf("Closed-loop functional targets (n=%d)", length(consensus)))
save_fig(p3, file.path(DAST,"closedloop_venn"), width=5.2, height=4.4)

cat("[fig] assets/: ligand_activity, ucell_violin, closedloop_venn (.pdf+.png)\n")
sink(file.path(DRES,"sessionInfo.txt")); print(sessionInfo()); sink()
