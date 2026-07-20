# =============================================================================
# 518 · beyondcell 单细胞药物响应 (drug response heterogeneity from scRNA)
# -----------------------------------------------------------------------------
# 复现 beyondcell 核心:用药物扰动签名(LINCS/PSc 的 up/down 基因集)对每个细胞算
# beyondcell score(BCS = UCell(up) − UCell(down)),据 BCS 谱把细胞分成"治疗簇"
# (therapeutic clusters),并按组间 BCS 差异给药物排序 → 找出能区分敏感/耐药亚群的药。
#
# ★实现方法核心(turnkey 免装 beyondcell 重包);接真数据时签名用 beyondcell 的
#   PSc/SSc 集合或自建 LINCS 签名(模块070 chemCPA / 015 可衔接)。
# Turnkey: Rscript 518_beyondcell_drug_response.R   (合成 scRNA+药物签名→results/+assets/)
# 复用 _framework/theme_pub.R;无条形图(heatmap + lollipop + UMAP)。
# =============================================================================
suppressWarnings(suppressMessages({ library(ggplot2) }))

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

ucell <- function(expr, sig, maxRank=800){               # UCell 式秩和打分
  sig<-intersect(sig,rownames(expr)); ns<-length(sig); if(ns<2) return(rep(0,ncol(expr)))
  apply(expr,2,function(x){ r<-rank(-x,ties.method="average"); r[r>maxRank]<-maxRank+1
    1-(sum(r[match(sig,rownames(expr))])-ns*(ns+1)/2)/(ns*maxRank) }) }

## ---- 1. 合成:scRNA(敏感/耐药两态)+ 6 药 up/down 签名 -------------------
ng<-220; nc<-300; genes<-sprintf("G%03d",1:ng)
state<-rep(c("sensitive","resistant"),each=nc/2)
E<-matrix(rpois(ng*nc,4),ng,nc,dimnames=list(genes,sprintf("C%03d",1:nc)))
sens_prog<-genes[1:30]; res_prog<-genes[31:60]
E[sens_prog,state=="sensitive"]<-E[sens_prog,state=="sensitive"]+rpois(30*sum(state=="sensitive"),6)
E[res_prog, state=="resistant"]<-E[res_prog, state=="resistant"]+rpois(30*sum(state=="resistant"),6)
drugs<-sprintf("Drug%d",1:6)
sig_up<-list(Drug1=sens_prog[1:20], Drug2=sens_prog[5:24],          # 1-2 命中敏感态
             Drug3=res_prog[1:20],  Drug4=res_prog[5:24],          # 3-4 命中耐药态
             Drug5=sample(genes,20), Drug6=sample(genes,20))        # 5-6 无特异
sig_dn<-lapply(sig_up,function(u) sample(setdiff(genes,u),20))
saveRDS(list(E=E,state=state,up=sig_up,dn=sig_dn), file.path(DDAT,"beyondcell_demo.rds"))
cat(sprintf("[gen] synthetic scRNA %dx%d (sensitive/resistant) + %d drug signatures (demo only)\n",ng,nc,length(drugs)))

## ---- 2. BCS = UCell(up) − UCell(down),每细胞每药 ------------------------
BCS<-sapply(drugs,function(d) ucell(E,sig_up[[d]])-ucell(E,sig_dn[[d]]))  # 细胞 x 药
rownames(BCS)<-colnames(E)

## ---- 3. 治疗簇(对 BCS 谱 kmeans)----------------------------------------
set.seed(42); tc<-kmeans(scale(BCS),centers=2,nstart=20)$cluster
tc<-paste0("TC",tc)
# 让 TC 命名稳定:含更多 sensitive 的簇记为 TC-sens
tab<-table(tc,state); sens_tc<-rownames(tab)[which.max(tab[,"sensitive"])]
tc<-ifelse(tc==sens_tc,"TC-sensitive","TC-resistant")

## ---- 4. 药物排序(组间 BCS 差异:敏感态 − 耐药态)------------------------
drug_diff<-data.frame(drug=drugs,
  delta=sapply(drugs,function(d) mean(BCS[state=="sensitive",d])-mean(BCS[state=="resistant",d])),
  p=sapply(drugs,function(d) wilcox.test(BCS[,d]~state)$p.value))
drug_diff<-drug_diff[order(-abs(drug_diff$delta)),]
write.csv(drug_diff, file.path(DRES,"drug_ranking.csv"), row.names=FALSE)
tc_mean<-aggregate(BCS, list(TC=tc), mean)
write.csv(tc_mean, file.path(DRES,"BCS_by_therapeutic_cluster.csv"), row.names=FALSE)
cat(sprintf("[beyondcell] top differential drug: %s (Δ=%.2f); TCs recovered: %s\n",
            drug_diff$drug[1], drug_diff$delta[1], paste(unique(tc),collapse="/")))

## ---- 5. 出图(heatmap + lollipop + UMAP;无条形图)-----------------------
col<-pal_pub(3,"npg")
# Fig1: BCS 热图(治疗簇均值,药 x TC,RdBu)
hd<-reshape(tc_mean,direction="long",varying=list(2:(length(drugs)+1)),v.names="BCS",
            times=drugs,timevar="drug",idvar="TC")
hd$drug<-factor(hd$drug,levels=drug_diff$drug)
p1<-ggplot(hd,aes(TC,drug,fill=BCS))+geom_tile(color="white",linewidth=0.6)+
  scale_fill_gradient2(low="#2166AC",mid="#F7F7F7",high="#B2182B",midpoint=0,name="Mean BCS")+
  labs(x=NULL,y=NULL,title="beyondcell score by therapeutic cluster")+theme_pub(base_size=11)
save_fig(p1,file.path(DAST,"bcs_heatmap"),width=4.2,height=3.8)

# Fig2: 药物差异 lollipop(敏感−耐药)
drug_diff$drug<-factor(drug_diff$drug,levels=rev(drug_diff$drug))
drug_diff$sig<-drug_diff$p<0.05
p2<-ggplot(drug_diff,aes(delta,drug))+
  geom_vline(xintercept=0,linetype=2,color="grey60")+
  geom_segment(aes(x=0,xend=delta,yend=drug,color=sig),linewidth=1.2,alpha=0.7)+
  geom_point(aes(color=sig),size=4)+
  scale_color_manual(values=c(`FALSE`="grey65",`TRUE`=col[1]),labels=c("ns","p<0.05"),name=NULL)+
  labs(x="ΔBCS (sensitive − resistant)",y=NULL,title="Differential drug response")+theme_pub(base_size=11)
save_fig(p2,file.path(DAST,"drug_lollipop"),width=5.2,height=3.4)

# Fig3: UMAP(PCA 2D 近似)着色 top 药 BCS
pca<-prcomp(scale(t(E)))$x[,1:2]
ud<-data.frame(D1=pca[,1],D2=pca[,2],BCS=BCS[,as.character(drug_diff$drug[1])],TC=tc)
p3<-ggplot(ud,aes(D1,D2,color=BCS))+geom_point(size=1.3,alpha=0.85)+
  scale_color_gradient2(low="#2166AC",mid="grey90",high="#B2182B",midpoint=0,name="BCS")+
  labs(x="PC1",y="PC2",title=sprintf("BCS of %s (top drug)",drug_diff$drug[1]))+
  theme_pub(base_size=11)+theme(axis.text=element_blank(),axis.ticks=element_blank())
save_fig(p3,file.path(DAST,"bcs_umap"),width=4.8,height=3.8)

cat("[fig] assets/: bcs_heatmap, drug_lollipop, bcs_umap (.pdf+.png)\n")
sink(file.path(DRES,"sessionInfo.txt")); print(sessionInfo()); sink()
