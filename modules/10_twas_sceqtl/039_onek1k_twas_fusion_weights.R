# =============================================================================
# 编号       : R039
# 脚本名     : Gen_weights.R
# 分类       : 10_twas_sceqtl
# 项目来源   : 论文配套GitHub_TWAS免疫细胞癌症_Qin2026
# 用途       : 按细胞类型生成 FUSION TWAS 所需 weights .RDat 文件（targetC / S_targetC / S_allC 三种权重矩阵）。
# 结果图     : 未检测到明确作图输出
# 主要 R 包  : data.table
# 整理时间   : 2026-05-10
# =============================================================================
#'
#' @title Generating weights files required for FUSION TWAS
#'
#' @description Generating weights files required for FUSION TWAS.
#'
#' @param Gname Add this prefix to all the saved results files.
#' @param in_dir Direction for files stored from about two steps.
#' @param out_dir Direction for output weights files.
#' @return This function generates weights files required for FUSION TWAS. 
#
#' @export
Weights_gen <- function(Gname, in_dir, out_dir){
  
  load(paste0(in_dir, Gname,"_beta"))
  load(paste0(in_dir, Gname, "_models_Betas_snps"))

  cv.performance <- data.frame(hom=c(1.0, 1.0), tiss=c(1.0, 1.0), full=c(1.0, 1.0), full_tiss=c(1.0, 1.0))
  rownames(cv.performance) <- c("rsq", "pval")
  cv.performance <- as.matrix(cv.performance)

  for (i in 1:length(full_tiss_beta_list)){

  	#allcells1=c("CD4_ET","NK","CD4_NC","CD8_S100B","CD8_ET","B_IN","CD8_NC","B_Mem","NK_R","Mono_NC","Mono_C","DC","Plasma","CD4_SOX4")
  	TissID<- substr(names(full_tiss_beta_list)[[i]], 8, nchar(names(full_tiss_beta_list)[[i]]))
  
    tiss_path <- paste0(out_dir, "/", TissID)
  	if (!dir.exists(tiss_path)) {
  		# If the folder doesn't exist, create it
  		dir.create(tiss_path, recursive = TRUE)
	  } 

    #if (file.exists(paste0(out_dir, "/weights/", TissID, "/BC_", Gname, "_wgt.RDat"))) {print(paste0("Skip (file exits) !!!!!")); next}
    snps$V2 <- paste0(snps$V1, "_", snps$V4)
    wgt.matrix1 <- full_tiss_beta_list[[i]]
    # snp_pos1 <- snp_pos_sub[paste0(wgt.matrix1$chr, "_", wgt.matrix1$pos), ]

    wgt.matrix1$snpid <- paste0(wgt.matrix1$chr, "_", wgt.matrix1$pos)

    # hom_beta_list[[i]]$snpid <- paste0(hom_beta_list[[i]]$chr, "_", hom_beta_list[[i]]$pos)
  	full_beta_list[[i]]$snpid <- paste0(full_beta_list[[i]]$chr, "_", full_beta_list[[i]]$pos)
  	tiss_beta_list[[i]]$snpid <- paste0(tiss_beta_list[[i]]$chr, "_", tiss_beta_list[[i]]$pos)
  	# full_tiss_beta_list[[i]]$snpid <- paste0(full_tiss_beta_list[[i]]$chr, "_", full_tiss_beta_list[[i]]$pos)

  	# wgt.matrix2 <- merge(wgt.matrix1, hom_beta_list[[i]][, c("snpid", "beta")], by="snpid" ,all.x=T)
  	wgt.matrix2 <- merge(wgt.matrix1, full_beta_list[[i]][, c("snpid", "beta")], by="snpid" ,all.x=T)
  	wgt.matrix2 <- merge(wgt.matrix2, tiss_beta_list[[i]][, c("snpid", "beta")], by="snpid" ,all.x=T)
  	colnames(wgt.matrix2)[c(2, 5:6)] <- c("S_allC","S_targetC", "targetC")
    rownames(wgt.matrix2) <- wgt.matrix2$snpid
  	wgt.matrix <- as.matrix(wgt.matrix2[, c(6,5,2)])
  	wgt.matrix[is.na(wgt.matrix)] <- 0

  	hsq <- c(0.1, 0.1)
  	hsq.pv <- 0.01

    #snps1 <- strsplit(snp_pos1$SNP, "_")
    #snps2 <- t(data.frame(snps1))
     
    # snps <- data.frame(chr=snp_pos1$Chr, snpid=snp_pos1$SNP, status=0, 
    #                    pos=snp_pos1$Pos, P0=snp_pos1$A2, P1=snp_pos1$A1)
    colnames(snps) <- paste0("V", 1:6)
    rownames(snps) <- snps$V2

  	save(cv.performance, snps, wgt.matrix, hsq, hsq.pv, file=paste0(out_dir, "/", TissID, "/", Gname, "_wgt.RDat"))
  }
  unlink(paste0(in_dir, "*"))
}


