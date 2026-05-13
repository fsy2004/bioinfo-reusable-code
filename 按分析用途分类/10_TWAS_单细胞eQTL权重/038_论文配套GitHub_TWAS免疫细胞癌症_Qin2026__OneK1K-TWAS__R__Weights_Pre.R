# =============================================================================
# 编号       : R038
# 脚本名     : Weights_Pre.R
# 分类       : 10_TWAS_单细胞eQTL权重
# 项目来源   : 论文配套GitHub_TWAS免疫细胞癌症_Qin2026
# 用途       : 合并前两步的系数，形成 FUSION TWAS 所需 weights 文件准备步骤。
# 结果图     : 未检测到明确作图输出
# 主要 R 包  : data.table
# 整理时间   : 2026-05-10
# =============================================================================
#'
#' @title Combing coefficients from above two steps and prepare for weights file
#'
#' @description Combing coefficients from above two steps and prepare for weights file.
#'
#' @param Gname Add this prefix to all the saved results files.
#' @param in_dir Direction for files stored from about two steps.
#' @return This function returns the combined coefficients from above two steps. 
#
#' @export
Weights_pre <- function(Gname, in_dir){
  load(paste0(in_dir, Gname, "_beta"))

  het_betas_snp <- het_tiss_betas
  het_int_snp <- het_tiss_ints
  tiss_betas_snp <- tiss_betas
  tiss_int_snp <- tiss_ints
  hom_betas_snp <- hom_beta_vals
  hom_int_snp <- hom_int
  snp_pos <- snps
  X_snp <- X

  colnames(snp_pos) <- c("chr", "snpid", "info", "pos", "A1", "A2")
  rm(het_tiss_betas)
  rm(het_tiss_ints)
  rm(tiss_betas)
  rm(tiss_ints)
  rm(hom_beta_vals)
  rm(hom_int)

  load(paste0(in_dir, Gname, "_predictors"))
  Yhats_het_snp <- Yhats_het
  Yhats_hom_snp <- Yhats_hom
  Yhats_tiss_snp <- Yhats_tiss

  load(paste0(in_dir, Gname, "_models_Betas"))
  load(paste0(in_dir, Gname, "_models_prediction"))

  tiss_betas <- tiss_betas
  tiss_int <- tiss_int
  full_betas <- full_betas
  full_int <- full_int
  full_tiss_betas <- full_tiss_betas
  full_tiss_int <- full_tiss_int

  Yhats_tiss_model <- Yh_tiss
  Yhats_full_model <- Yh_full
  Yhats_full_tiss_model <- Yh_full_tiss1
  rm(Yh_tiss)
  rm(Yh_full)
  rm(Yh_full_tiss1)

  p=length(Yhats_het_snp)
  tiss_beta_list <-vector("list",p)
  tiss_int_list <-vector("list",p)
  full_beta_list <-vector("list",p)
  full_int_list <-vector("list",p)
  full_tiss_beta_list <-vector("list",p)
  full_tiss_int_list <-vector("list",p)
  names(tiss_beta_list) <- names(Yhats_het_snp)
  names(tiss_int_list) <- names(Yhats_het_snp)
  names(full_beta_list) <- names(Yhats_het_snp)
  names(full_int_list) <- names(Yhats_het_snp)
  names(full_tiss_beta_list) <- names(Yhats_het_snp)
  names(full_tiss_int_list) <- names(Yhats_het_snp)


  for (i in 1:p){
    print(paste0(i, "-th tiss"))

    Yhats_het_snp1 <- Yhats_het_snp[[i]]
    Yhats_hom_snp1 <- Yhats_hom_snp[[i]]  
    Yhats_tiss_snp1 <- Yhats_tiss_snp[[i]]  
 
    het_betas_snp_1 <- het_betas_snp[[i]] 
    het_int_snp_1 <- het_int_snp[[i]]
    tiss_betas_snp_1 <- tiss_betas_snp[[i]] 
    tiss_int_snp_1 <- tiss_int_snp[[i]]
    hom_betas_snp <- hom_betas_snp
    hom_int_snp <- hom_int_snp
  
    # Yhats_het_check <- (as.matrix(X_snp[rownames(Yhats_het_snp1),])) %*% het_betas_snp_1 + het_int_snp_1
    # Yhats_hom_check <- (as.matrix(X_snp[rownames(Yhats_hom_snp1),])) %*% hom_betas_snp   + hom_int_snp
    # Yhats_tiss_check <- (as.matrix(X_snp[rownames(Yhats_tiss_snp1),])) %*% tiss_betas_snp_1 + tiss_int_snp_1

    het_beta_combi_val <- data.frame(beta=as.numeric(het_betas_snp_1), 
                                               snpid=snp_pos$snpid,
                                               chr=snp_pos$chr, 
                                               pos=snp_pos$pos)
    tiss_beta_combi_val <- data.frame(beta=as.numeric(tiss_betas_snp_1), 
                                               snpid=snp_pos$snpid,
                                               chr=snp_pos$chr, 
                                               pos=snp_pos$pos)
    hom_beta_combi_val <- data.frame(beta=as.numeric(hom_betas_snp),
                                               snpid=snp_pos$snpid,
                                               chr=snp_pos$chr, 
                                               pos=snp_pos$pos)
    het_int_combi_val <- het_int_snp_1
    tiss_int_combi_val <- tiss_int_snp_1
    hom_int_combi_val <- hom_int_snp
  
    # Yhats_het_check <- (as.matrix(X_snp[rownames(Yhats_het_snp1),])) %*% as.matrix(het_beta_combi_val[,1]) + het_int_combi_val
    # Yhats_tiss_check <- (as.matrix(X_snp[rownames(Yhats_tiss_snp1),])) %*% as.matrix(tiss_beta_combi_val[,1]) + tiss_int_combi_val
    # Yhats_hom_check <- (as.matrix(X_snp[rownames(Yhats_hom_snp1),])) %*% as.matrix(hom_beta_combi_val[,1]) + hom_int_combi_val
  
    tiss_beta_list[[i]] <- tiss_beta_combi_val
    tiss_int_list[[i]] <- tiss_int_combi_val

    ## Full model
    full_betas_1 <- full_betas[[i]]
    full_int_1 <- full_int[[i]] 
    Yhats_full_model_1 <- as.matrix(Yhats_full_model[[i]][,"pred"])
  
    coef_full_model <- matrix(rep(full_betas_1, each=length(het_betas_snp_1)), nrow=length(het_betas_snp_1))
    coef_full_snp <- matrix(c(hom_betas_snp, het_betas_snp_1), nrow=length(het_betas_snp_1))
  
    full_beta_combi_val <- apply(coef_full_model * coef_full_snp, 1, sum)
    full_int_combi_val <- sum(full_betas_1 * c(hom_int_snp, het_int_snp_1)) + full_int_1
  
    # Yhats_full_model_check <- as.matrix(X_snp[rownames(Yhats_het_snp1),]) %*% full_beta_combi_val + full_int_combi_val
    #Yhats_full_model_check1 <- Yhats_hom_1 * full_betas_1[1] + Yhats_het_snp1 * full_betas_1[2] + full_int_1
  
    full_beta_combi_val <- data.frame(beta=full_beta_combi_val,
                                               snpid=snp_pos$snpid,
                                               chr=snp_pos$chr, 
                                               pos=snp_pos$pos)
    # Yhats_full_model_check <- as.matrix(X_snp[rownames(Yhats_full_model_1),]) %*% as.matrix(full_beta_combi_val[,1]) + full_int_combi_val
  
  
    full_beta_list[[i]] <- full_beta_combi_val
    full_int_list[[i]] <- full_int_combi_val
  
    ## Full + tissues model
    full_tiss_betas_1 <- full_tiss_betas[[i]]
    full_tiss_int_1 <- full_tiss_int[[i]]
    Yhats_full_tiss_model_1 <- as.matrix(Yhats_full_tiss_model[[i]][,"pred"])
  
    coef_full_tiss_model <- matrix(full_tiss_betas_1)
    coef_full_tiss_snp <- cbind(hom_betas_snp, do.call(cbind, het_betas_snp))
  
    ## dim(coef_full_tiss_snp)=c(6489, 15)    dim(coef_full_tiss_model)=c(15, 1)
    full_tiss_beta_combi_val <- as.numeric(coef_full_tiss_snp %*% coef_full_tiss_model)
  
    int_full_tiss_snp <- c(hom_int_snp, do.call(cbind, het_int_snp))
    full_tiss_int_combi_val <- sum(full_tiss_betas_1 * int_full_tiss_snp) + full_tiss_int_1
  
    # Yhats_full_tiss_model_check <- as.matrix(X_snp[rownames(Yhats_full_tiss_model_1),]) %*% full_tiss_beta_combi_val + full_tiss_int_combi_val
  
    full_tiss_beta_combi_val <- data.frame(beta=full_tiss_beta_combi_val, 
                                               snpid=snp_pos$snpid,
                                               chr=snp_pos$chr, 
                                               pos=snp_pos$pos)

    # Yhats_full_tiss_model_check <- as.matrix(X_snp[rownames(Yhats_full_tiss_model_1),]) %*% as.matrix(full_tiss_beta_combi_val[,1]) + full_tiss_int_combi_val
  
    full_tiss_beta_list[[i]] <- full_tiss_beta_combi_val
    full_tiss_int_list[[i]] <- full_tiss_int_combi_val
  
  
  }


  save(tiss_beta_list, tiss_int_list, 
     full_beta_list, full_int_list, 
     full_tiss_beta_list, full_tiss_int_list, 
     file = paste0(in_dir,Gname,"_models_Betas_snps"))

}
