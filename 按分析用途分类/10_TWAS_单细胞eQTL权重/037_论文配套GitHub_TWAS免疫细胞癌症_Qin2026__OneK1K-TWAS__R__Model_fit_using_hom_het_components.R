# =============================================================================
# 编号       : R037
# 脚本名     : Model_fit_using_hom_het_components.R
# 分类       : 10_TWAS_单细胞eQTL权重
# 项目来源   : 论文配套GitHub_TWAS免疫细胞癌症_Qin2026
# 用途       : 基于第一步得到的同/异质成分进一步拟合预测模型（弹性网），输出预测残差与系数，用于后续 FUSION TWAS。
# 结果图     : 未检测到明确作图输出
# 主要 R 包  : data.table; doParallel; foreach; glmnet
# 整理时间   : 2026-05-10
# =============================================================================
#'
#' @title Fitting models using shared and specific components generated in the first step
#'
#' @description Fitting models using shared and specific components generated in the first step.
#'
#' @param Gname Add this prefix to all the saved results files.
#' @return This function returns the estimated coefficients in the second step
#
#' @export
Model_fit <- function(Gname){
  out_dir <- "results/"
  X_file <- paste0("Exp/", Gname, "/snpexp")
  Y_file_dir <- paste0(out_dir, Gname, "_Y_residual")
  Y_hat_dir <- paste0(out_dir, Gname, "_predictors")

  Tissues <-  c("B_IN","B_Mem","CD4_ET","CD4_NC","CD4_SOX4","CD8_ET","CD8_NC","CD8_S100B","DC","Mono_C","Mono_NC","NK","NK_R","Plasma")

  suppressWarnings(expr = {X<-data.table::fread(file = X_file, sep='\t', data.table=F)})
  X<-as.matrix(data.frame(X, row.names=1, check.names = F))

  hom_expr_mat <- as.data.frame(data.table::fread(Y_file_dir))
  row.names(hom_expr_mat) <- hom_expr_mat$sampleid
  hom_expr_mat <- hom_expr_mat[,-1]
  
  q<-ncol(hom_expr_mat)
  m<-ncol(X)
  N<-nrow(X)
  if (q==2) {print(paste0("Only two tissues here, skip this gene!!!!!")); next}

  load(Y_hat_dir)
  Yhats_het <- Yhats_het
  Yhats_hom <- Yhats_hom
  Yhats_tiss <- Yhats_tiss
  
  content_alpha=0.5
  
  Yh_tiss <- vector("list", q)
  Yh_full<-vector("list",q)
  Yh_full_tiss1<-vector("list",q)

  for(i in 1:q){
      Yh_tiss[[i]]<-matrix(NA, ncol=2, nrow=nrow(Yhats_tiss[[i]]), 
                         dimnames = list(rownames(Yhats_tiss[[i]]), c("pred", "Y")))
    Yh_full[[i]]<-matrix(NA, ncol=2, nrow=nrow(Yhats_hom[[i]]), 
                         dimnames = list(rownames(Yhats_hom[[i]]), c("pred", "Y")))   
    Yh_full_tiss1[[i]]<-matrix(NA, ncol=2, nrow=nrow(Yhats_hom[[i]]), 
                               dimnames = list(rownames(Yhats_hom[[i]]), c("pred", "Y")))
  }

  tiss_betas <- vector("list", q)
  tiss_int <- vector("list", q)

  full_betas<-vector("list", q)
  full_int<-vector("list", q)

  full_tiss_betas<-vector("list", q)
  full_tiss_int <-vector("list", q) 
   
  R_abs <- matrix(NA, q, 3)
  colnames(R_abs) <- c("Tiss", "Full", "Full_tiss")
  rownames(R_abs) <- substr(names(Yhats_hom),8,nchar(names(Yhats_hom)))
  for (i in 1:q){
    
    fold_hom_expr_mat <- hom_expr_mat        
    test_hom_expr_mat <- hom_expr_mat

    fold_Yhats_hom <- Yhats_hom[[i]]
    fold_Yhats_het <- Yhats_het[[i]]
    fold_Yhats_tiss <- Yhats_tiss[[i]]
            
    fold_test_Yhats_hom<-Yhats_hom[[i]]
    fold_test_Yhats_het<-Yhats_het[[i]]    
    fold_test_Yhats_tiss<-Yhats_tiss[[i]]

    ## Tissue specific
    m2 <- lm(fold_hom_expr_mat[rownames(fold_Yhats_tiss),i] ~ as.numeric(fold_Yhats_tiss))
    Yhat1_tiss <- predict(m2, data.frame(fold_Yhats_tiss = as.numeric(fold_test_Yhats_tiss)))
    R_tiss <- abs(cor(Yhat1_tiss, test_hom_expr_mat[rownames(fold_Yhats_tiss),i]))

    tiss_betas[[i]]<- coef(m2)[-1]
    names(tiss_betas[[i]]) <- "tiss"
    tiss_int[[i]]<-coef(m2)[1]


    ## Content method (shared component + tissue specific)
    X4 <- data.frame(hom=as.numeric(fold_Yhats_hom), het=as.numeric(fold_Yhats_het))
    m4 <- lm(fold_hom_expr_mat[rownames(fold_Yhats_hom),i] ~ ., X4)
    Yhat1_full <- predict(m4, data.frame(hom = as.numeric(fold_test_Yhats_hom), het= as.numeric(fold_test_Yhats_het)))
    R_full <- abs(cor(Yhat1_full, test_hom_expr_mat[rownames(fold_test_Yhats_het),i]))
    
    full_betas[[i]]<-coef(m4)[-1]
    full_int[[i]]<-coef(m4)[1]

    Yh_tiss[[i]][rownames(fold_test_Yhats_hom),"pred"] <- Yhat1_tiss
    Yh_tiss[[i]][rownames(fold_test_Yhats_hom),"Y"] <- test_hom_expr_mat[rownames(fold_test_Yhats_hom),i]
 
    Yh_full[[i]][rownames(fold_test_Yhats_hom),"pred"] <- Yhat1_full
    Yh_full[[i]][rownames(fold_test_Yhats_hom),"Y"] <- test_hom_expr_mat[rownames(fold_test_Yhats_hom),i]
    

    # full model with across tissue information
    ##############################################
    Yhat_all <- data.frame(sampleid=rownames(Yhats_hom[[i]]), pred=Yhats_hom[[i]])
    colnames(Yhat_all)[2] <- "Hom"
    for (s in 1:length(Yhats_het)){
      Yhat1 <- as.data.frame(Yhats_het[[s]])
      Yhat1$sampleid <- rownames(Yhat1)
      colnames(Yhat1)[1] <- substr(names(Yhats_het[s]), 8, nchar(names(Yhats_het[s])))
      Yhat_all <- merge(Yhat_all, Yhat1, by="sampleid", all=T)
    }
    
    Y_value <- as.data.frame(hom_expr_mat[,i])
    Y_value$sampleid <- rownames(hom_expr_mat)
    colnames(Y_value)[1] <- "Y"
    Yhat_all_whole <- merge(Y_value, Yhat_all, by="sampleid", all.x=T)
    rownames(Yhat_all_whole) <- Yhat_all_whole$sampleid
    Yhat_all_whole <- Yhat_all_whole[is.na(apply(Yhat_all_whole[,-1], 1, mean))==F,-1]
    
    explanatory1=bigstatsr::as_FBM(Yhat_all_whole[,-1], backingfile=paste0(out_dir, Gname, "_", i,"_content1_tmp"))
    #print(paste0("Number of common noNA individuals: ", sum(is.na(apply(Yhat_all_whole, 1, mean))==F)))
    
    set.seed(i)
    full_tiss_fit1<-bigstatsr::big_spLinReg(X = explanatory1, 
                                ind.train = match(rownames(Yhat_all_whole), rownames(Yhat_all_whole)),
                                y.train =Yhat_all_whole[,"Y"], K=10, alphas = c(content_alpha),warn=F)
    Yhat1_full_tiss1 <- predict(full_tiss_fit1, explanatory1, 
                               ind.row = match(rownames(Yhat_all_whole), rownames(Yhat_all_whole)))
    Y1 <- test_hom_expr_mat[rownames(Yhat_all_whole),]
    
    R_full_tiss1 <- abs(cor(Yhat1_full_tiss1, Y1[,i]))
    Yh_full_tiss1[[i]][rownames(Y1), "pred"] <- Yhat1_full_tiss1
    Yh_full_tiss1[[i]][rownames(Y1), "Y"] <- Y1[,i]
    
    full_tiss_beta_vals1 <- unlist(summary(full_tiss_fit1)$beta[
      which.min(summary(full_tiss_fit1)$validation_loss)])[1:ncol(explanatory1)]
    full_tiss_int1 <- unlist(summary(full_tiss_fit1)$intercept[
      which.min(summary(full_tiss_fit1)$validation_loss)])
    if(sum(is.na(full_tiss_beta_vals1))>0){
      full_tiss_beta_vals1[which(is.na(full_tiss_beta_vals1))]<-0
    }
    if(length(attr(full_tiss_fit1, "ind.col"))<dim(explanatory1)[2]){
      new_betas<-rep(0, dim(explanatory1)[2])
      new_betas[attr(full_tiss_fit1, "ind.col")] <- full_tiss_beta_vals1[1:length(attr(full_tiss_fit1, "ind.col"))]
      full_tiss_beta_vals1<-new_betas
    }

    full_tiss_betas[[i]] <- full_tiss_beta_vals1
    full_tiss_int[[i]] <- full_tiss_int1


    names(full_tiss_betas[[i]]) <- colnames(Yhat_all_whole)[-1]

  }
  
  system(paste0("rm ", out_dir, Gname, "*.bk"))
  
  for (i in 1:q){
    R_abs[i,1] <- abs(cor(Yh_tiss[[i]], use="complete.obs")[1,2])
    R_abs[i,2] <- abs(cor(Yh_full[[i]], use="complete.obs")[1,2])
    R_abs[i,3] <- abs(cor(Yh_full_tiss1[[i]], use="complete.obs")[1,2])
  }
  
  save(tiss_betas, tiss_int, 
       full_betas, full_int, 
       full_tiss_betas, full_tiss_int, 
       file = paste0(out_dir,Gname,"_models_Betas"))

  save(Yh_tiss, Yh_full, Yh_full_tiss1, 
       file = paste0(out_dir,Gname,"_models_prediction"))

  #write.csv(R_abs, file=paste0(out_dir, Gname, "_abs_R_Genes.csv"))

}



