# =============================================================================
# 编号       : R036
# 脚本名     : hom_het_components_fit.R
# 分类       : 10_twas_sceqtl
# 项目来源   : 论文配套GitHub_TWAS免疫细胞癌症_Qin2026
# 用途       : OneK1K 单细胞 eQTL 共享/特异成分模型拟合：通过弹性网回归对每个细胞类型估计同质（hom）+ 异质（het）成分系数，作为 TWAS 第一步。
# 结果图     : 未检测到明确作图输出
# 主要 R 包  : data.table; doParallel; foreach; glmnet
# 整理时间   : 2026-05-10
# =============================================================================
#'
#' @title Fitting shared and specific components for each tissue or cell type
#'
#' @description Fitting shared and specific components for each tissue or cell type.
#'
#' @param X A genotype matrix, individuals are row names, no column names. 
#' @param Ys Expression data matrix, individuals are row names, no column names.
#' @param covs Covariates data.
#' @param X_file Direction containing genotype matrix, individuals are row names, no column names.
#' @param Y_file_dir Direction containing only the expression files, individuals are row names, no column names.
#' @param cov_file_dir Specify this directory if you want to residualize the expression. Should share the same prefix as the expression files
#' @param verbose Provide more information. 
#' @param out_dir Output CV predictors and other statistics.
#' @param snps A tab-delimited .txt file containing snp information. six columns.
#' @param signif Nominally significant p-value threshold (alpha). default=0.1
#' @param gene_name Add this prefix to all the saved results files. default=NULL
#' @param content_alpha Regularization constant, 1 means LASSO regression.
#' @param tissue_alpha Regularization constant for the tissue by tissue apprach.
#' @param num_folds Number of folds for cross-validation.
#' @param scale_exp Center and scale the gene expression
#' @param seed Seed for analysis
#' 
#' @return This function returns the estimated shared and specific components and other coefficients
#
#' @export
hom_het_model <- function(X=NULL, 
                    Ys=NULL, 
                    covs=NULL, 
                    X_file=NULL,  # Genotype matrix, individuals are row names, no column names.
                    Y_file_dir=NULL, # Direction containing only the expression files, individuals are row names, no column names.
                    cov_file_dir=NULL, # Cov file direction, should share the same prefix as the expression files
                    verbose=T,
                    out_dir=NULL,  # output CV predictors and other statistics
                    snps=NULL,  # A tab-delimited txt file containing snp information. six columns.
                    signif=0.1,
                    gene_name=NULL, # Add this prefix to all the saved results files
                    content_alpha=0.5,  # Regularization constant, 1 means LASSO regression.
                    tissue_alpha=0.5,  # Regularization constant for the tissue by tissue apprach.
                    num_folds=10,  # Number of folds for cross-validatation
                    scale_exp=F,  # Center and scale the gene expression
                    seed=9000){
  
  #print(paste0("content_alpha=", content_alpha))
  ## Open data and get directories set up
  # if(verbose){
  #   message("Saving cross-validated predictors and performance metrics in ", out_dir)
  # }
  # if(is.null(twas_dir)){
  #   twas_dir=paste0(out_dir,"TWAS_weights")
  #   if(!dir.exists(twas_dir)){
  #     system(paste0("mkdir ", twas_dir)) 
  #   }
  #   if(verbose){
  #     message("Saving TWAS weights in ", twas_dir) 
  #   }
  # }else{
  #   if(verbose){
  #     message("Saving TWAS weights in ", twas_dir) 
  #   }
  # }
  ## By default, try to run content, but if there is only one context, do not
  runcontent=T
  
  # Check if X is a matrix or a path
  ## If matrix, just make sure it's good to go
  if(class(X) %in% c("matrix", "data.frame")){
    if(verbose){
      message("Converting objects to matrix...")
    }
    X<-as.matrix(X)
    Ys<-lapply(Ys, as.matrix)
    if(length(Ys)<3){
      runcontent=F
    }
    covs<-lapply(covs, as.matrix)
    ## If X is a path, then everything is else treated as such
    }else if(class(X_file) == "character"){
      if(verbose){
      message("Reading in files...") 
    }
    suppressWarnings(expr = {X<-data.table::fread(file = X_file, sep='\t', data.table=F)})
    X<-as.matrix(data.frame(X, row.names=1, check.names = F))
    if(!is.null(cov_file_dir)){
      covs<-vector("list", length = length(list.files(cov_file_dir)))
      names(covs)=list.files(cov_file_dir)
      for(i in 1:length(covs)){
        suppressWarnings(expr = {covs[[i]]<-data.table::fread(file = paste0(cov_file_dir,list.files(cov_file_dir)[i]), 
                                                  sep='\t', data.table=F)})
        dts=sapply(covs[[i]], class)
        if (any(dts == "character") & min(which(dts == "character") != 1)) {
          stop("Covariate(s) ",
               paste0(which(dts == "character"), collapse=", "),
               " is(are) character(s). Please convert to numeric using one-hot
             encoding or as an integer.")
        }
        covs[[i]]<-as.matrix(data.frame(covs[[i]], row.names=1, check.names = F))
      }
    }
    Ys<-vector("list", length = length(list.files(Y_file_dir)))
    if(length(Ys)<3){
      runcontent=F
      return()
      message("Only 1 context found. CONTENT will not be run. Skip this gene")
    }
    names(Ys)<-list.files(Y_file_dir)
    for(i in 1:length(Ys)){
      suppressWarnings(expr={ Ys[[i]]<-data.table::fread(file = 
                                               paste0(Y_file_dir,list.files(Y_file_dir)[i]), sep='\t', data.table=F)})
      Ys[[i]]<-as.matrix(data.frame(Ys[[i]], row.names=1, check.names = F))
      if(any(is.na(Ys[[i]])) | any(is.nan(Ys[[i]]))){
        remove=unique(c( which(is.na(Ys[[i]])), which(is.nan(Ys[[i]])) ))
        Ys[[i]]=Ys[[i]][-remove,,drop=F]
      }
      if(scale_exp){
        Ys[[i]]=scale(Ys[[i]])
      }
    }
  }else{
    message("X is neither a matrix nor a path. Please supply X as a matrix, dataframe, 
            or a path.")
  }
  
  if(is.null(snps)){
    message("snps file is null, TWAS weights will need to be updated after CONTENT is run...")
  }else if(class(snps) == "character"){
    snps=read.table(snps, sep='\t', check.names = F, stringsAsFactors = F)
    storesnps=snps
  }
  if(is.null(gene_name)){
    message("gene_name not supplied, using 'gene' as gene name.
          results will not be easily distinguishable.")
  }
  
  # Weight hom Y based on cells from all CTs in each individual
  # Yw <- fread(file=Y_W_file_dir, sep="\t", data.table=F)
  # Yw <- as.matrix(data.frame(Yw, row.names=1, check.names=F))
  # if(any(is.na(Yw)) | any(is.nan(Yw))){
  #   remove=unique(c( which(is.na(Yw)), which(is.nan(Yw)) ))
  #   Yw=Yw[-remove,,drop=F]
  # }
   
  # snps1 <- snps  
  # snps <- snps1[paste0(snps1$V1,"_",snps1$V4, "_", snps1$V5, "_", snps1$V6) %in% GWAS_snps$V1,]
  message(paste0("Number of SNPs in the model:", nrow(snps)))

  if (nrow(snps)<=1) {message(paste0("Skip this gene (nSNPs<=1)!!!!!!!!!")); return()}

  # X <- X[, paste0(snps1$V1,"_",snps1$V4, "_", snps1$V5, "_", snps1$V6) %in% GWAS_snps$V1]

  # define some commonly-used variables
  q<-length(Ys)
  m<-ncol(X)
  N<-nrow(X)
  Yhats_tiss<-vector("list",q)
  Yhats_hom<-vector("list",q)
  Yhats_het<-vector("list", q)

  for(i in 1:q){
    Yhats_het[[i]]<-matrix(NA, ncol=1, nrow=nrow(Ys[[i]]), 
                           dimnames = list(rownames(Ys[[i]]), "pred"))
    Yhats_hom[[i]]<-matrix(NA, ncol=1, nrow=nrow(Ys[[i]]), 
                           dimnames = list(rownames(Ys[[i]]), "pred"))
    Yhats_tiss[[i]]<-matrix(NA, ncol=1, nrow=nrow(Ys[[i]]), 
                            dimnames = list(rownames(Ys[[i]]), "pred"))
  }
  
  if(is.null(covs) & is.null(cov_file_dir)){
    if(verbose){
      message("Continuing without covariates")
    }
    ## The Ys are already residualized
    hom_expr_mat<-matrix(NA, nrow = nrow(X), ncol=q)
    rownames(hom_expr_mat)<-rownames(X)
    for(i in 1:q){
      hom_expr_mat[rownames(Ys[[i]]),i]<-Ys[[i]]
    }
  }else{
    if(verbose){
      message("Residualizing over covariates.")
    }
    
    
    ## If there are covariates, generate residuals:
    hom_expr_mat<-matrix(NA, nrow = nrow(X), ncol=q)
    rownames(hom_expr_mat)<-rownames(X)
    
    for(i in 1:q){
      #### Regress out covariates for each tissue,
      ## If user supplies one directory for covariates, make sure to match
      covidx=match(names(Ys)[i], names(covs))
      tmp<-lm(formula = Ys[[i]] ~ covs[[covidx]][rownames(Ys[[i]]), ])
      ## check the names of residuals, if there's not a perfect overlap
      ### with rownames(Ys[[i]]) they had missing covs and we remove
      if(all(rownames(Ys[[i]]) %in% names(tmp$residuals))){
        hom_expr_mat[rownames(Ys[[i]]),i]<-tmp$residuals
      }else{
        use_inds=rownames(Ys[[i]])[which(rownames(Ys[[i]]) %in% names(tmp$residuals))]
        Ys[[i]]=Ys[[i]][use_inds,,drop=F]
        hom_expr_mat[rownames(Ys[[i]]),i]<-tmp$residuals[rownames(Ys[[i]])]
      } 
    }
    
    # tmp<-lm(formula = Yw ~ covs[[which(names(covs)=="genexp_all_CTs")]][rownames(Yw),])
    # ### with rownames(Yw) they had missing covs and we remove
    # if(all(rownames(Yw) %in% names(tmp$residuals))){
    #   hom_expr_mat[rownames(Yw),q+1]<-tmp$residuals
    # }else{
    #   use_inds=rownames(Yw)[which(rownames(Yw) %in% names(tmp$residuals))]
    #   Yw=Yw[use_inds,,drop=F]
    #   hom_expr_mat[rownames(Yw),q+1]<-tmp$residuals[rownames(Yw)]
    # }    
  }
  
  X_all <- X
  all_missing<-names(rowMeans(hom_expr_mat, na.rm = T)[which(is.nan(rowMeans(hom_expr_mat, na.rm = T)))])
  remove_inds<-which(rownames(hom_expr_mat) %in% all_missing)
  # These people are not in any tissues, we can just remove them from the data
  if(length(remove_inds)>0){
    if(verbose){
      message("These individuals have NA in all Ys or covs: ")
      message(paste0(all_missing, collapse = ","))
    }
    hom_expr_mat<-hom_expr_mat[-remove_inds,,drop=F]
    X<-X[-remove_inds,]
  }
 
 
  if(verbose){
    message("Starting cross-validation")
    message("CONTENT temporary file is ", paste0(out_dir,gene_name, "_content_tmp.bk"))
  }
  if(file.exists(paste0(out_dir,gene_name, "_content_tmp.bk"))){
    system(paste0("rm ", paste0(out_dir,gene_name, "_content_tmp.bk")))
  }
  explanatory=bigstatsr::as_FBM(X, backingfile=paste0(out_dir,gene_name, "_content_tmp"))
  
  
  #fit the model on all of the data for TWAS
  #####################################################
  #####################################################
  ## Fit the homogeneous component
  set.seed(seed)

    hom_fit<-bigstatsr::big_spLinReg(X = explanatory, ind.train = match(rownames(hom_expr_mat), rownames(X)),
                          y.train=rowMeans(x=hom_expr_mat, na.rm = T), alphas=c(content_alpha),K=10,warn=F)
    
    hom_beta_vals<-unlist(summary(hom_fit)$beta[
      which.min(summary(hom_fit)$validation_loss)])[1:ncol(explanatory)]
    hom_int<-unlist(summary(hom_fit)$intercept[
      which.min(summary(hom_fit)$validation_loss)])
    if(sum(is.na(hom_beta_vals))>0){
      hom_beta_vals[which(is.na(hom_beta_vals))]<-0
    }
    if(length(attr(hom_fit, "ind.col"))<dim(explanatory)[2]){
      new_betas<-rep(0, dim(explanatory)[2])
      new_betas[attr(hom_fit, "ind.col")] <- hom_beta_vals[1:length(attr(hom_fit, "ind.col"))]
      hom_beta_vals<-new_betas
    }

    
    names(hom_beta_vals) <- snps$V2
 
  
  ## Fit the heterogeneous components
  het_tiss_betas<-vector("list", q)
  het_tiss_ints<-vector("list", q)
  
  tiss_betas<-vector("list", q)
  tiss_ints<-vector("list", q) 

  for(j in 1:q){

      ### Find out which individuals were not present in other tissues
      nan_names<-names(rowMeans(hom_expr_mat[rownames(Ys[[j]]),-j],
                                na.rm = T)[which(is.nan(rowMeans(hom_expr_mat[rownames(Ys[[j]]),-j], na.rm = T)))])
      ### It may happen that there are individuals in this training set who
      ### do not appear in any other tissues, or they only appear in the test
      ### set of other tissues. If that's the case, remove them
      if(length(nan_names) == 0){
        subset_inds<-intersect(rownames(Ys[[j]]), rownames(hom_expr_mat))
        cur_hom_expr_mat<-hom_expr_mat[subset_inds,,drop=F]
      }else{
        subset_inds<-rownames(Ys[[j]])[!(rownames(Ys[[j]]) %in% nan_names)]
        subset_inds<-intersect(subset_inds, rownames(hom_expr_mat))
        cur_hom_expr_mat<-hom_expr_mat[subset_inds,,drop=F]
      }
      ### First, keep only the train individuals that are in this tissue+other tissues
      
      ## Fit the tissue by tissue approach as well
      tiss_betas[[j]]<-rep(NA, m)
      tiss_ints[[j]]<-NA # 0
      tiss_response=hom_expr_mat[rownames(Ys[[j]]),j]
      if(length(tiss_response) < 15){
        next
      }

       set.seed(j)
      tiss_fit<-bigstatsr::big_spLinReg(X = explanatory, match(rownames(Ys[[j]]), rownames(X)),
                               y.train = tiss_response,K=10, alphas = c(tissue_alpha),warn=F)
       
      tiss_beta_vals<-unlist(summary(tiss_fit)$beta[
          which.min(summary(tiss_fit)$validation_loss)])[1:m]
      tiss_tiss_int<-unlist(summary(tiss_fit)$intercept[
          which.min(summary(tiss_fit)$validation_loss)])
      if(sum(is.na(tiss_beta_vals))>0){
        tiss_beta_vals[which(is.na(tiss_beta_vals))]<-0
      }
      if(length(attr(tiss_fit, "ind.col"))<dim(explanatory)[2]){
        new_betas<-rep(0, dim(explanatory)[2])
        new_betas[attr(tiss_fit, "ind.col")] <- tiss_beta_vals[1:length(attr(tiss_fit, "ind.col"))]
        tiss_beta_vals<-new_betas
      }

      tiss_betas[[j]]<-tiss_beta_vals
      tiss_ints[[j]]<-tiss_tiss_int
        
      names(tiss_betas[[j]]) <- snps$V2



      ### tissue_j expression - mean(all tissues except j expression)
      ##############################################################
      het_response<-cur_hom_expr_mat[, j]-rowMeans(x=cur_hom_expr_mat,  na.rm = T)
   
      if(length(het_response) < 15){
        het_tiss_betas[[j]]<- rep(NA, m)
        het_tiss_ints[[j]]<-NA #0
        next
      }
   
      set.seed(j)
      het_fit<-bigstatsr::big_spLinReg(X = explanatory, ind.train = match(rownames(cur_hom_expr_mat), rownames(X)),
                          y.train = het_response, K=10, alphas = c(content_alpha),warn=F)
  
      het_beta_vals<-unlist(summary(het_fit)$beta[
        which.min(summary(het_fit)$validation_loss)])[1:m]
      het_tiss_int<-unlist(summary(het_fit)$intercept[
        which.min(summary(het_fit)$validation_loss)])
      if(sum(is.na(het_beta_vals))>0){
        het_beta_vals[which(is.na(het_beta_vals))]<-0
      }
      if(length(attr(het_fit, "ind.col"))<dim(explanatory)[2]){
        new_betas<-rep(0, dim(explanatory)[2])
        new_betas[attr(het_fit, "ind.col")] <- het_beta_vals[1:length(attr(het_fit, "ind.col"))]
        het_beta_vals<-new_betas
      }

      het_tiss_betas[[j]]<-het_beta_vals
      het_tiss_ints[[j]]<-het_tiss_int

      names(het_tiss_betas[[j]]) <- snps$V2
   
  }
  
  for(i in 1:q){
    Yhats_hom[[i]][rownames(Ys[[i]]),]<-X[rownames(Ys[[i]]),] %*% hom_beta_vals + hom_int
    Yhats_het[[i]][rownames(Ys[[i]]),]<-X[rownames(Ys[[i]]),] %*% het_tiss_betas[[i]] + het_tiss_ints[[i]]
    Yhats_tiss[[i]][rownames(Ys[[i]]),]<-X[rownames(Ys[[i]]),] %*% tiss_betas[[i]] + tiss_ints[[i]]
  }
  
  names(Yhats_het)<-names(Ys)
  names(Yhats_hom)<-names(Ys)
  names(Yhats_tiss)<-names(Ys)
  
  names(het_tiss_betas) <- names(Ys)
  names(het_tiss_ints) <- names(Ys)
  names(tiss_betas) <- names(Ys)
  names(tiss_ints) <- names(Ys)

  hom_expr_mat <- as.data.frame(hom_expr_mat)
  hom_expr_mat$sampleid <- rownames(hom_expr_mat)
  hom_expr_mat <- hom_expr_mat[,c(ncol(hom_expr_mat), 1:(ncol(hom_expr_mat)-1))]
  CTnames <- substr(names(Ys), 8, nchar(names(Ys)))
  colnames(hom_expr_mat)[2:ncol(hom_expr_mat)] <- CTnames
  X=X_all

  data.table::fwrite(hom_expr_mat, file=paste0(out_dir,gene_name, "_Y_residual"), sep="\t", quote=F)

  save(Yhats_het, Yhats_hom, Yhats_tiss, file = paste0(out_dir,gene_name,"_predictors"))
  save(hom_beta_vals, hom_int, tiss_betas, tiss_ints, het_tiss_betas, het_tiss_ints, 
       snps, X, 
       file = paste0(out_dir,gene_name,"_beta"))
  
  #system(paste0("rm ", paste0(out_dir, gene_name, "_content_tmp.bk")))
  message("Done!")
  
}  





