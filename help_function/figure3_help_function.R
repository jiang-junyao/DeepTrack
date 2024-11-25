gene_enrich <- function(gene,org.db,enrich.db,use_internal_data = TRUE,
                        organism = NULL,pvalueCutoff = 0.05){
  library(clusterProfiler)
  gene1 <- clusterProfiler::bitr(gene, fromType = "SYMBOL",
                                 toType = c("SYMBOL", "ENTREZID"),
                                 OrgDb = org.db)
  if (enrich.db =='KEGG') {
    k1 <- clusterProfiler::enrichKEGG(gene = gene1$ENTREZID,
                                      pvalueCutoff = pvalueCutoff
                                      ,use_internal_data = use_internal_data,
                                      minGSSize = 2)
  }else if(enrich.db =='GO'){
    k1 = clusterProfiler::enrichGO(gene = gene1$ENTREZID,
                                   OrgDb = org.db,
                                   keyType = "ENTREZID",
                                   ont = "BP",
                                   pvalueCutoff = pvalueCutoff,minGSSize = 2)
  }
  return(k1)
}

extract_expressed_TFs <- function(seurat_object,TFs,cells_quantile = 0.05){
  TFs <- TFs[TFs%in%rownames(seurat_object)]
  matrix_tf <- seurat_object@assays$RNA@counts[TFs,]
  if (cells_quantile==0) {
    TfExp <- matrix_tf[rowSums(as.matrix(matrix_tf))>0,]
  }else{
    quantile_exp <- ncol(as.matrix(matrix_tf))/(1/cells_quantile)
    TfExp <- matrix_tf[ncol(as.matrix(matrix_tf))-rowSums(as.matrix(matrix_tf==0))>quantile_exp,]}
  return(TfExp)
}
diffgenetest_pseudotime <- function(monocle_object){
  mo <- detectGenes(monocle_object, min_expr = 1)
  #mo <- estimateSizeFactors(mo)
  diff1 <- monocle::differentialGeneTest(mo,
                                         fullModelFormulaStr = "~Pseudotime",
                                         relative_expr = TRUE,cores=20
  )
  ed <- c()
  for (i in diff1$num_cells_expressed) {
    a <- log(i) * 0.95 - log(i) * 0.05
    ed <- c(ed, a)
  }
  diff1$expression_difference <- ed
  return(diff1)
}
sparse.cor <- function(x){
  n <- nrow(x)
  cMeans <- colMeans(x)
  cSums <- colSums(x)
  # Calculate the population covariance matrix.
  # There's no need to divide by (n-1) as the std. dev is also calculated the same way.
  # The code is optimized to minize use of memory and expensive operations
  covmat <- tcrossprod(cMeans, (-2*cSums+n*cMeans))
  crossp <- as.matrix(crossprod(x))
  covmat <- covmat+crossp
  sdvec <- sqrt(diag(covmat)) # standard deviations of columns
  covmat/crossprod(t(sdvec)) # correlation matrix
}
min_max_scale <- function(predicted_values){
  v1=(predicted_values - min(predicted_values)) / (max(predicted_values) - min(predicted_values))
  return(v1)
}
get_pseudotime_de <- function(data, reverse = FALSE,gene.use = NULL) {
  pd <- data.frame(Pseudotime=1:ncol(data),sample='meta',row.names = colnames(data))
  pd <- new("AnnotatedDataFrame", data = pd)
  fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
  fd <- new("AnnotatedDataFrame", data = fData)
  data1=as(as.matrix(data), "sparseMatrix")
  monocle_cds <- monocle::newCellDataSet(data1,
                                         phenoData = pd,
                                         featureData = fd,
                                         lowerDetectionLimit = 0.1,
                                         expressionFamily = VGAM::negbinomial.size()
  )
  cds <- BiocGenerics::estimateSizeFactors(monocle_cds)
  deg=diffgenetest_pseudotime(cds)
  return(deg)
}
get_pseudotime_de_seurat <- function(data, reverse = FALSE,gene.use = NULL) {
  pd <- data.frame(Pseudotime=data$scale_pseudotime,
                   sample='meta',row.names = colnames(data))
  pd <- new("AnnotatedDataFrame", data = pd)
  fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
  fd <- new("AnnotatedDataFrame", data = fData)
  #data1=as(data@assays$RNA@counts, "sparseMatrix")
  monocle_cds <- monocle::newCellDataSet(data@assays$RNA@data,
                                         phenoData = pd,
                                         featureData = fd,
                                         lowerDetectionLimit = 0.1,
                                         expressionFamily = VGAM::negbinomial.size()
  )
  cds <- BiocGenerics::estimateSizeFactors(monocle_cds)
  deg=diffgenetest_pseudotime(cds)
  return(deg)
}
gene_cor = function(candidate_gene,obj,col1){
  library(foreach)
  library(doParallel)
  registerDoParallel(cores=10)
  pval = c()
  corv = c()
  results <- foreach(i = candidate_gene, .combine=rbind) %dopar% {
    cor_test <- cor.test(obj@assays$RNA@data[i,], obj@meta.data[,col1])
    cbind(gene=i, pval=cor_test$p.value, corv=cor_test$estimate)
  }
  stopImplicitCluster()
  cor_df <- as.data.frame(results)
  colnames(cor_df) <- c("candidate_gene", "pval", "corv")
  return(cor_df)
}

smooth_pseudotime<-function(exp,metadata,feature_use,bin_num=100,group=NULL,
                            order1){
  library(tidyverse)
  if (is.null(group)) {
    smooth_exp = get_SmoothByBin_PseudotimeExp(exp,metadata,bin_num)
    y=1:bin_num
  }else{
    smooth_list = list()
    for (i in order1) {
      meta_use1 = metadata[metadata[,group]==i,]
      exp_use = exp[,rownames(meta_use1)]
      smooth_exp = get_SmoothByBin_PseudotimeExp(exp_use,meta_use1,bin_num)
      colnames(smooth_exp) = paste0(i,1:ncol(smooth_exp))
      smooth_list[[i]] = as.data.frame(smooth_exp)
    }
    smooth_exp = do.call(bind_cols,smooth_list)
    y = 1:(bin_num*length(unique(metadata[,group])))
  }
  print(dim(smooth_exp))
  smooth_exp = as.matrix(smooth_exp)
  fitted_list = list()
  all_exp_diff = c()
  for (i in feature_use) {
    predicted_values=loess(smooth_exp[i,]~y)$fitted
    exp_diff = quantile(smooth_exp[i,],0.95)-quantile(smooth_exp[i,],0.05)
    all_exp_diff = c(all_exp_diff,exp_diff)
    # if (exp_diff>0.1) {
    #   scaled_values <- (predicted_values - min(predicted_values)) / (max(predicted_values) - min(predicted_values))
    #   fitted_list[[i]]=predicted_values
    # }
    fitted_list[[i]]=predicted_values
  }
  names(all_exp_diff) = feature_use
  fitted_mt = t(as.data.frame(fitted_list))
  rownames(fitted_mt) = feature_use
  smooth_list = list(smooth_exp[feature_use,],fitted_mt,all_exp_diff)
  return(smooth_list)
}
scale_vector <- function(val){
  scaled_values <- (val - min(val)) / (max(val) - min(val))
  return(scaled_values)
}
scale_mt <- function(mt){
  scale_list = list()
  for (i in 1:nrow(mt)) {
    value_use = mt[i,]
    scaled_values <- (value_use - min(value_use)) / (max(value_use) - min(value_use))
    scale_list[[i]] = scaled_values
  }
  scale_df = t(data.frame(scale_list))
  rownames(scale_df) = rownames(mt)
  colnames(scale_df) = colnames(mt)
  return(scale_df)
}
smoothByBin <- function(Exp1, Pseudotime1, PseudotimeRange1 = NULL,
                        SmoothLength1 = 100, ByBin1 = c("Equal.Pseudotime", "Equal.Cells")) {
  if (ncol(Exp1) != nrow(Pseudotime1)) {
    stop("The length of pseudotime is not equal to the number of cells")
  } else if (!is.element("Pseudotime", colnames(Pseudotime1))) {
    stop("No pseudotime inforamtion in variable Pseudotime1")
  }
  
  if (ByBin1[1] == "Equal.Pseudotime") {
    if (is.null(PseudotimeRange1)) {
      PseudotimeBin1 <- seq(min(Pseudotime1$Pseudotime), max(Pseudotime1$Pseudotime),
                            length.out = SmoothLength1 + 1)
    } else {
      PseudotimeBin1 <- seq(PseudotimeRange1[1], PseudotimeRange1[2],
                            length.out = SmoothLength1 + 1)
    }
  } else {
    Pseudotime1 <- Pseudotime1[order(Pseudotime1$Pseudotime), ]
    Bin1 <- ceiling(nrow(Pseudotime1) / SmoothLength1)
    PseudotimeBin1 <- seq(1, nrow(Pseudotime1), by = Bin1)
    PseudotimeBin1[length(PseudotimeBin1) + 1] <- nrow(Pseudotime1)
  }
  
  Exp2 <- array(0, dim = c(nrow(Exp1), SmoothLength1))
  for (i in 1:(length(PseudotimeBin1) - 1)) {
    if (ByBin1[1] == "Equal.Pseudotime") {
      if (i == length(PseudotimeBin1) - 1) {
        Cells1 <- Pseudotime1[Pseudotime1$Pseudotime >= PseudotimeBin1[i] &
                                Pseudotime1$Pseudotime <= PseudotimeBin1[i + 1], ]
      } else {
        Cells1 <- Pseudotime1[Pseudotime1$Pseudotime >= PseudotimeBin1[i] &
                                Pseudotime1$Pseudotime < PseudotimeBin1[i + 1], ]
      }
    } else {
      if (i == length(PseudotimeBin1) - 1) {
        Cells1 <- Pseudotime1[PseudotimeBin1[i]:PseudotimeBin1[i + 1], ]
      } else {
        Cells1 <- Pseudotime1[PseudotimeBin1[i]:(PseudotimeBin1[i + 1] - 1), ]
      }
    }
    
    if (nrow(Cells1) > 1) {
      Exp2[, i] <- rowMeans(Exp1[, match(rownames(Cells1), colnames(Exp1))])
    } else if (nrow(Cells1) == 1) {
      Exp2[, i] <- Exp1[, match(rownames(Cells1), colnames(Exp1))]
    }
  }
  rownames(Exp2) <- rownames(Exp1)
  
  return(Exp2)
}




smoothByState <- function(seurat_obj, Pseudotime1, each_state_bin1){
  ### order cells
  Pseudotime1 <- Pseudotime1[order(Pseudotime1$Pseudotime),]
  Pseudotime1$State <- as.character(Pseudotime1$State)
  seurat_exp <- as.matrix(seurat_obj@assays$RNA@data)
  seurat_exp <- seurat_exp[,rownames(Pseudotime1)]
  ### order states
  order_state <- as.character(Pseudotime1$State[!duplicated(Pseudotime1$State)])
  ### divide cells in each state to bins
  each_state_bin <- each_state_bin1
  for (j in order_state) {
    
    Pseudotime2 <- Pseudotime1[Pseudotime1$State==j,]
    seurat_exp2_state <- seurat_exp[,rownames(Pseudotime2)]
    PseudotimeBin1 <- round(seq(1,nrow(Pseudotime2),length.out=each_state_bin+1))
    Exp2 <- array(0, dim = c(nrow(seurat_exp), each_state_bin))
    for (i in 1:(length(PseudotimeBin1) - 1)) {
      Exp2[, i] <- rowMeans(seurat_exp2_state[, PseudotimeBin1[i]:PseudotimeBin1[i+1]])
    }
    if (j==1) {
      Exp3 = Exp2
    }else{Exp3 <-cbind(Exp3,Exp2)}
  }
  rownames(Exp3) <- rownames(seurat_exp)
  return(Exp3)
}


get_SmoothByBin_PseudotimeExp <- function(mt,metadata,Bin = 50) {
  
  
  
  TotalBin1 <- Bin
  Exp1 <- smoothByBin(mt, 
                      metadata,
                      SmoothLength1 = TotalBin1, ByBin1 = c("Equal.Pseudotime"))
  
  
  return(Exp1)
  
}

gene_cor_df = function(mt){
  library(foreach)
  library(doParallel)
  registerDoParallel(cores=10)
  pval = c()
  corv = c()
  results <- foreach(i = rownames(mt), .combine=rbind) %dopar% {
    cor_test <- cor.test(mt[i,], 1:ncol(mt))
    cbind(gene=i, pval=cor_test$p.value, corv=cor_test$estimate)
  }
  stopImplicitCluster()
  cor_df <- as.data.frame(results)
  colnames(cor_df) <- c("candidate_gene", "pval", "corv")
  cor_df[,3] = as.numeric(cor_df[,3])
  cor_df[,2] = as.numeric(cor_df[,2])
  return(cor_df)
}