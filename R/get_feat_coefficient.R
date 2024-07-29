#' @title getsolvedmatrix
#' @description Calculates the solution matrix using ridge regression.
#' @param Xm Matrix of predictors.
#' @param Ym Matrix of responses.
#' @param LAMBDA Regularization parameter,default = 0.01.
#' @return Matrix Amat_g, the solution matrix.
#' @export
#' @examples
#' Xm <- matrix(rnorm(100), nrow=10, ncol=10)
#' Ym <- matrix(rnorm(10), nrow=10, ncol=1)
#' getsolvedmatrix(Xm, Ym)

getsolvedmatrix <- function(Xm, Ym, LAMBDA = 0.01) {
  # Amat=solve(Xmat,Ymat) # solve AX=B, or Xmat * A =Ymat
  TMmat_g = (t(Xm) %*% Xm) + LAMBDA * diag(ncol(Xm))

  Amat_g = solve(TMmat_g) %*% t(Xm) %*% Ym
  return(Amat_g)
}

#' @title getsolvedmatrix_with_permutation_cell_label
#' @description Computes a solution matrix using ridge regression with the option to permute cell labels.
#' @param Xm Matrix of predictors.
#' @param Ym Matrix of responses.
#' @param LAMBDA Regularization parameter,default = 0.01.
#' @param n_permutation Number of permutations to perform, default = 10,000.
#' @return A list containing two elements:
#'   \item{Amat_ret}{Matrix of coefficients from the original regression.}
#'   \item{Amat_ret_higher}{Matrix indicating the proportion of permutations where the absolute value of
#'   coefficients exceeds those in Amat_ret.}
#' @export
#' @examples
#' Xm <- matrix(rnorm(100), nrow=10, ncol=10)
#' rownames(Xm) <- paste0("cell_", 1:nrow(Xm))
#' Ym <- matrix(rnorm(10), nrow=10, ncol=2)
#' rownames(Ym) <- paste0("cell_", 1:nrow(Ym))
#' getsolvedmatrix_with_permutation_cell_label(Xm, Ym,n_permutation = 10)

getsolvedmatrix_with_permutation_cell_label <- function(Xm, Ym, LAMBDA = 0.01, n_permutation = 10000) {
  #set.seed(000)
  Amat_ret = getsolvedmatrix(Xm, Ym, LAMBDA = LAMBDA)
  Amat_ret_higher = matrix(rep(0, ncol(Amat_ret) * nrow(Amat_ret)), nrow = nrow(Amat_ret))
  rownames(Amat_ret_higher) = rownames(Amat_ret)
  colnames(Amat_ret_higher) = colnames(Amat_ret)
  # permute N times randomly shuffle cell labels
  for (npm in 1:n_permutation) {
    if (npm%%5000 == 0) {
      message(paste("Permutation:", npm, "/", n_permutation, "..."))
    }
    cells_shu = sample(rownames(Ym), nrow(Ym))
    Xm_s = Xm[cells_shu, ]
    Ym_s = Ym  # [cells_shu,]
    rownames(Ym_s) = cells_shu
    Amat_random = getsolvedmatrix(Xm_s, Ym_s, LAMBDA = LAMBDA)

    Amat_ret_higher = Amat_ret_higher + (abs(Amat_random) > abs(Amat_ret)) * 1
    # browser()
  }
  Amat_ret_higher = Amat_ret_higher/n_permutation
  return(list(Amat_ret, Amat_ret_higher))
}


#' @title Get the Weighted Feature Coefficient by Users
#' @description This function allows users to obtain their own weighted feature coefficient matrix.
#' @param cell_identity A numeric matrix indicating cell identities,with spots in rows and cell identities in columns.
#' @param raw_exprs Raw gene expression matrix, with gene symbols in rows and spot names in columns.
#' @param n_permutation An integer indicating the number of permutations for randomly shuffling the cell
#'     labels of the input matrix to create a null distribution. Default is 10,000.
#' @param LAMBDA A numeric value for the regularization parameter. Default is 0.01.
#' @param sam.name Sample name
#' @param LOW.DR Minimal population fractions of genes for smoothing for copykat. Default is 0.03.
#' @param UP.DR Minimal population fractions of genes for segmentation for copykat. Default is 0.2.
#' @param ngene.chr Minimal number of genes per chromosome for cell filtering for copykat. Default is 5.
#' @param KS.cut Segmentation parameters, input 0 to 1; larger looser criteria for copykat. Default is 0.1.
#' @param win.size Minimal window sizes for segmentation for copykat. Default is 25.
#' @param pvalue Significance threshold for weighted feature coefficients. Calculated by permutation.
#' @param seed Seed for random number generation for reproducibility.
#' @return A matrix of weighted feature coefficients.
#' @importFrom dplyr %>%
#' @importFrom dplyr filter
#' @export
#' @examples
#' data(example_raw_count)
#' raw_exprs <- example_raw_count[,1:10]
#' cell_identity <- matrix(0, nrow = 10, ncol = 9)
#' rownames(cell_identity) <- colnames(raw_exprs)
#' colnames(cell_identity) <- paste0("cell_type", 1:9)
#' for (i in 1:10) {
#'  n_ones <- sample(0:3, 1)
#'   if (n_ones > 0) {
#'     ones_positions <- sample(1:9, n_ones)
#'     cell_identity[i, ones_positions] <- 1
#'   }
#' }
#' LAMBDA <- 0.01
#' n_permutation <- 10000
#' pvalue <- 0.05
#' sam.name <- "sample1"
#' LOW.DR <- 0.03
#' UP.DR <- 0.2
#' ngene.chr <- 10
#' KS.cut <- 0.1
#' win.size <- 25
#' seed <- 1234
#' get_feat_coefficient(cell_identity,
#'                      raw_exprs,
#'                      LAMBDA,
#'                      n_permutation,
#'                      pvalue,
#'                      sam.name,
#'                      LOW.DR,
#'                      UP.DR,
#'                      ngene.chr,
#'                      KS.cut,
#'                      win.size,
#'                      seed)

get_feat_coefficient<- function(cell_identity,
                                raw_exprs=raw_exprs,
                                LAMBDA=0.01,
                                n_permutation=10000,
                                pvalue=pvalue,
                                sam.name=sam.name,
                                LOW.DR=0.03,
                                UP.DR=0.2,
                                ngene.chr=5,
                                KS.cut=0.1,
                                win.size=25,
                                seed) {
  set.seed(seed)
  raw_exprs<-as.matrix(raw_exprs)
  message("start gene marker enrichment")
  gsva_score<-as.data.frame(t(GSVA::gsva(raw_exprs, markerlist)))
  message("start cnv detection")
  norm.cell.names<- colnames(normal_count)
  spatial_count<-data.frame(raw_exprs,check.names = FALSE)
  spatial_count$gene<-rownames(spatial_count)
  normal_count_tmp<-normal_count
  normal_count_tmp$gene<-rownames(normal_count_tmp)
  all_count<-merge(spatial_count,normal_count_tmp,by="gene")
  rownames(all_count)<-all_count$gene
  all_count<-all_count[,-1]
  id.type<-"S"
  load(system.file("extdata", "full.anno.RData", package = "SpatialDeX"))
  load(system.file("extdata", "cyclegenes.RData", package = "SpatialDeX"))
  load(system.file("extdata", "DNA.hg20.RData", package = "SpatialDeX"))
  cnv_data <- copykat::copykat(
    rawmat = all_count,
    id.type = id.type,
    LOW.DR = LOW.DR,
    UP.DR = UP.DR,
    ngene.chr = ngene.chr,
    sam.name = sam.name,
    KS.cut = KS.cut,
    win.size = win.size,
    norm.cell.names = norm.cell.names,
    plot.genes = "FALSE"
  )
  copykat.pred <- "copykat.pred"
  cnv_result<-cnv_data$prediction%>%filter(copykat.pred!="not.defined")
  gsva_score$cell.names<-rownames(gsva_score)
  #cnv_result$cell.names<-rownames(cnv_result)
  combined_data<-merge(cnv_result,gsva_score,by="cell.names")
  combined_data$cnv_info<-ifelse(combined_data$copykat.pred=="aneuploid",1,0)
  rownames(combined_data)<-combined_data$cell.names
  combined_data<-combined_data[,-1]
  feature_names<-intersect(colnames(combined_data),c(names(markerlist),"cnv_info"))

  data_feature<-combined_data[,feature_names]
  #if( identical(c(names(markerlist),"cnv_info"),colnames(data_feature))==TRUE){

  data_feature[is.na(data_feature)]<-0
  data_feature<-apply(data_feature,c(1,2),as.numeric)
  data_feature[,-dim(data_feature)[2]]<-apply(data_feature[,-dim(data_feature)[2]],2,function(x) (x-min(x))/(max(x)-min(x)))

  Amat_pm_lst <- getsolvedmatrix_with_permutation_cell_label(cell_identity, data_feature, LAMBDA = LAMBDA, n_permutation = n_permutation)
  Amat_s <- Amat_pm_lst[[1]]
  Amat_pval <- Amat_pm_lst[[2]]

  Amat_s <- t(scale(t(Amat_s), center = TRUE, scale = TRUE))
  Amat_s[is.na(Amat_s)] <- 0

  Amat_s2 <- matrix(rep(0, nrow(Amat_s) * ncol(Amat_s)), nrow = nrow(Amat_s))
  rownames(Amat_s2) <- rownames(Amat_s)
  colnames(Amat_s2) <- colnames(Amat_s)

  for (rn in 1:nrow(Amat_s)) {
    for (cn in 1:ncol(Amat_s)) {
      if (Amat_pval[rn, cn] < pvalue & Amat_s[rn, cn] > 0) {
        Amat_s2[rn, cn] <- Amat_s[rn, cn]
      } else {
        Amat_s2[rn, cn] <- Amat_s2[rn, cn]
      }
    }
  }

  return(Amat_s2)
}
