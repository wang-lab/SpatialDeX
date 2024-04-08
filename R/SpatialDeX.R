#' @title replace_smallest_with_zero
#' @description This function replaces the smallest values in a vector with zero.
#' @param row a numeric vector
#' @param top_cell_type defind the number of cell types to retain based on the proportion of each cell type. Choose an integer from 0 to 9.
#' @return a numeric vector with the smallest values replaced by zero.
#' @importFrom dplyr %>%
#' @importFrom dplyr filter
#' @export
#' @examples
#' row <- c(0.2, 0.1, 0.4, 0.3)
#' top_cell_type <-2
#' replace_smallest_with_zero(row,top_cell_type)

replace_smallest_with_zero <- function(row,top_cell_type) {
  smallest_indices <- order(row)[1:(9 - top_cell_type)]  # Indices of the 5 smallest values
  #smallest_indices <- which(row<0.1)
  row[smallest_indices] <- 0          # Replace with 0
  return(row)
}


#' @title SpatialDeX
#' @name SpatialDeX
#' @description Reference-free cell deconvolution for spatial transcriptome
#' @param raw_exprs raw gene expression matrix, with gene symbols in rows and spot names in columns.
#' @param lambda lambda
#' @param sam.name sample name
#' @param LOW.DR minimal population fractions of genes for smoothing for copykat.
#' @param UP.DR minimal population fractions of genes for segmentation for copykat.
#' @param ngene.chr minimal number of genes per chromosome for cell filtering for copykat.
#' @param KS.cut segmentation parameters, input 0 to 1; larger looser criteria for copykat.
#' @param win.size minimal window sizes for segmentation for copykat.
#' @param top_cell_type defind the number of cell types to retain based on the proportion of each cell type. Choose an integer from 0 to 9.
#' @importFrom dplyr %>%
#' @importFrom dplyr filter
#' @return prediction_result_filter
#' @export


SpatialDeX = function(raw_exprs,
                      lambda,
                      sam.name,
                      LOW.DR,
                      UP.DR,
                      ngene.chr,
                      KS.cut,
                      win.size,
                      top_cell_type=NULL) {
  raw_exprs<-as.matrix(raw_exprs)
  print("start gene marker enrichment")
  gsva_score<-as.data.frame(t(GSVA::gsva(raw_exprs, markerlist)))
  print("start cnv detection")
  norm.cell.names= colnames(normal_count)
  spatial_count<-data.frame(raw_exprs,check.names = FALSE)
  spatial_count$gene<-rownames(spatial_count)
  normal_count_tmp<-normal_count
  normal_count_tmp$gene<-rownames(normal_count_tmp)
  all_count<-merge(spatial_count,normal_count_tmp,by="gene")
  rownames(all_count)<-all_count$gene
  all_count<-all_count[,-1]
  id.type="S"
  load(system.file("extdata", "full.anno.RData", package = "SpatialDeX"))
  load(system.file("extdata", "cyclegenes.RData", package = "SpatialDeX"))
  load(system.file("extdata", "DNA.hg20.RData", package = "SpatialDeX"))
  cnv_data<-copykat::copykat(rawmat=all_count,id.type=id.type, LOW.DR=LOW.DR, UP.DR=UP.DR,ngene.chr=ngene.chr,sam.name=sam.name,
                             KS.cut=KS.cut,win.size=win.size, norm.cell.names= norm.cell.names,plot.genes="FLASE")

  cnv_result<-cnv_data$prediction%>%filter(copykat.pred!="not.defined")
  gsva_score$cell.names<-rownames(gsva_score)
  #cnv_result$cell.names<-rownames(cnv_result)
  combined_data<-merge(cnv_result,gsva_score,by="cell.names")
  combined_data$cnv_info<-ifelse(combined_data$copykat.pred=="aneuploid",1,0)
  rownames(combined_data)<-combined_data$cell.names
  combined_data<-combined_data[,-1]
  feature_names<-intersect(colnames(combined_data),colnames(weighted_matrix))
  weighted_feature<-weighted_matrix[,feature_names]
  data_feature<-combined_data[,feature_names]
  if( identical(colnames(weighted_feature),colnames(data_feature))==TRUE){

    data_feature[is.na(data_feature)]<-0
    data_feature<-apply(data_feature,c(1,2),as.numeric)
    data_feature<-apply(data_feature,2,function(x) (x-min(x))/(max(x)-min(x)))
    weighted_feature<-as.matrix(weighted_feature)
    prediction_result<-data_feature%*%t(weighted_feature)%*%solve((weighted_feature %*% t(weighted_feature)) +lambda * diag(nrow(weighted_feature)))
    prediction_result<-as.data.frame(prediction_result)
    prediction_result[prediction_result<0]=0
    proportion_max <- apply(prediction_result, 1, max)
    rm_spot_1<-which(proportion_max <0.15)
    rm_spot_2<-which(rowSums(prediction_result>0.3)>2)
    rm_spot<-c(rm_spot_1,rm_spot_2)
    if(length(rm_spot)>0){
      prediction_result_rm<-data.frame(prediction_result[!rownames(prediction_result)%in% names(rm_spot),])
      uniden_spot <- data.frame(matrix(rep("unidentified", length(names(rm_spot)) * 9), nrow = length(names(rm_spot))))
      rownames(uniden_spot) <- names(rm_spot)
      colnames(uniden_spot)<-colnames(prediction_result)
    }else{
      prediction_result_rm<-prediction_result
      uniden_spot <-data.frame()
    }


    if (!is.null(top_cell_type)) {
      prediction_result_filter <- as.data.frame(t(apply(prediction_result_rm, 1, replace_smallest_with_zero, top_cell_type=top_cell_type)))
    } else {
      prediction_result_filter <- prediction_result_rm
    }

    prediction_result_filter<-t(apply(prediction_result_filter,1,function(x) x/sum(x)))
    prediction_result_filter[is.nan(prediction_result_filter)] <- 0
    prediction_result_filter<-rbind(prediction_result_filter,uniden_spot)

    return(prediction_result_filter)
  }}

