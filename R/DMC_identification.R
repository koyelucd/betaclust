#' @title The DMC identification function
#' @export
#' @description A function to identify the most differentially
#' methylated clusters from \eqn{K} clusters.
#'
#' @details This function selects the most diffentially methylated clusters
#' based on AUC and WD metric calculated and returns the CpG sites
#' belonging to those clusters.
#'
#'

#' @seealso \code{\link{beta_kr}}
#' @seealso \code{\link{pca.methylation.data}}
#' @seealso \code{\link{plot.betaclust}}
#' @seealso \code{\link{summary.betaclust}}
#' @seealso \code{\link{betaclust}}
#'
#' @param object A betaclust object
#' @param data A dataframe of dimension \eqn{C \times NR} containing methylation
#' values for \eqn{C} CpG sites from \eqn{R} samples collected from \eqn{N}
#' patients which was passed as an argument to the
#' \code{\link[betaclust:betaclust]{betaclust}} function.
#' @param CpG_site_list The IlmnID of all the CpG sites analysed by
#' \code{\link[betaclust:betaclust]{betaclust}} function.
#' @param threshold The threshold value/s for selecting the most differentially
#'  methylated clusters, default= 0.65
#' @param metric The metric (AUC or WD selected). default="AUC"
#'
#' @return The function returns a dataframe of CpG sites and
#'  methylation values identified to belong to the most
#'   differentially methylated clusters
#'
#' @examples
#' \donttest{
#' my.seed <- 190
#' M <- 3
#' N <- 4
#' R <- 2
#' data_output <- betaclust(pca.methylation.data[1:30,2:9], M, N, R,
#'             model_names = "K.R",
#'             parallel_process = FALSE, seed = my.seed)
#' dmc_df <-DMC_identification(data_output,pca.methylation.data[1:30,2:9],
#' pca.methylation.data[1:30,1],
#'  threshold = 0.65, metric = "AUC")
#'
#'}
#' @importFrom foreach %dopar%
#' @importFrom stats C
#' @importFrom utils flush.console
#' @importFrom utils txtProgressBar

DMC_identification<-function(object, data,CpG_site_list, threshold = 0.65, metric = "AUC")
{
  AUC=object$optimal_model_results$DM$AUC
  WD=object$optimal_model_results$DM$WD
  dmc_cluster=vector()
  if(metric=="AUC")
  {
    for(j in 1:object$K)
    {
      dmc_cluster[j]=ifelse(AUC[j]>=threshold,1,0)
    }
  }else{
    for(j in 1:object$K)
    {
      dmc_cluster[j]=ifelse(WD[j]>=threshold,1,0)
    }
  }
  dmc_cluster_index=which(dmc_cluster==1)
  classification=object$optimal_model_results$classification
  df=as.data.frame(cbind(CpG_site_list,data,classification))
  colnames(df)[1]="IlmnID"
  dmc_df=df[df$classification %in% dmc_cluster_index,]

  return(dmc_df)

}
