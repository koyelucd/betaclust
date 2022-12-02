#' @title Summarizing the beta mixture model fits
#' @description Summary method for a \code{\link[betaclust:betaclust]{betaclust}} object containing the results under the optimal model selected.
#' @rdname summary.betaclust
#' @export
#' @param object A \code{\link[betaclust:betaclust]{betaclust}} object.
#' @param ... Further arguments passed to or from other methods.
#' @return An object of class \code{\link[betaclust:summary.betaclust]{summary.betaclust}} which contains the following list:
#' \itemize{
#' \item C - The number of CpG sites analysed using the beta mixture models.
#' \item N - The number of patients analysed using the beta mixture models.
#' \item R - The number of samples analysed using the beta mixture models.
#' \item K - The number of methylation states in R DNA samples.
#' \item modelName - The optimal model selected.
#' \item loglik - The log-likelihood value for the selected optimal model.
#' \item information_criterion - The information criterion used to select the optimal model.
#' \item ic_output - This stores the information criterion value calculated for each model.
#' \item classification - The total number of CpG sites in each cluster.
#' \item prop_data - The proportion of CpG sites in each cluster.}
#' @examples
#' \donttest{
#' my.seed <- 190
#' M <- 3
#' N <- 4
#' R <- 2
#' data_output <- betaclust(pca.methylation.data[1:30,2:9], M, N, R,
#'             model_names=c("K..","KN.","K.R"), model_selection="BIC",
#'             parallel_process = FALSE, seed=my.seed)
#' summary(data_output)}
#' @seealso \code{\link{betaclust}}

summary.betaclust<-function(object,...)
{
  title <- paste("Multivariate beta mixture model fitted by EM algorithm")
  ic_value=min(object$ic_output)
  clust_count=length(object$optimal_model_results$cluster_size)
  loglik=object$optimal_model_results$llk[length(object$optimal_model_results$llk)]
  #classification<-as.factor(object$optimal_model_results$data[,ncol(object$optimal_model_results$data)])
  #clustering<-as.factor(object$optimal_model_results$data[,ncol(object$optimal_model_results$data)])
 # classification<-table(as.factor(object$optimal_model_results$data[,ncol(object$optimal_model_results$data)]))
  classification<-table(as.factor(object$optimal_model_results$classification))
  prop_data<-round((as.numeric(object$optimal_model_results$cluster_size)/object$C),3)
  obj <- list(title = title,
              C=object$C,
              N=object$N,
              R=object$R,
              K = clust_count,
              modelName = object$optimal_model,
              loglik = loglik,
              information_criterion=object$information_criterion,
              ic_output=ic_value,
              classification=classification,
              prop_data=prop_data)
  class(obj) <- "summary.betaclust"
  print(obj)
  invisible(obj)
}

