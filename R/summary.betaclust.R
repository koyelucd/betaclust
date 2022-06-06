#' @title Summary statistics of betaclust output
#' @description Calculates and prints the summary statistics of the optimal model selected for printing.
#' @export
#' @param x betaclust object
#' @return An object of class "summary.betaclust".
#' The object returns the following list of values:
#' \itemize{
#' \item CpG_sites - The number of CpG sites analysed using the beta mixture models.
#' \item patients - The number of patients analysed using the beta mixture models.
#' \item samples - The numder of samples analysed using the beta mixture models.
#' \item cluster_count - The number of groups, the data is clustered into.
#' \item modelName - The optimal model selected.
#' \item loglik - The log-likelihood value for the selected optimal model.
#' \item Information_criterion - The information criterion used to select the optimal model.
#' \item ic_output - This stores the information criterion value calculated for each model.
#' \item classification - The total number of CpG sites identified in each cluster. }

summary.betaclust<-function(object)
{
  title <- paste("Multivariate Beta mixture model fitted by EM algorithm")
  ic_value=min(object$ic_output)
  clust_count=length(object$best_model$cluster_count)
  loglik=object$best_model$llk[length(object$best_model$llk)]
  #classification<-as.factor(object$best_model$data[,ncol(object$best_model$data)])
  clustering<-as.factor(object$best_model$data[,ncol(object$best_model$data)])
  classification<-table(clustering)
  obj <- list(title = title,
              CpG_sites=object$CpG_sites,
              patients=object$patients,
              samples=object$samples,
              cluster_count = clust_count,
              modelName = object$optimal_model,
              loglik = loglik,
              information_criterion=object$information_criterion,
              ic_output=ic_value,
              classification=classification)
  class(obj) <- "summary.betaclust"
  return(obj)
}
