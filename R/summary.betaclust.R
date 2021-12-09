#' ICL to compare the methods
#' @export
#' @param x betaclust object
#' @return obj

summary.betaclust<-function(object)
{
  #K<-object$K
  C=nrow(object$best_model$data)
  d=ncol(object$best_model$data)-1
  title <- paste("Multivariate Beta mixture model fitted by EM algorithm")
  ic_value=min(object$ic_output)
  clust_count=length(object$best_model$cluster_count)
  loglik=object$best_model$llk[length(object$best_model$llk)]
  obj <- list(title = title,
              C = C,
              d = d,
              cluster_count = clust_count,
              modelName = object$optimal_model,
              loglik = loglik,
              information_criterion=object$information_criterion,
              ic_output=ic_value)
  class(obj) <- "summary.betaclust"
  return(obj)
}
