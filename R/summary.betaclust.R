#' ICL to compare the methods
#' @export
#' @param x betaclust object
#' @return obj

summary.betaclust<-function(object)
{
  K<-object$K
  title <- paste("Multivariate Beta mixture model fitted by EM algorithm")
  obj <- list(title = title,
              C = nrow(object$best_model$data),
              d = ncol(object$best_model$data),
              cluster_count = length(object$best_model$cluster_count),
              modelName = object$optimal_model,
              loglik = object$best_model$llk[length(object$best_model$llk)],
              information_criterion=object$model_selection,ic_output=min(object$ic_op))
  class(obj) <- "summary.betaclust"
  return(obj)
}
