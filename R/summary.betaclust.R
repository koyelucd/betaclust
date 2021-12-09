#' ICL to compare the methods
#' @export
#' @param x betaclust object
#' @return obj

summary.betaclust<-function(object)
{
  K<-object$K
  title <- paste("Multivariate Beta mixture model fitted by EM algorithm")
  obj <- list(title = title,
              C = nrow(object$final_output$data),
              d = ncol(object$final_output$data),
              cluster_count = length(object$final_output$cluster_count),
              modelName = object$optimal_model,
              loglik = object$final_output$llk[length(object$final_output$llk)],
              information_criterion=object$model_selection,ic_output=min(object$ic_op))
  class(obj) <- "summary.betaclust"
  return(obj)
}
