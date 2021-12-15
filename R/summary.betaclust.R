#' ICL to compare the methods
#' @export
#' @param x betaclust object
#' @return obj

summary.betaclust<-function(object)
{
  #K<-object$K
  #C=nrow(object$best_model$data)
  #d=ncol(object$best_model$data)-1
  title <- paste("Multivariate Beta mixture model fitted by EM algorithm")
  ic_value=min(object$ic_output)
  clust_count=length(object$best_model$cluster_count)
  loglik=object$best_model$llk[length(object$best_model$llk)]
  classification<-as.factor(object$best_model$data[,ncol(object$best_model$data)])
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
