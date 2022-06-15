#' @title Summarizing the beta mixture model fits
#' @description Summary method for a "betaclust" object containing the results under the optimal model selected.
#' @export
#' @param object A betaclust object.
#' @return An object of class "summary.betaclust which contains the following list of values:
#' \itemize{
#' \item C - the number of CpG sites analysed using the beta mixture models.
#' \item N - the number of patients analysed using the beta mixture models.
#' \item R - the numder of samples analysed using the beta mixture models.
#' \item K - the number of methylation profiles identified.
#' \item modelName - The optimal model selected.
#' \item loglik - The log-likelihood value for the selected optimal model.
#' \item Information_criterion - The information criterion used to select the optimal model.
#' \item ic_output - This stores the information criterion value calculated for each model.
#' \item classification - The total number of CpG sites identified in each cluster. }
#' @examples
#' \dontrun{
#' data_output=betaclust(pca.methylation.data[,2:9],K,patients,samples,
#'             model_names=c("C..","CN.","C.R"),model_selection="BIC",seed=my.seed)
#' summary(data_output)}
#' @seealso \code{\link{betaclust}}

summary.betaclust<-function(object)
{
  title <- paste("Multivariate Beta mixture model fitted by EM algorithm")
  ic_value=min(object$ic_output)
  clust_count=length(object$optimal_model_results$cluster_size)
  loglik=object$optimal_model_results$llk[length(object$optimal_model_results$llk)]
  #classification<-as.factor(object$optimal_model_results$data[,ncol(object$optimal_model_results$data)])
  clustering<-as.factor(object$optimal_model_results$data[,ncol(object$optimal_model_results$data)])
  classification<-table(clustering)
  obj <- list(title = title,
              C=object$C,
              N=object$N,
              R=object$R,
              K = clust_count,
              modelName = object$optimal_model,
              loglik = loglik,
              information_criterion=object$information_criterion,
              ic_output=ic_value,
              classification=classification)
  class(obj) <- "summary.betaclust"
  return(obj)
}

print.summary.betaclust <- function(x, digits = getOption("digits"), ...)
{
  txt <- paste(rep("-", min(nchar(x$title), getOption("width"))), collapse = "")
  cat(txt)
  cat("\n")
  cat(x$title)
  cat("\n")
  cat(txt)
  #
  cat("\n")
  if(x$K == 0)
  {
    cat("betaclust model with only a noise component:")
  } else
  {
    cat(paste0("betaclust ", x$modelName, " model with ",
               x$K, ifelse(x$K > 1, " components", " component"),
               ":"))
  }
  cat("\n")



  tab <- data.frame("log-likelihood" = x$loglik,  "Information-criterion"=x$information_criterion,"IC-value" = x$ic_output,
                    "CpG-sites" = x$C, "Patients" = x$N, "Samples" = x$R,
                    row.names = "", check.names = FALSE)
  print(tab)
  #
  cat("\nClustering table:")
  #print(table(x$classification))
  print(x$classification)
  #
  invisible(x)
}
