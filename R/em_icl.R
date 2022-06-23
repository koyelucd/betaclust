#' @title Integrated Complete-data Likelihood (ICL) Criterion
#' @description Compute the ICL value for the optimal model.
#'
#' @details Computes the ICL for a specified model given the log-likelihood, the dimension of the data, and the model specification. This criterion penalises the BIC by including the entropy term favouring the
#' well separated clusters.
#' @export
#' @seealso \code{\link{em_aic}}
#' @seealso  \code{\link{em_bic}}
#' @param llk log-likelihood value.
#' @param C number of CpG sites.
#' @param M number of methylation profiles identified in a DNA sample.
#' @param patients number of patients.
#' @param samples number of DNA samples collected from each patient.
#' @param model_name fitted mixture model (method=c("C..","CN.","C.R")).
#' @param z z matrix used for computing the complete-data log-likelihood function.
#' @return The ICL value for the selected model.

em_icl<-function(llk,C,M,patients=4,samples=1,model_name="C..",z){

  R=samples
  N=patients
  K=M
  mod_len=length(model_name)
  icl=vector("numeric",mod_len)

  bic<-em_bic(llk,C,K,model_name,patients,samples)
  for(i in 1:mod_len)
  {


    D <- matrix(0, C, ncol(z[[i]]))
    z_i<-z[[i]]
    for(j in 1:C)
      D[j, which.max(z_i[j,])] <- 1


    if(is.na(llk[i]))
    {
      icl[i]=NA
    }else{
      icl[i]=bic[i] + 2*sum(D * ifelse(z_i > 0, log(z_i), 0))
    }
  }

  return(icl)
}
