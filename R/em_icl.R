globalVariables(c("patients","samples"))
#' @title Integrated Complete-data Likelihood (ICL) Criterion
#' @description Compute the ICL value for the optimal model.
#'
#' @details Computes the ICL for a specified model given the log-likelihood, the dimension of the data, and the model specification. This criterion penalises the BIC by including an entropy term favouring the
#' well separated clusters.
#' @export
#' @seealso \code{\link{em_aic}}
#' @seealso  \code{\link{em_bic}}
#' @param llk log-likelihood value.
#' @param C number of CpG sites.
#' @param M number of methylation states identified in a DNA sample.
#' @param N number of patients.
#' @param R number of DNA samples collected from each patient.
#' @param model_name fitted mixture model (model_name = c("K..","KN.","K.R")).
#' @param z a matrix of posterior probability of cluster membership, stored as z in the object from \code{\link[betaclust]{beta_k}}/\code{\link[betaclust]{beta_kn}}/\code{\link[betaclust]{beta_kr}} functions.
#' @return The ICL value for the selected model.

em_icl<-function(llk,C,M,N,R,model_name="K..",z){

  #R=samples
  #N=patients
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
