#' @title Bayesian Information Criterion
#' @description Compute the BIC value for the optimal model.
#' @details Computes the BIC for a specified model given the log-likelihood, the dimension of the data, and the model specification.
#' @export
#' @seealso \code{\link{em_aic}}
#' @seealso  \code{\link{em_icl}}
#' @param llk log-likelihood value.
#' @param C number of CpG sites.
#' @param K number of methylation profiles to be identified.
#' @param patients number of patients in the study.
#' @param samples number of samples collected from each patient for study.
#' @param model_name mixture model (method=c("C..","CN.","C.R")).
#' @return The BIC value for the selected model.

em_bic<-function(llk,C,K,patients=4,samples=1,model_name="C.."){

  R=samples
  N=patients
  mod_len=length(model_name)
  bic=vector("numeric",mod_len)
  #print(model_names)
  num_par=0

  for(i in 1:mod_len)
  {
    num_par=0
    if(model_name[i] ==  "C..") ##C
    {
      num_par = (K*2)+(K-1)
      #print(num_par)
    }else if(model_name[i] == "CN.") ##CN
    {
      num_par = (K*N*2)+(K-1)
      #print(num_par)
    }else if(model_name[i] == "C.R") ##CR
    {
      num_par = (K*R*2)+(K-1)
      #print(num_par)
    }

    if(is.na(llk[i]))
    {
      bic[i]=NA
    }else{
      bic[i]= -(2*llk[i])+(num_par * log(C))
    }
  }

  return(bic)
}
