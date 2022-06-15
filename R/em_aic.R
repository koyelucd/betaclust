#' @title Akaike Information Criterion
#' @description Compute the AIC value for the optimal model.
#' @details Computes the AIC for a specified model given the log-likelihood, the dimension of the data, and the model specification.
#' @export
#' @seealso \code{\link{em_bic}}
#' @seealso \code{\link{em_icl}}
#' @param llk log-likelihood value.
#' @param C number of CpG sites.
#' @param K number of methylation profiles to be identified.
#' @param patients number of patients in the study.
#' @param samples number of samples collected from each patient for study.
#' @param model_name mixture model (method=c("C..","CN.","C.R")).
#' @return The AIC value for the selected model.

em_aic<-function(llk,C,K,patients=4,samples=1,model_name="C.."){

  R=samples
  N=patients
  mod_len=length(model_name)
  aic=vector("numeric",mod_len)

  for(i in 1:mod_len)
  {
    if(model_name[i] ==  "C..") ##C
    {
      num_par = (K*2)+(K-1)
    }else if(model_name[i] == "CN.") ##CN
    {
      num_par = (K*N*R*2)+(K-1)
    }else if(model_name[i] == "C.R") ##CR
    {
      num_par = (K*R*2)+(K-1)
    }

    if(is.na(llk[i]))
    {
      aic[i]=NA
    }else{
      aic[i]= -(2*llk[i])+(2*num_par)
    }
  }

  return(aic)
}
