#' @title Akaike Information Criterion
#' @description Compute the AIC value for the optimal model.
#' @details Computes the AIC for a specified model given the log-likelihood, the dimension of the data, and the model specification.
#' @export
#' @seealso \code{\link{em_bic}}
#' @seealso \code{\link{em_icl}}
#' @param llk log-likelihood value.
#' @param C number of CpG sites.
#' @param M number of methylation profiles identified in a DNA sample.
#' @param patients number of patients.
#' @param samples number of DNA samples collected from each patient.
#' @param model_name fitted mixture model (method=c("K..","KN.","K.R")).
#' @return The AIC value for the selected model.

em_aic<-function(llk,C,M,patients=4,samples=1,model_name="K.."){

  R=samples
  N=patients
  mod_len=length(model_name)
  aic=vector("numeric",mod_len)

  K=M
  for(i in 1:mod_len)
  {
    if(model_name[i] ==  "K..") ##K
    {
      num_par = (K*2)+(K-1)
    }else if(model_name[i] == "KN.") ##KN
    {
      num_par = (K*N*R*2)+(K-1)
    }else if(model_name[i] == "K.R") ##KR
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
