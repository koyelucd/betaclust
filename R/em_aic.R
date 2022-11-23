#' @title Akaike Information Criterion
#' @description Compute the AIC value for the optimal model.
#' @details Computes the AIC for a specified model given the log-likelihood, the dimension of the data, and the model names.
#' @export
#' @seealso \code{\link{em_bic}}
#' @seealso \code{\link{em_icl}}
#' @param llk Log-likelihood value.
#' @param C Number of CpG sites.
#' @param M Number of methylation states identified in a DNA sample.
#' @param N Number of patients.
#' @param R Number of DNA samples collected from each patient.
#' @param model_name Fitted mixture model. Options are "K..", "KN." and/or "K.R" (default = "K..").
#' @return The AIC value for the selected model.

em_aic<-function(llk,C,M,N,R,model_name="K.."){

  #R=samples
  #N=patients
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
