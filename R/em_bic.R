#' @title Bayesian Information Criterion
#' @description Compute the BIC value for the optimal model.
#' @details Computes the BIC for a specified model given the log-likelihood, the dimension of the data, and the model names.
#' @export
#' @seealso \code{\link{em_aic}}
#' @seealso  \code{\link{em_icl}}
#' @param llk Log-likelihood value.
#' @param C Number of CpG sites.
#' @param M Number of methylation states identified in a DNA sample.
#' @param N Number of patients.
#' @param R Number of DNA samples collected from each patient.
#' @param model_name Fitted mixture model. Options are "K..", "KN." and/or "K.R" (default = "K..").
#' @return The BIC value for the selected model.

em_bic<-function(llk,C,M,N,R,model_name="K.."){

  #R=samples
  #N=patients
  K=M
  mod_len=length(model_name)
  bic=vector("numeric",mod_len)
  #print(model_names)
  num_par=0

  for(i in 1:mod_len)
  {
    num_par=0
    if(model_name[i] ==  "K..") ##K
    {
      num_par = (K*2)+(K-1)
      #print(num_par)
    }else if(model_name[i] == "KN.") ##KN
    {
      num_par = (K*N*2)+(K-1)
      #print(num_par)
    }else if(model_name[i] == "K.R") ##KR
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
