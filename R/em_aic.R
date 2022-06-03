#' @title Akaike Information Criterion
#' @description The AIC value used to select the optimal model
#' @export
#' @param llk log-likelihood value
#' @param C number of CpG sites
#' @param K number of clusters
#' @param patients number of patients
#' @param samples no. of samples
#' @param model_names mixture model (method=c("C..","CN.","C.R"))
#' @return The AIC value for the selected model

em_aic<-function(llk,C,K,patients=4,samples=1,model_names="C.."){

  R=samples
  N=patients
  mod_len=length(model_names)
  aic=vector("numeric",mod_len)

  for(i in 1:mod_len)
  {
    if(model_names[i] ==  "C..") ##C
    {
      num_par = (K*2)+(K-1)
    }else if(model_names[i] == "CN.") ##CN
    {
      num_par = (K*N*R*2)+(K-1)
    }else if(model_names[i] == "C.R") ##CR
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
