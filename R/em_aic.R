#' AIC to compare the methods
#' @export
#' @param llk log-likelihood value
#' @param C number of CpG sites
#' @param K number of clusters
#' @param patients number of patients
#' @param samples no. of samples
#' @param mixm mixture model (method=c("C..","CN.","C.R"))
#' @return aic

em_aic<-function(llk,C,K,patients=4,samples=1,mixm="C.."){

  R=samples
  N=patients
  mod_len=length(mixm)
  aic=vector("numeric",mod_len)

  for(i in 1:mod_len)
  {
    if(mixm[i] ==  "C..") ##C
    {
      num_par = (K*2)+(K-1)
    }else if(mixm[i] == "CN.") ##CN
    {
      num_par = (K*N*R*2)+(K-1)
    }else if(mixm[i] == "C.R") ##CR
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
