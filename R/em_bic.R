#' BIC to compare the methods
#' @export
#' @param llk log-likelihood value
#' @param C number of CpG sites
#' @param K number of clusters
#' @param patients number of patients
#' @param samples no. of samples
#' @param mixm mixture model (method=c("C..","CN.","C.R"))
#' @return bic

em_bic<-function(llk,C,K,patients=4,samples=1,mixm="C.."){

  R=samples
  N=patients
  mod_len=length(mixm)
  bic=vector("numeric",mod_len)

  for(i in 1:mod_len)
  {
    num_par=0
    if(mixm[i] ==  "C..") ##C
    {
      num_par = (K*2)+(K-1)
      print(num_par)
    }else if(mixm[i] == "CN.") ##CN
    {
      num_par = (K*N*2)+(K-1)
      print(num_par)
    }else if(mixm[i] == "C.R") ##CR
    {
      num_par = (K*R*2)+(K-1)
      print(num_par)
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
