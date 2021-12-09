#' BIC to compare the methods
#' @export
#' @param llk log-likelihood value
#' @param C number of CpG sites
#' @param K number of clusters
#' @param patients number of patients
#' @param samples no. of samples
#' @param model_names mixture model (method=c("C..","CN.","C.R"))
#' @return bic

em_bic<-function(llk,C,K,patients=4,samples=1,model_names="C.."){

  R=samples
  N=patients
  mod_len=length(model_names)
  bic=vector("numeric",mod_len)
  #print(model_names)
  num_par=0

  for(i in 1:mod_len)
  {
    num_par=0
    if(model_names[i] ==  "C..") ##C
    {
      num_par = (K*2)+(K-1)
      #print(num_par)
    }else if(model_names[i] == "CN.") ##CN
    {
      num_par = (K*N*2)+(K-1)
      #print(num_par)
    }else if(model_names[i] == "C.R") ##CR
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
