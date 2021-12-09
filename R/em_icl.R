#' ICL to compare the methods
#' @export
#' @param llk log-likelihood value
#' @param C number of CpG sites
#' @param K number of clusters
#' @param patients number of patients
#' @param samples no. of samples
#' @param model_names mixture model (method=c("C..","CN.","C.R"))
#' @param z z matrix for each output
#' @return icl

em_icl<-function(llk,C,K,patients=4,samples=1,model_names="C..",z){

  R=samples
  N=patients
  mod_len=length(model_names)
  icl=vector("numeric",mod_len)

  bic<-em_bic(llk,C,K,model_names,patients,samples)
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
