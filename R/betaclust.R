#' betaclust wrapper function
#' @export
#' @param  X methylation values for CpG sites frpm R samples collected from N patients
#' @param K number of methylation groups to be identified (default=3)
#' @param patients number of patients in the study
#' @param samples number of samples collected from each patient for study
#' @param mixm mixture model to run (Models= c(C..,CN.,C.R), default=C..)
#' @param model_selection optimal model selection based on information criterion. (Methods=AIC,BIC,ICL,default=BIC)
#' @param seed seed for reproducible work
#' @param register setting for parallelization
#' @importFrom foreach %dopar%
#' @importFrom stats C

betaclust<-function(X,K=3,patients,samples,mixm="C..",model_selection="BIC",seed,register=NULL){
  len=length(mixm)
  llk<-vector()
  C=nrow(X)
  z<-vector(mode = "list", length = len)
  #print(X[1:3,])
  if(length(mixm))
  {
    for(i in 1:len)
    {
      ## Call for C.. function
      if(mixm[i] == "C..")
      {
        #print(mixm[i])
        ## check if samples>1
        if(samples>1)
        {
          warning("C.. method can run for single sample and not multiple samples.
                  Method run for 1st sample only.", call. = FALSE)
        }

        call_data=X[,1:patients]
        #print(call_data[1:4,])
        c_out<-beta_c(call_data,K,seed,register)
        len_llk<-length(c_out$llk)
        llk<-c(llk,c_out$llk[len_llk])
        z[[i]]<-c_out$z
        #print(c_out$alpha)
        #print(c_out$beta)
      }
      ## Call for CN. function
      if(mixm[i] == "CN.")
      {
        ## check if samples>1
        if(samples>1)
          warning("CN. method can run for single sample and not multiple samples.
                  Method run for 1st sample only.", call. = FALSE)
        call_data<-X[,1:patients]
        cn_out<-beta_cn(call_data,K,seed,register)
        len_llk<-length(cn_out$llk)
        llk<-c(llk,cn_out$llk[len_llk])
        z[[i]]<-cn_out$z
        #print(cn_out$alpha)
        #print(cn_out$beta)
      }
      ## Call for C.R function
      if(mixm[i] == "C.R")
      {
        ## check if samples>1
        if(samples<=1){
          warning("C.R can only run for multiple samples. Please pass DNA methylation values for more than 1 sample.", call. = FALSE)
          llk<-c(llk,NA)
          }else {
          cr_out<-beta_cr(X,K,patients,samples,seed,register)
          len_llk<-length(cr_out$llk)
          llk<-c(llk,cr_out$llk[len_llk])
          z[[i]]<-cr_out$z
          #print(cr_out$alpha)
          #print(cr_out$beta)
        }
      }
    }

  }

  # print(llk)
  # print(mixm)
  # print(model_selection)
  ### Optimal model selection and output

  if(model_selection=="BIC")
  {
    ## compare bic value
    ic_op<-em_bic(llk,C,K,patients,samples,mixm)
    min_index<-which.min(ic_op)
    min_method<-mixm[min_index]

  }else if(model_selection=="AIC")
  {
    ## compare aic value
    ic_op<-em_aic(llk,C,K,patients,samplesmixm)
    min_index<-which.min(ic_op)
    min_method<-mixm[min_index]
  }else if(model_selection=="ICL")
  {
    ## compare icl value
    ic_op<-em_icl(llk,C,K,patients,samples,mixm,z)
    min_index<-which.min(ic_op)
    min_method<-mixm[min_index]
  }

  call_function<-match.call()

  if(min_method=="C..")
  {
    final_output<-c_out
  }else if(min_method=="CN.")
  {
    final_output<-cn_out
  }else if(min_method=="C.R")
  {
    final_output<-cr_out
  }

  print("Execution is complete")

  return(list(information_criterion=model_selection,ic_output=ic_op,optimal_model=min_method,function_call=call_function,best_model=final_output))
}


