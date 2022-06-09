#' @title The betaclust wrapper function
#' @export
#' @description A family of model based clustering techniques to identify the methylation profiles of the beta valued DNA methylation data
#'
#' @details This is a wrapper function which can be used to run all three models (C.., CN., C.R) together.
#' The C.. and CN. models are used to analyse a single DNA sample and cluster the CpG sites into the 3 methylation profiles (hypomethylation, hemimethylation, hypermethylation).
#' The thresholds can be objectively identified from the clustering solution. The C. R model is used to analyse R samples to
#' the differentially methylated CpG sites between R DNA samples.
#'
#' @seealso \code{\link{beta_c}}
#' @seealso \code{\link{beta_cn}}
#' @seealso \code{\link{beta_cr}}
#' @seealso \code{\link{pca.methylation.data}}
#' @seealso \code{\link{plot.betaclust}}
#' @seealso \code{\link{summary.betaclust}}
#'
#' @param X methylation values for CpG sites frpm R samples collected from N patients
#' @param K number of methylation groups to be identified (default=3)
#' @param patients number of patients in the study
#' @param samples number of samples collected from each patient for study
#' @param model_names mixture model to run (Models= c(C..,CN.,C.R), default=C..)
#' @param model_selection optimal model selection based on information criterion. (Methods=AIC,BIC,ICL,default=BIC)
#' @param seed seed for reproducible work
#' @param register setting for parallelization
#'
#' @return The function returns an object of "betaclust" class.
#' The class object contains following values:
#' \itemize{
#' \item Information_criterion - The information criterion used to select the optimal model.
#' \item ic_output - This stores the information criterion value calculated for each model.
#' \item optimal_model - The model selected as optimal.
#' \item function_call - The parameters passed as arguments to the function betaclust.
#' \item CpG_sites - The number of CpG sites analysed using the beta mixture models.
#' \item patients - The number of patients analysed using the beta mixture models.
#' \item samples - The numder of samples analysed using the beta mixture models.
#' \item best_model - This contains the final results for the optimal model selected. Thus this contains the following values:
#'    \itemize{
#'    \item cluster_count - The total number of CpG sites identified in each cluster.
#'    \item llk - The vector containing log-likelihood values calculated for each step of parameter estimation.
#'    \item data - This contains the methylation dataset along with the cluster label as determined by the mixture model.
#'    \item alpha - This contains the shape parameter 1 for the beta mixtures for \eqn{K^R} groups.
#'    \item beta - This contains the shape parameter 2 for the beta mixtures for \eqn{K^R} groups.
#'    \item tau - The proportion of CpG sites in each cluster.
#'    \item z - The matrix contains the probability calculated for each CpG site belonging to the \eqn{K^R} clusters.
#'    \item uncertainty - The uncertainty of a CpG site belonging to the identified cluster.
#'    }
#' }
#'
#' @examples
#' \dontrun{
#' data(pca.methylation.data)
#' my.seed=190
#' K=3
#' patients=4
#' samples=2
#' data_output=betaclust(pca.methylation.data[,2:9],K,patients,samples,
#'             model_names=c("C..","CN.","C.R"),model_selection="BIC",seed=my.seed)
#' }
#'
#'
#' @importFrom foreach %dopar%
#' @importFrom stats C
#' @importFrom utils txtProgressBar
#' @references {Silva, R., Moran, B., Russell, N.M., Fahey, C., Vlajnic, T., Manecksha, R.P., Finn, S.P., Brennan, D.J., Gallagher, W.M., Perry, A.S.: Evaluating liquid biopsies for methylomic profiling of prostate cancer. Epigenetics 15(6-7), 715-727 (2020). doi:10.1080/15592294.2020.1712876.}
#' @references {Fraley, C., Raftery, A.E.: How many clusters? which clustering method? answers via model-based cluster analysis. The computer journal 41, 578-588 (1998). doi: 10.1093/comjnl/41.8.578.}
#' @references {Dempster, A., Laird, N., Rubin, D.: Maximum likelihood from incomplete data via the em algorithm. Journal of the royal statistical society 39(1), 1-38 (1977). doi:10.1111/j.2517-6161.1977.tb01600.x.}
#' @references {Diamond, H.G., Straub, A.: Bounds for the logarithm of the euler gamma function and its derivatives. Journal of mathematical analysis and applications 433(2), 1072-1083 (2016).doi:10.1016/j.jmaa.2015.08.034.}


betaclust<-function(data,K=3,patients,samples,model_names="C..",model_selection="BIC",seed,register=NULL){

  X=data
  len=length(model_names)
  llk<-vector()
  C=nrow(X)
  z<-vector(mode = "list", length = len)

  K_len=K
  model_len=length(model_names)

  ## Progress bar
  cat("fitting ...\n")
  #flush.console()
  pbar <- utils::txtProgressBar(min = 0, max = model_len, style = 3)
  on.exit(close(pbar))


  #print(X[1:3,])
  if(length(model_names))
  {
    for(i in 1:len)
    {
      ## Call for C.. function
      if(model_names[i] == "C..")
      {
        #print(model_names[i])
        ## check if samples>1
        if(samples>1)
        {
          warning("C.. method can run for a single sample and not multiple samples.
                  Method is run for 1st sample only.", call. = FALSE)
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
      if(model_names[i] == "CN.")
      {
        ## check if samples>1
        if(samples>1)
          warning("CN. method can run for a single sample and not multiple samples.
                  Method is run for 1st sample only.", call. = FALSE)
        call_data<-X[,1:patients]
        cn_out<-beta_cn(call_data,K,seed,register)
        len_llk<-length(cn_out$llk)
        llk<-c(llk,cn_out$llk[len_llk])
        z[[i]]<-cn_out$z
        #print(cn_out$alpha)
        #print(cn_out$beta)
      }
      ## Call for C.R function
      if(model_names[i] == "C.R")
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
      ##Set Progress bar
      utils::setTxtProgressBar(pbar,i)
    }

  }

  # print(llk)
  # print(model_names)
  # print(model_selection)
  ### Optimal model selection and output

  if(model_selection=="BIC")
  {
    ## compare bic value
    ic_op<-em_bic(llk,C,K,patients,samples,model_names)
    min_index<-which.min(ic_op)
    min_method<-model_names[min_index]

  }else if(model_selection=="AIC")
  {
    ## compare aic value
    ic_op<-em_aic(llk,C,K,patients,samples,model_names)
    min_index<-which.min(ic_op)
    min_method<-model_names[min_index]
  }else if(model_selection=="ICL")
  {
    ## compare icl value
    ic_op<-em_icl(llk,C,K,patients,samples,model_names,z)
    min_index<-which.min(ic_op)
    min_method<-model_names[min_index]
  }

  call_function<-match.call()

  N=patients
  R=0
  if(min_method=="C..")
  {
    final_output<-c_out
    R=1
  }else if(min_method=="CN.")
  {
    final_output<-cn_out
    R=1
  }else if(min_method=="C.R")
  {
    final_output<-cr_out
    R=samples
  }



  #print("Execution is complete")

  beta_out<-list(information_criterion=model_selection,ic_output=ic_op,
                 optimal_model=min_method,function_call=call_function,CpG_sites=C,patients=N,samples=R,
                 best_model=final_output)
  class(beta_out)<-"betaclust"
  return(beta_out)
}




