#' @title The betaclust wrapper function
#' @export
#' @description A family of model based clustering techniques to identify methylation profiles in beta valued DNA methylation data.
#'
#' @details This is a wrapper function which can be used to fit all three models (C.., CN., C.R) together.
#' The C.. and CN. models are used to analyse a single DNA sample (\eqn{R = 1}) and cluster the \eqn{C} CpG sites into the \eqn{K} methylation profiles. As each CpG site can belong to either of the \eqn{M=3} methylation profiles (hypomethylation, hemimethylation and hypermethylation), the default value for \eqn{K=M=3}.
#' The thresholds between methylation profiles can be objectively identified from the clustering solution.
#' The C.R model is used to analyse \eqn{R} independent samples collected from \eqn{N} patients, where each sample contains \eqn{C} CpG sites, and cluster
#' the dataset to identify the differentially methylated CpG sites between the \eqn{R} DNA samples.
#'
#' @seealso \code{\link{beta_c}}
#' @seealso \code{\link{beta_cn}}
#' @seealso \code{\link{beta_cr}}
#' @seealso \code{\link{pca.methylation.data}}
#' @seealso \code{\link{plot.betaclust}}
#' @seealso \code{\link{summary.betaclust}}
#'
#' @param data Methylation values for \eqn{C} CpG sites from \eqn{R} samples collected from \eqn{N} patients.
#' @param K Number of methylation profiles to be identified.
#' @param patients Number of patients in the study.
#' @param samples Number of samples collected from each patient for study.
#' @param model_names Models to run from the set of models, C.., CN. and C.R, default = C.. . See details.
#' @param model_selection Information criterion used for model selection. Options are AIC/BIC/ICL/default=BIC.
#' @param seed Seed to allow for reproducibility.
#' @param register Setting for registering the parallel backend with the "foreach" package. To start parallel execution of R code on machine with multiple cores, "NULL" value needs to be assigned to this parameter.
#'
#' @return The function returns an object of "betaclust" class which contains the following values:
#' \itemize{
#' \item information_criterion - the information criterion used to select the optimal model.
#' \item ic_output - this stores the information criterion value calculated for each model.
#' \item optimal_model - the model selected as optimal.
#' \item function_call - the parameters passed as arguments to the function betaclust.
#' \item C - the number of CpG sites analysed using the beta mixture models.
#' \item N - the number of patients analysed using the beta mixture models.
#' \item R - the numder of samples analysed using the beta mixture models.
#' \item optimal_model_results - this contains information from the optimal model. Specifically,
#'    \itemize{
#'    \item cluster_size - the total number of CpG sites identified in each cluster.
#'    \item llk - a vector containing the log-likelihood value at each step of the EM algorithm.
#'    \item data - this contains the methylation dataset along with the cluster label for each CpG site.
#'    \item alpha - this contains the shape parameter 1 for the beta mixture model.
#'    \item delta - this contains the shape parameter 2 for the beta mixture model.
#'    \item tau - the proportion of CpG sites in each cluster.
#'    \item z - a matrix containing the probability for each CpG site of belonging to each of the \eqn{K} clusters.
#'    \item uncertainty - the uncertainty of each CpG site's clustering.
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
#' @references {Majumdar, K., Silva, R., Perry, A.S., Watson, R.W., Murphy, T.B., Gormley, I.C.: betaclust: a family of mixture models for beta valued DNA methylation data.}


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
          warning("C.. model only considers a single sample and not multiple samples.
                  Model is fitted to 1st sample only.", call. = FALSE)
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
          warning("CN. model only considers a single sample and not multiple samples.
                  Model is fitted to 1st sample only.", call. = FALSE)
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
          warning("C.R considers only multiple samples. Please pass DNA methylation values for more than 1 sample.", call. = FALSE)
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
                 optimal_model=min_method,function_call=call_function,C=C,N=N,R=R,
                 optimal_model_results=final_output)
  class(beta_out)<-"betaclust"
  return(beta_out)
}




