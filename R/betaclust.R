#' @title The betaclust wrapper function
#' @export
#' @description A family of model-based clustering techniques to identify methylation states in beta-valued DNA methylation data.
#'
#' @details This is a wrapper function which can be used to fit all three models (K.., KN., K.R) within a single function.
#'
#' The K.. and KN. models are used to analyse a single DNA sample type (\eqn{R = 1}) and cluster the \eqn{C} CpG sites into the \eqn{K} clusters which represent the different methylation states in a DNA sample type. As each CpG site can belong to any of the \eqn{M=3} methylation states (hypomethylation, hemimethylation and hypermethylation), the default value for \eqn{K=M=3}.
#' The thresholds between methylation states are objectively inferred from the clustering solution.
#'
#' The K.R model is used to analyse \eqn{R} independent sample types collected from \eqn{N} patients, where each sample contains \eqn{C} CpG sites, and cluster
#' the dataset into \eqn{K=M^R} clusters to identify the differentially methylated CpG (DMC) sites between the \eqn{R} DNA sample types.
#'
#' @seealso \code{\link{beta_k}}
#' @seealso \code{\link{beta_kn}}
#' @seealso \code{\link{beta_kr}}
#' @seealso \code{\link{pca.methylation.data}}
#' @seealso \code{\link{plot.betaclust}}
#' @seealso \code{\link{summary.betaclust}}
#' @seealso \code{\link{threshold}}
#'
#' @param data A dataframe of dimension \eqn{C \times NR} containing methylation values for \eqn{C} CpG sites from \eqn{R} sample types collected from \eqn{N} patients.
#' Samples are grouped together in the dataframe such that the columns are ordered as Sample1_Patient1, Sample1_Patient2, Sample2_Patient1, Sample2_Patient2, etc.
#' @param M Number of methylation states to be identified in a DNA sample type.
#' @param N Number of patients in the study.
#' @param R Number of sample types collected from each patient for the study.
#' @param model_names Models to run from the set of models, K.., KN. and K.R, default = K.. . See details.
#' @param model_selection Information criterion used for model selection. Options are AIC, BIC or ICL (default = BIC).
#' @param parallel_process The "TRUE" option results in parallel processing of the models for increased computational efficiency. The default option has been set as "FALSE" due to package testing limitations.
#' @param seed Seed to allow for reproducibility (default = NULL).
#'
#' @return The function returns an object of the \code{\link[betaclust:betaclust]{betaclust}} class which contains the following values:
#' \itemize{
#' \item information_criterion - The information criterion used to select the optimal model.
#' \item ic_output - The information criterion value calculated for each model.
#' \item optimal_model - The model selected as optimal.
#' \item function_call - The parameters passed as arguments to the function \code{\link[betaclust:betaclust]{betaclust}}.
#' \item K - The number of clusters identified using the beta mixture models.
#' \item C - The number of CpG sites analysed using the beta mixture models.
#' \item N - The number of patients analysed using the beta mixture models.
#' \item R - The number of sample types analysed using the beta mixture models.
#' \item optimal_model_results - Information from the optimal model. Specifically,
#'    \itemize{
#'    \item cluster_size - The total number of CpG sites in each of the K clusters.
#'    \item llk - A vector containing the log-likelihood value at each step of the EM algorithm.
#'    \item alpha - This contains the first shape parameter for the beta mixture model.
#'    \item delta - This contains the second shape parameter for the beta mixture model.
#'    \item tau - The proportion of CpG sites in each cluster.
#'    \item z - A matrix of dimension \eqn{C \times K} containing the posterior probability of each CpG site belonging to each of the \eqn{K} clusters.
#'    \item classification - The classification corresponding to z, i.e. map(z).
#'    \item uncertainty - The uncertainty of each CpG site's clustering.
#'    \item thresholds - Threshold points calculated under the K.. or the KN. model.
#'    \item DM - The AUC and WD metric for distribution similarity in each cluster.
#'    }
#' }
#'
#' @examples
#' \donttest{
#' my.seed <- 190
#' M <- 3
#' N <- 4
#' R <- 2
#' data_output <- betaclust(pca.methylation.data[1:30,2:9], M, N, R,
#'             model_names = c("K..","KN.","K.R"), model_selection = "BIC",
#'             parallel_process = FALSE, seed = my.seed)
#'
#'}
#' @importFrom foreach %dopar%
#' @importFrom stats C
#' @importFrom utils flush.console
#' @importFrom utils txtProgressBar
#' @references {Silva, R., Moran, B., Russell, N.M., Fahey, C., Vlajnic, T., Manecksha, R.P., Finn, S.P., Brennan, D.J., Gallagher, W.M., Perry, A.S.: Evaluating liquid biopsies for methylomic profiling of prostate cancer. Epigenetics 15(6-7), 715-727 (2020). \doi{10.1080/15592294.2020.1712876}.}
#' @references {Majumdar, K., Silva, R., Perry, A.S., Watson, R.W., Murphy, T.B., Gormley, I.C.: betaclust: a family of mixture models for beta valued DNA methylation data. arXiv [stat.ME] (2022). \doi{10.48550/ARXIV.2211.01938}.}


betaclust<-function(data,M=3,N,R,model_names="K..",model_selection="BIC",parallel_process=FALSE,seed=NULL){

  X=as.data.frame(data)
  len=length(model_names)
  llk<-vector()
  C=nrow(X)
  z<-vector(mode = "list", length = len)

  M_len=M
  model_len=length(model_names)

  ## Progress bar

  message("fitting ...\n")
  flush.console()
    pbar <- utils::txtProgressBar(min = 0, max = model_len, style = 3)
    on.exit(close(pbar))

  if(length(model_names))
  {
    for(i in 1:len)
    {
      ## Call for K.. function
      if(model_names[i] == "K..")
      {
        ## check if R>1
        if(R>1)
        {
          warning("K.. model only considers a single sample type and not multiple sample types.
                  Model is fitted to 1st sample type only.", call. = FALSE)
        }

        call_data=X[,1:N]
        k_out<-beta_k(call_data,M,parallel_process,seed)
        len_llk<-length(k_out$llk)
        llk<-c(llk,k_out$llk[len_llk])
        z[[i]]<-k_out$z
      }
      ## Call for KN. function
      if(model_names[i] == "KN.")
      {
        ## check if R>1
        if(R>1)
          warning("KN. model only considers a single sample type and not multiple sample types.
                  Model is fitted to 1st sample type only.", call. = FALSE)
        call_data<-X[,1:N]
        kn_out<-beta_kn(call_data,M,parallel_process,seed)
        len_llk<-length(kn_out$llk)
        llk<-c(llk,kn_out$llk[len_llk])
        z[[i]]<-kn_out$z
      }
      ## Call for K.R function
      if(model_names[i] == "K.R")
      {
        ## check if R>1
        if(R<=1){
          warning("K.R considers only multiple sample types. Please pass DNA methylation values for more than 1 sample type.", call. = FALSE)
          llk<-c(llk,NA)
        }else if(ncol(X)<(N*R)){
          warning("K.R considers only multiple sample types. Please pass DNA methylation values for more than 1 sample type.", call. = FALSE)
          llk<-c(llk,NA)
          }else{
          kr_out<-beta_kr(X,M,N,R,parallel_process,seed)
          len_llk<-length(kr_out$llk)
          llk<-c(llk,kr_out$llk[len_llk])
          z[[i]]<-kr_out$z
        }
      }
      ##Set Progress bar
      utils::setTxtProgressBar(pbar,i)
    }

  }


  ### Optimal model selection and output

  if(model_selection=="BIC")
  {
    ## compare bic value
    ic_op<-em_bic(llk,C,M,N,R,model_names)
    min_index<-which.min(ic_op)
    min_method<-model_names[min_index]

  }else if(model_selection=="AIC")
  {
    ## compare aic value
    ic_op<-em_aic(llk,C,M,N,R,model_names)
    min_index<-which.min(ic_op)
    min_method<-model_names[min_index]
  }else if(model_selection=="ICL")
  {
    ## compare icl value
    ic_op<-em_icl(llk,C,M,N,R,model_names,z)
    min_index<-which.min(ic_op)
    min_method<-model_names[min_index]
  }

  call_function<-match.call()


  R_old=R
  K=0
  if(min_method=="K..")
  {
    final_output<-k_out
    thresholds<-threshold(k_out,call_data,min_method)
    final_output$thresholds<-thresholds
    K=M
    R=1
  }else if(min_method=="KN.")
  {
    final_output<-kn_out
    thresholds<-threshold(kn_out,call_data,min_method)
    final_output$thresholds<-thresholds
    K=M
    R=1
  }else if(min_method=="K.R")
  {
    final_output<-kr_out
    K=M^R_old
    R=R_old
  }


  beta_out<-list(information_criterion=model_selection,ic_output=ic_op,
                 optimal_model=min_method,function_call=call_function,K=K,C=C,N=N,R=R,
                 optimal_model_results=final_output)
  class(beta_out)<-"betaclust"
  return(beta_out)
}




