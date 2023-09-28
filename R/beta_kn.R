globalVariables(c("k"))
#' @title Fit the KN. model
#' @description Fit the KN. model from the \code{\link[betaclust:betaclust]{betaclust}} family of beta mixture models for DNA methylation data.
#'              The KN. model analyses a single DNA sample type and identifies the thresholds between the different methylation states.
#'
#' @export
#'
#' @details The KN. model clusters each of the \eqn{C} CpG sites into one of \eqn{K} methylation states, based on data from \eqn{N} patients for one DNA sample type (i.e. \eqn{R = 1}).
#' As each CpG site can belong to any of the \eqn{M = 3} methylation states (hypomethylated, hemimethylated or hypermethylated), the default value of \eqn{K = M = 3}.
#' The KN. model differs from the K.. model as it is less parsimonious, allowing cluster and patient-specific shape parameters. The returned object can be passed as an input parameter to the
#' \code{\link[betaclust:threshold]{threshold}} function available in this package to calculate the thresholds between the methylation states.
#'
#' @seealso \code{\link{beta_k}}
#' @seealso \code{\link{betaclust}}
#' @seealso \code{\link{threshold}}
#'
#' @param data A dataframe of dimension \eqn{C \times N} containing methylation values for \eqn{C} CpG sites from \eqn{R = 1} sample type collected from \eqn{N} patients.
#' Samples are grouped together in the dataframe such that the columns are ordered as Sample1_Patient1, Sample1_Patient2, etc.
#' @param M Number of methylation states to be identified in a DNA sample type.
#' @param parallel_process The "TRUE" option results in parallel processing of the models for increased computational efficiency. The default option has been set as "FALSE" due to package testing limitations.
#' @param seed Seed to allow for reproducibility (default = NULL).
#'
#' @return A list containing:
#' \itemize{
#'    \item cluster_size - The total number of CpG sites in each of the K clusters.
#'    \item llk - A vector containing the log-likelihood value at each step of the EM algorithm.
#'    \item alpha - The first shape parameter for the beta mixture model.
#'    \item delta - The second shape parameter for the mixture model.
#'    \item tau - The estimated mixing proportion for each cluster.
#'    \item z - A matrix of dimension \eqn{C \times K} containing the posterior probability of each CpG site belonging to each of the \eqn{K} clusters.
#'    \item classification - The classification corresponding to z, i.e. map(z).
#'    \item uncertainty - The uncertainty of each CpG site's clustering.    }
#'
#' @examples
#' my.seed <- 190
#' M <- 3
#' data_output <- beta_kn(pca.methylation.data[1:30,2:5], M,
#'                        parallel_process = FALSE, seed = my.seed)
#' thresholds <- threshold(data_output, pca.methylation.data[1:30,2:5], "KN.")
#' @importFrom foreach %dopar%
#' @importFrom stats C
#' @importFrom utils txtProgressBar

beta_kn<-function(data,M=3,parallel_process=FALSE,seed=NULL){

  X=as.data.frame(data)
  ##### KN. Model #######


  ## select the # of cores on which the parallel code is to run
  register=NULL
  if(is.null(register)){


    if (parallel_process == FALSE) {
      # use 2 cores in CRAN/Travis/AppVeyor
      ncores <- 2L
    } else {
      # use all cores in devtools::test()
      ncores <- parallel::detectCores()
    }
    #ncores = parallel::detectCores()
    my.cluster <- parallel::makeCluster(ncores-1)
    doParallel::registerDoParallel(cl = my.cluster)}


  ## declare aliases for dopar command
  `%dopar%` <- foreach::`%dopar%`
  `%do%` <- foreach::`%do%`

  X<-as.data.frame(X)
  if(any(is.na(X)))
  {
    stop("Missing values observed in dataset")
  }

  flag_uni=TRUE

  while(flag_uni){

    ## set the seed for reproducible work
    if (!is.null(seed))
      set.seed(seed)


    flag_uni=FALSE
    K=M
    ## Initial clustering using k-means
    k_cluster<-stats::kmeans(X,K)
    mem <- k_cluster$cluster
    data_clust<-cbind(X,mem)

    ## Get values for parameters C (no. of CpG sites), N (No. of patients)
    x=as.data.frame(data_clust)
    mem=x[,ncol(x)]
    data_full=x
    x=x[,-ncol(x)]
    if(is.vector(x)){
      C=length(x)
      N=1

    }else{
      C=nrow(x)
      N=ncol(x)
    }

    ##  starting values from initial clustering
    mu=matrix(NA,ncol=N,nrow = K)
    sigma=matrix(NA,ncol=N,nrow = K)
    alpha=matrix(NA,ncol=N,nrow = K)
    delta=matrix(NA,ncol=N,nrow = K)
    term=matrix(NA,ncol=N,nrow = K)
    tau=vector(length=K)
    mu_new=matrix(NA,ncol=N,nrow = K)
    sigma_sq=matrix(NA,ncol=N,nrow = K)
    alpha_new=matrix(NA,ncol=N,nrow = K)
    delta_new=matrix(NA,ncol=N,nrow = K)
    term_new=matrix(NA,ncol=N,nrow = K)
    tau_new=vector(length=K)

    x=as.data.frame(x)

    ## function for calculating starting values for beta distribution parameters
    ## we consider the models unimodal in nature which restricts the parameters
    ## to be greater than 1.
    ini_theta_estimation = function(data,N,mem,C)
    {
      data=as.data.frame(data)
      mu=vector()
      sigma=vector()
      term=vector()
      al=vector()
      de=vector()
      tau=vector()
      for(n in 1:N)
      {
        mu[n] <- mean(data[,n])
        sigma[n] <- stats::sd(data[,n])
        term[n]=(mu[n]*(1-mu[n])/(sigma[n]^2))-1
        al[n]=mu[n]*term[n]
        if(al[n]<1)
        {
          al[n]=2
        }
        de[n]=(1-mu[n])*term[n]
        if(de[n]<1)
        {
          de[n]=2
        }
        tau <- length(mem)/C
      }
      return(list(al=al,de=de,tau=tau))
    }

    ## parallelization of starting value calculation
    ini_theta = foreach::foreach(k = 1:K, .combine=rbind) %dopar%
      {ini_theta_estimation(x[mem==k,],N,mem[mem==k],C)}

    ## unlisting the values returned from parallelized code
    mid=K+1
    len_theta=K*2
    tau_len=len_theta+1
    end=length(ini_theta)
    alpha=matrix(unlist(ini_theta[1:K]),nrow=K,ncol=N,byrow=T)
    delta=matrix(unlist(ini_theta[mid:len_theta]),nrow=K,ncol=N,byrow=T)
    tau=unlist(ini_theta[tau_len:end])


    ##Initialise complete data log-likelihood
    z <- matrix(NA,ncol = K,nrow=C)
    z_2 <- matrix(NA,ncol = K,nrow=C)
    z_new<-matrix(NA,ncol=K,nrow=C)
    flag=TRUE
    tol=1e-5
    l0=1e-7
    llk_iter<-vector(mode="numeric")
    llk_iter[1]=l0

    while(flag)
    {
      #E-step
      z_estimation = function(x,al,de,p,N)
      {
        l3=1
        for(n in 1:N) ##no. of samples
        {

          l3=l3*stats::dbeta(x[,n],al[n],de[n])
        }

        return(p*l3)
      }

      ## parallelization of E-step
      z = foreach::foreach(k = 1:K, .combine=cbind) %dopar%
        { z_estimation(x,alpha[k,],delta[k,],tau[k],N)}

      z <- sweep(z, 1, STATS = rowSums(z), FUN = "/")

      ## M-Step

      theta_estimation = function(x,z1,N)
      {
        y11=vector()
        y22=vector()
        al_new=vector()
        de_new=vector()
        for(n in 1:N)
        {
          y11[n]=(sum(z1*log(x[,n])))/(sum(z1))
          y22[n]=(sum(z1*log(1-x[,n])))/(sum(z1))
          term=0
          term=((exp(-y11[n])-1)*(exp(-y22[n])-1))-1
          al_new[n]=0.5+(0.5*exp(-y22[n])/term)
          de_new[n]= (0.5*exp(-y22[n])*(exp(-y11[n])-1))/term
        }
        return(list(al_new=al_new,de_new=de_new))
      }

      ## parallelization of M-step
      theta = foreach::foreach(k = 1:K, .combine=rbind) %dopar%
        { theta_estimation(x,z[,k],N)}

      ## unlisting the values returned from parallelized code
      mid=K+1
      len_theta=K*2
      alpha_new=matrix(unlist(theta[1:K]),nrow=K,ncol=N,byrow=T)
      delta_new=matrix(unlist(theta[mid:len_theta]),nrow=K,ncol=N,byrow=T)


      ### Update z using new alpha and delta

      z_new = foreach::foreach(k = 1:K, .combine=cbind) %dopar%
        {z_estimation(x,alpha_new[k,],delta_new[k,],tau[k],N)}


      z.total=rowSums(z_new)
      z_new <- sweep(z_new, 1, STATS = z.total, FUN = "/")
      s=apply(z_new,2,sum)
      tau_new=s/C

      alpha=alpha_new
      delta=delta_new
      tau=tau_new

      if(any(alpha <1)){
        flag_uni=TRUE
      }else if(any(delta <1)){
        flag_uni=TRUE
      }else{
        flag_uni=FALSE
      }

      if(flag_uni==TRUE){
        if(!is.null(seed))
        {seed = seed+1}else seed=1
        break
      }

      l1=sum(log(z.total))
      llk_iter=c(llk_iter,l1)

      ##check convergence
      diff=abs(l1-l0)/abs(l1)
      flag<- diff>tol
      l0=l1

    }

  }

  ## Clustering
  mem_final<-matrix(NA,C,1)
  mem_final<-apply(z_new, 1, which.max)
  classification=mem_final
  cluster_count=table(mem_final)

  ### uncertainty
  cert=apply(z_new,1,max)
  uc=1-cert

  parallel::stopCluster(cl=my.cluster)
  #### Return data
  return(list(cluster_size=cluster_count,llk=llk_iter,alpha=alpha,delta=delta,
              tau=tau,z=z_new,classification=classification,uncertainty=uc))

}
