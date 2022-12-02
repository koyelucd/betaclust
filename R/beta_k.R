globalVariables(c("k"))
#' @title Fit the K.. model
#' @description Fit the K.. model from the \code{\link[betaclust:betaclust]{betaclust}} family of beta mixture models for DNA methylation data.
#'              The K.. model analyses a single DNA sample and identifies the thresholds between the different methylation states.
#'
#' @export
#'
#' @details The K.. model clusters each of the \eqn{C} CpG sites into one of \eqn{K} methylation states, based on data from \eqn{N} patients for one DNA sample (i.e. \eqn{R = 1}).
#' As each CpG site can belong to any of the \eqn{M = 3} methylation states (hypomethylated, hemimethylated or hypermethylated), the default value of \eqn{K = M = 3}.
#' Under the K.. model the shape parameters of each cluster are constrained to be equal for each patient. The returned object from this function can be passed as an input parameter to the
#' \code{\link[betaclust:threshold]{threshold}} function available in this package to calculate the thresholds between the methylation states.
#'
#' @seealso \code{\link{beta_kn}}
#' @seealso \code{\link{betaclust}}
#' @seealso \code{\link{threshold}}
#'
#' @param data A dataframe of dimension \eqn{C \times N} containing methylation values for \eqn{C} CpG sites from \eqn{R = 1} samples collected from \eqn{N} patients.
#' Samples are grouped together in the dataframe such that the columns are ordered as Sample1_Patient1, Sample1_Patient2, etc.
#' @param M Number of methylation states to be identified in a DNA sample.
#' @param parallel_process The "TRUE" option results in parallel processing of the models for increased computational efficiency. The default option has been set as "FALSE" due to package testing limitations.
#' @param seed Seed to allow for reproducibility (default = NULL).
#'
#'
#' @return A list containing:
#' \itemize{
#'    \item cluster_size - The total number of CpG sites in each of the K clusters.
#'    \item llk - A vector containing the log-likelihood value at each step of the EM algorithm.
#'    \item alpha - The first shape parameter for the beta mixture model.
#'    \item delta - The second shape parameter for the beta mixture model.
#'    \item tau - The proportion of CpG sites in each cluster.
#'    \item z - A matrix of dimension \eqn{C \times K} containing the posterior probability of each CpG site belonging to each of the \eqn{K} clusters.
#'    \item classification - The classification corresponding to z, i.e. map(z).
#'    \item uncertainty - The uncertainty of each CpG site's clustering.    }
#'
#' @examples
#' my.seed <- 190
#' M <- 3
#' data_output <- beta_k(pca.methylation.data[1:30,2:5], M,
#'                       parallel_process = FALSE, seed = my.seed)
#' thresholds <- threshold(data_output, pca.methylation.data[1:30,2:5], "K..")
#' @importFrom foreach %dopar%
#' @importFrom stats C
#' @importFrom utils txtProgressBar


beta_k<-function(data,M=3,parallel_process=FALSE,seed = NULL){


  ######### K.. Model ##########
  X=data


  ## select the # of cores on which the parallel code is to run
  register=NULL
  if(is.null(register)){

    if (parallel_process==FALSE) {
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
    x=as.matrix(data_clust)
    mem=x[,ncol(x)]
    data_full=x
    x=x[,-ncol(x)]
    C=nrow(x)
    N=ncol(x)

    ##  starting values from initial clustering
    mu=vector("numeric",K)
    sigma=vector("numeric",K)
    mu_new=vector("numeric",K)
    sigma_sq=vector("numeric",K)
    alpha=vector("numeric",K)
    delta=vector("numeric",K)
    alpha_new=vector("numeric",K)
    delta_new=vector("numeric",K)
    term=vector("numeric",K)
    term_new=vector("numeric",K)
    tau=vector(length=K)
    tau_new=vector(length=K)


    ## function for calculating starting values for beta distribution parameters
    ## we consider the models unimodal in nature which restricts the parameters
    ## to be greater than 1.
    ini_theta_estimation = function(data,N,mem,C)
    {
      mu=vector()
      sigma=vector()
      term=vector()
      al=vector()
      de=vector()
      pi=vector()
      mu<-mean(data)
      sigma <- stats::sd(data)
      term=(mu*(1-mu)/(sigma^2))-1
      al=mu*term
      if(al<=1)
      {
        al=2
      }
      de=(1-mu)*term
      if(de<=1)
      {
        de=2
      }
      pi <- length(mem)/C
      return(list(al=al,de=de,tau=pi))
    }

    ## parallelization of starting value calculation
    ini_theta = foreach::foreach(k = 1:K, .combine=rbind) %dopar%
      {ini_theta_estimation(x[mem==k,],N,mem[mem==k],C)}

    ## unlisting the values returned from parallelized code
    mid=K+1
    len_theta=K*2
    tau_len=len_theta+1
    end=length(ini_theta)
    alpha=unlist(ini_theta[1:K])
    delta=unlist(ini_theta[mid:len_theta])
    tau=unlist(ini_theta[tau_len:end])


    ## Initialise complete data log-likelihood
    z <- matrix(NA,ncol = K,nrow=C)
    z_new <- matrix(NA,ncol = K,nrow=C)
    tol=1e-04
    l0=1e-07
    llk_iter<-vector(mode="numeric")
    llk_iter[1]=l0
    flag=TRUE

    while(flag)
    {
      ## E-step
      z_estimation = function(x,al,de,pi,N)
      {

        l3=1
        for(n in 1:N) ##no. of samples
        {

          l3=l3*stats::dbeta(x[,n],al,de)
        }

        return(pi*l3)
      }

      ## parallelization of E-step
      z = foreach::foreach(k = 1:K, .combine=cbind) %dopar%
        { z_estimation(x,alpha[k],delta[k],tau[k],N)}

      z <- sweep(z, 1, STATS = rowSums(z), FUN = "/")

      ## M-Step

      theta_estimation = function(x,z1,N)
      {
        y11=0
        y22=0
        al_new=0
        de_new=0
        y11=(sum(z1*rowSums(log(x))))/(N*sum(z1))
        y22=(sum(z1*rowSums(log(1-x))))/(N*sum(z1))
        term=0
        term=((exp(-y11)-1)*(exp(-y22)-1))-1
        al_new=0.5+(0.5*exp(-y22)/term)
        de_new= (0.5*exp(-y22)*(exp(-y11)-1))/term
        return(list(al_new=al_new,de_new=de_new))
      }


      ## parallelization of M-step
      theta = foreach::foreach(k = 1:K, .combine=rbind) %dopar%
        {theta_estimation(x,z[,k],N)}

      ## unlisting the values returned from parallelized code
      mid=K+1
      len_theta=K*2
      alpha_new=unlist(theta[1:K])
      delta_new=unlist(theta[mid:len_theta])


      ### Update z using new alpha and delta

      z_new = foreach::foreach(k = 1:K, .combine=cbind) %dopar%
        {z_estimation(x,alpha_new[k],delta_new[k],tau[k],N)}


      z.total=rowSums(z_new)
      z_new=z_new/z.total
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
  #complete_data<-matrix(NA,C,(N+1))
  mem_final<-apply(z_new, 1, which.max)
  classification=mem_final
  #complete_data<-cbind(x,mem_final)
  cluster_count=table(mem_final)

  ### Uncertainty
  cert=apply(z_new,1,max)
  uc=1-cert

  tau=round((as.numeric(cluster_count)/C),3)

  parallel::stopCluster(cl=my.cluster)

  #### Return data
  return(list(cluster_size=cluster_count,llk=llk_iter,alpha=alpha,delta=delta,tau=tau,z=z_new,classification=classification,uncertainty=uc))

}
