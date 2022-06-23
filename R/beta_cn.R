#' @title The CN. model
#' @description Fit the CN. model from the family of beta mixture models for DNA methylation data.
#'              The CN. model analyses a single DNA sample and identifies the thresholds for the different methylation profiles.
#'
#' @export
#'
#' @details This model clusters each of the \eqn{C} CpG sites into one of \eqn{K} methylation profiles, based on data from \eqn{N} patients for one DNA sample (i.e. \eqn{R=1}).
#' As each CpG site can belong to either of the \eqn{M=3} methylation profiles (hypomethylated, hemimethylated or hypermethylated), the default value of \eqn{K=M=3}.
#' The CN. model differs from the C.. model as it is less parsimonious, allowing cluster and patient-specific shape parameters.
#'
#' @seealso \code{\link{beta_c}}
#' @seealso \code{\link{betaclust}}
#'
#' @param data Methylation values for \eqn{C} CpG sites from \eqn{R=1} samples collected from \eqn{N} patients.
#' @param M Number of methylation profiles to be identified in a DNA sample.
#' @param seed Seed to allow for reproducibility.
#' @param register Setting for registering the parallel backend with the 'foreach' package. To start parallel execution of R code on machine with multiple cores, 'NULL' value needs to be assigned to this parameter.
#'
#' @return A list containing:
#' \itemize{
#'    \item cluster_size - the total number of CpG sites identified in each cluster.
#'    \item llk - a vector containing the log-likelihood value at each step of the EM algorithm.
#'    \item data - this contains the methylation dataset along with the cluster label for each CpG site.
#'    \item alpha - this contains the shape parameter 1 for the beta mixture model.
#'    \item delta - this contains the shape parameter 2 for the mixture model.
#'    \item tau - the proportion of CpG sites in each cluster.
#'    \item z - a matrix containing the probability for each CpG site of belonging to each of the \eqn{K} clusters.
#'    \item uncertainty - the uncertainty of each CpG site's clustering.    }
#'
#' @examples
#' \dontrun{
#' data(pca.methylation.data)
#' my.seed=190
#' M=3
#' data_output=beta_cn(pca.methylation.data[,2:5],M,seed=my.seed)
#' }
#' @importFrom foreach %dopar%
#' @importFrom stats C
#' @importFrom utils txtProgressBar
#' @references {Microsoft, Weston, S. (2022): foreach: Provides Foreach Looping Construct. R package version 1.5.2. https://CRAN.R-project.org/package=foreach.}

beta_cn<-function(data,M=3,seed,register=NULL){

  X=data
  ##### CN Model #######
  ## set the seed for reproducible work
  if (!missing(seed))
    set.seed(seed)

  ## select the # of cores on which the parallel code is to run
  if(is.null(register)){
    ncores = parallel::detectCores()
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

    l1=sum(log(z.total))
    llk_iter=c(llk_iter,l1)

    ##check convergence
    diff=abs(l1-l0)/abs(l1)
    flag<- diff>tol
    l0=l1

  }

  ## Clustering
  mem_final<-matrix(NA,C,1)
  complete_data<-matrix(NA,C,(N+1))
  mem_final<-apply(z_new, 1, which.max)
  complete_data<-cbind(x,mem_final)
  cluster_count=table(mem_final)

  ### uncertainty
  cert=apply(z_new,1,max)
  uc=1-cert


  parallel::stopCluster(cl=my.cluster)
  #### Return data
  return(list(cluster_size=cluster_count,llk=llk_iter,data=complete_data,alpha=alpha,delta=delta,tau=tau,z=z_new,uncertainty=uc))

}
