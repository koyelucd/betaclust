#' C.. Model from the family of beta mixture models for DNA methylation data
#' @export
#' @param  X methylation values for CpG sites frpm R samples collected from N patients
#' @param K number of methylation groups to be identified (default=3)
#' @param seed seed for reproducible work
#' @param register setting for parallelization
#' @importFrom foreach %dopar%
#' @importFrom stats C

beta_c<-function(X,K=3,seed,register=NULL){

  ## set the seed for reproducible work
  if (!missing(seed))
    set.seed(seed)

  ## select the # of cores on which the parallel code is to run
  if(is.null(register)){
    ncores = parallel::detectCores()
    doParallel::registerDoParallel(ncores-1)}

  ## declare aliases for dopar command
  `%dopar%` <- foreach::`%dopar%`
  `%do%` <- foreach::`%do%`

  X<-as.data.frame(X)
  if(any(is.na(X)))
  {
    stop("Missing values observed in dataset")
  }

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
  beta=vector("numeric",K)
  alpha_new=vector("numeric",K)
  beta_new=vector("numeric",K)
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
    be=vector()
    pi=vector()
    mu<-mean(data)
    sigma <- stats::sd(data)
    term=(mu*(1-mu)/(sigma^2))-1
    al=mu*term
    if(al<=1)
    {
      al=2
    }
    be=(1-mu)*term
    if(be<=1)
    {
      be=2
    }
    pi <- length(mem)/C
    return(list(al=al,be=be,tau=pi))
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
  beta=unlist(ini_theta[mid:len_theta])
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
    z_estimation = function(x,al,be,pi,N)
    {

      l3=1
      for(n in 1:N) ##no. of samples
      {

        l3=l3*stats::dbeta(x[,n],al,be)
      }

      return(pi*l3)
    }

    ## parallelization of E-step
    z = foreach::foreach(k = 1:K, .combine=cbind) %dopar%
     { z_estimation(x,alpha[k],beta[k],tau[k],N)}

    z <- sweep(z, 1, STATS = rowSums(z), FUN = "/")

    ## M-Step

    theta_estimation = function(x,z1,N)
    {
      y11=0
      y22=0
      al_new=0
      be_new=0
      y11=(sum(z1*rowSums(log(x))))/(N*sum(z1))
      y22=(sum(z1*rowSums(log(1-x))))/(N*sum(z1))
      term=0
      term=((exp(-y11)-1)*(exp(-y22)-1))-1
      al_new=0.5+(0.5*exp(-y22)/term)
      be_new= (0.5*exp(-y22)*(exp(-y11)-1))/term
      return(list(al_new=al_new,be_new=be_new))
    }


    ## parallelization of M-step
    theta = foreach::foreach(k = 1:K, .combine=rbind) %dopar%
      {theta_estimation(x,z[,k],N)}

    ## unlisting the values returned from parallelized code
    mid=K+1
    len_theta=K*2
    alpha_new=unlist(theta[1:K])
    beta_new=unlist(theta[mid:len_theta])


    ### Update z using new alpha and beta

    z_new = foreach::foreach(k = 1:K, .combine=cbind) %dopar%
      {z_estimation(x,alpha_new[k],beta_new[k],tau[k],N)}


    z.total=rowSums(z_new)
    z_new=z_new/z.total
    s=apply(z_new,2,sum)
    tau_new=s/C

    alpha=alpha_new
    beta=beta_new
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

  ### Uncertainity
  cert=apply(z_new,1,max)
  uc=1-cert


  doParallel::stopImplicitCluster()

  #### Return data
  return(list(data=complete_data,alpha=alpha,beta=beta,tau=tau,z=z_new,uncertainity=uc,llk=llk_iter))

}
