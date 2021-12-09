#' CN. Model from the family of beta mixture models for DNA methylation data
#' @export
#' @param  X methylation values for CpG sites frpm R samples collected from N patients
#' @param K number of methylation groups to be identified (default=3)
#' @param seed seed for reproducible work
#' @param register setting for parallelization
#' @importFrom foreach %dopar%
#' @importFrom stats C
#' @importFrom utils txtProgressBar

beta_cn<-function(data,K=3,seed,register=NULL){

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
  beta=matrix(NA,ncol=N,nrow = K)
  term=matrix(NA,ncol=N,nrow = K)
  tau=vector(length=K)
  mu_new=matrix(NA,ncol=N,nrow = K)
  sigma_sq=matrix(NA,ncol=N,nrow = K)
  alpha_new=matrix(NA,ncol=N,nrow = K)
  beta_new=matrix(NA,ncol=N,nrow = K)
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
    be=vector()
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
      be[n]=(1-mu[n])*term[n]
      if(be[n]<1)
      {
        be[n]=2
      }
      tau <- length(mem)/C
    }
    return(list(al=al,be=be,tau=tau))
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
  beta=matrix(unlist(ini_theta[mid:len_theta]),nrow=K,ncol=N,byrow=T)
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
    z_estimation = function(x,al,be,p,N)
    {
      l3=1
      for(n in 1:N) ##no. of samples
      {

        l3=l3*stats::dbeta(x[,n],al[n],be[n])
      }

      return(p*l3)
    }

    ## parallelization of E-step
    z = foreach::foreach(k = 1:K, .combine=cbind) %dopar%
     { z_estimation(x,alpha[k,],beta[k,],tau[k],N)}

    z <- sweep(z, 1, STATS = rowSums(z), FUN = "/")

    ## M-Step

    theta_estimation = function(x,z1,N)
    {
      y11=vector()
      y22=vector()
      al_new=vector()
      be_new=vector()
      for(n in 1:N)
      {
        y11[n]=(sum(z1*log(x[,n])))/(sum(z1))
        y22[n]=(sum(z1*log(1-x[,n])))/(sum(z1))
        term=0
        term=((exp(-y11[n])-1)*(exp(-y22[n])-1))-1
        al_new[n]=0.5+(0.5*exp(-y22[n])/term)
        be_new[n]= (0.5*exp(-y22[n])*(exp(-y11[n])-1))/term
      }
      return(list(al_new=al_new,be_new=be_new))
    }

    ## parallelization of M-step
    theta = foreach::foreach(k = 1:K, .combine=rbind) %dopar%
      { theta_estimation(x,z[,k],N)}

    ## unlisting the values returned from parallelized code
    mid=K+1
    len_theta=K*2
    alpha_new=matrix(unlist(theta[1:K]),nrow=K,ncol=N,byrow=T)
    beta_new=matrix(unlist(theta[mid:len_theta]),nrow=K,ncol=N,byrow=T)


    ### Update z using new alpha and beta

    z_new = foreach::foreach(k = 1:K, .combine=cbind) %dopar%
      {z_estimation(x,alpha_new[k,],beta_new[k,],tau[k],N)}


    z.total=rowSums(z_new)
    z_new <- sweep(z_new, 1, STATS = z.total, FUN = "/")
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
  cluster_count=table(mem_final)

  ### uncertainty
  cert=apply(z_new,1,max)
  uc=1-cert


  parallel::stopCluster(cl=my.cluster)
  #### Return data
  return(list(cluster_count=cluster_count,llk=llk_iter,data=complete_data,alpha=alpha,beta=beta,tau=tau,z=z_new,uncertainty=uc))

}
