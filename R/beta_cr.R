#' C.R Model from the family of beta mixture models for DNA methylation data
#' @export
#' @param  X methylation values for CpG sites frpm R samples collected from N patients
#' @param K number of methylation groups to be identified (default=3)
#' @param patients number of patients in the study
#' @param samples number of samples collected from each patient for study
#' @param seed seed for reproducible work
#' @param register setting for parallelization
#' @importFrom foreach %dopar%
#' @importFrom stats C


beta_cr<-function(X,K=3,patients,samples,seed,register=NULL){

  ##### CR Model #######

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

  ## Check for missing data
  X<-as.data.frame(X)
  if(any(is.na(X)))
  {
    stop("Missing values observed in dataset")
  }
  ## Initial clustering using k-means
  ## K= 3 for hypo,hyper and hemi methylation. For R sample types there can be
  ## 3^R combinations of these profiles hence K^R clusters

  clusters=K^samples
  K=clusters
  k_cluster<-stats::kmeans(X,K)
  mem <- k_cluster$cluster
  data_clust<-cbind(X,mem)

  ## Get values for parameters C (no. of CpG sites), N (No. of patients), R (No. of samples)
  x=as.matrix(data_clust)
  mem=x[,ncol(x)]
  data_full=x
  x=x[,-ncol(x)]
  C=nrow(x)
  R=samples
  N=patients

  ##  starting values from initial clustering
  mu=matrix(NA,ncol=R,nrow = K)
  sigma=matrix(NA,ncol=R,nrow = K)
  alpha=matrix(NA,ncol=R,nrow = K)
  beta=matrix(NA,ncol=R,nrow = K)
  term=matrix(NA,ncol=R,nrow = K)
  tau=vector(length=K)
  alpha_new=matrix(NA,ncol=R,nrow = K)
  beta_new=matrix(NA,ncol=R,nrow = K)
  tau_new=vector(length=K)


  ## function for calculating starting values for beta distribution parameters
  ## we consider the models unimodal in nature which restricts the parameters
  ## to be greater than 1.
  ini_theta_estimation = function(data,R,N,mem,C)
  {
    mu=vector()
    sigma=vector()
    term=vector()
    al=vector()
    be=vector()
    tau=vector()
    for(r in 1:R)
    {
      mu[r] <- mean(data[,(((r-1)*N)+1):(((r-1)*N)+N)])
      sigma[r] <- stats::sd(data[,(((r-1)*N)+1):(((r-1)*N)+N)])
      term[r]=(mu[r]*(1-mu[r])/(sigma[r]^2))-1
      al[r]=mu[r]*term[r]
      if(al[r]<1)
      {
        al[r]=2
      }
      be[r]=(1-mu[r])*term[r]
      if(be[r]<1)
      {
        be[r]=2
      }
      tau <- length(mem)/C
    }
    return(list(al=al,be=be,tau=tau))
  }

  ## parallelization of starting value calculation
  ini_theta = foreach::foreach(k = 1:K, .combine=rbind) %dopar%
    ini_theta_estimation(x[mem==k,],R,N,mem[mem==k],C)

  ## unlisting the values returned from parallelized code
  mid=K+1
  len_theta=K*2
  tau_len=len_theta+1
  end=length(ini_theta)
  alpha=matrix(unlist(ini_theta[1:K]),nrow=K,ncol=R,byrow=T)
  beta=matrix(unlist(ini_theta[mid:len_theta]),nrow=K,ncol=R,byrow=T)
  tau=unlist(ini_theta[tau_len:end])




  ##Initialise complete data log-likelihood
  z <- matrix(NA,ncol = K,nrow=C)
  z_new<-matrix(NA,ncol=K,nrow=C)
  flag=TRUE
  tol=1e-5
  l0=1e-7
  llk_iter<-vector(mode="numeric")
  llk_iter[1]=l0
  i=0

  while(flag)
  {
    #E-step

    z_estimation = function(x,al,be,p,R,N)
    {
      l3=1
      for(r in 1:R) ##no. of samples
      {

        l3=l3*stats::dbeta(x[,(((r-1)*N)+1):(((r-1)*N)+N)],al[r],be[r])
      }

      l4=apply(l3,1,prod)
      return(p*l4)
    }

    ## parallelization of E-step
    z = foreach::foreach(k = 1:K, .combine=cbind) %dopar%
      z_estimation(x,alpha[k,],beta[k,],tau[k],R,N)

    z <- sweep(z, 1, STATS = rowSums(z), FUN = "/")


    ## M-step
    theta_estimation = function(x,z1,R,N)
    {
      y11=vector()
      y22=vector()
      al_new=vector()
      be_new=vector()
      for(r in 1:R)
      {
        y11[r]=(sum(z1*rowSums(log(x[,(((r-1)*N)+1):(((r-1)*N)+N)]))))/(N*sum(z1))
        y22[r]=(sum(z1*rowSums(log(1-x[,(((r-1)*N)+1):(((r-1)*N)+N)]))))/(N*sum(z1))
        term=0
        term=((exp(-y11[r])-1)*(exp(-y22[r])-1))-1
        al_new[r]=0.5+(0.5*exp(-y22[r])/term)
        be_new[r]= (0.5*exp(-y22[r])*(exp(-y11[r])-1))/term
      }
      return(list(al_new=al_new,be_new=be_new))
    }

    ## parallelization of M-step
    theta = foreach::foreach(k = 1:K, .combine=rbind) %dopar%
      theta_estimation(x,z[,k],R,N)

    ## unlisting the values returned from parallelized code
    mid=K+1
    len_theta=K*2
    alpha_new=matrix(unlist(theta[1:K]),nrow=K,ncol=R,byrow=T)
    beta_new=matrix(unlist(theta[mid:len_theta]),nrow=K,ncol=R,byrow=T)


    ### Update z using new alpha and beta

    z_new = foreach::foreach(k = 1:K, .combine=cbind) %dopar%
      z_estimation(x,alpha_new[k,],beta_new[k,],tau[k],R,N)

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
  complete_data<-matrix(NA,C,(N*R+1))
  mem_final<-apply(z_new, 1, which.max)
  complete_data<-cbind(x,mem_final)


  ### Uncertainity
  cert=apply(z_new,1,max)
  uc=1-cert


  #### Return data
  return(list(llk=llk_iter,data=complete_data,alpha=alpha,beta=beta,tau=tau,z=z_new,uncertainity=uc))

}
