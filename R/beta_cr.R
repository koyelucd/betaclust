#' @title The C.R Model
#' @export
#' @description A beta mixture model for identifying differentially methylated CpG sites between \eqn{R} DNA samples collected from \eqn{N} patients.
#'
#' @seealso \code{\link{betaclust}}
#'
#' @param data Methylation values for \eqn{C} CpG sites from \eqn{R} samples collected from \eqn{N} patients.
#' @param M Number of methylation profiles to be identified.
#' @param patients Number of patients in the study.
#' @param samples Number of samples collected from each patient for study.
#' @param seed Seed to allow for reproducibility.
#' @param register Setting for registering the parallel backend with the "foreach" package. To start parallel execution of R code on machine with multiple cores, "NULL" value needs to be assigned to this parameter.
#'
#' @details
#' The C.R model allows identification of the differentially methylated CpG sites between the \eqn{R} DNA samples collected from each of \eqn{N} patients.
#' As each CpG site in a DNA sample can belong to either of \eqn{M} methylation profiles, there can be \eqn{K=M^R} methylation profile changes between \eqn{R} DNA samples.
#' The parameters vary for each DNA sample but are constrained to be equal for each patient. An initial clustering using K-means is performed to identify \eqn{K} clusters. The resulting clustering solution is provided as
#' starting values to the Expectation-Maximisation algorithm. A digamma approximation is used to obtain the maximised
#' parameters in the M-step instead of a computationally inefficient numerical optimisation step.
#' @return A list containing:
#' \itemize{
#'    \item cluster_size - the total number of CpG sites identified in each of the K clusters.
#'    \item llk - a vector containing the log-likelihood value at each step of the EM algorithm.
#'    \item data - this contains the methylation dataset along with the cluster label for each CpG site.
#'    \item alpha - this contains the shape parameter 1 for the beta mixture model.
#'    \item delta - this contains the shape parameter 2 for the beta mixture model.
#'    \item tau - the proportion of CpG sites in each cluster.
#'    \item z - a matrix containing the probability for each CpG site of belonging to each of the \eqn{K} clusters.
#'    \item uncertainty - the uncertainty of each CpG site's clustering.    }
#'
#' @examples
#' \dontrun{
#' data(pca.methylation.data)
#' my.seed=190
#' M=3
#' patients=4
#' samples=2
#' data_output=beta_cr(pca.methylation.data[,2:5],M,patients,samples,seed=my.seed)
#' }
#' @importFrom foreach %dopar%
#' @importFrom stats C
#' @importFrom utils txtProgressBar
#' @references {Microsoft, Weston, S. (2022): foreach: Provides Foreach Looping Construct. R package version 1.5.2. https://CRAN.R-project.org/package=foreach.}




beta_cr<-function(data,M=3,patients,samples,seed,register=NULL){

  X=data
  ##### CR Model #######

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

  ## Check for missing data
  X<-as.data.frame(X)
  if(any(is.na(X)))
  {
    stop("Missing values observed in dataset")
  }
  ## Initial clustering using k-means
  ## K= 3 for hypo,hyper and hemi methylation. For R sample types there can be
  ## 3^R combinations of these profiles hence K^R clusters

  clusters=M^samples
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
  delta=matrix(NA,ncol=R,nrow = K)
  term=matrix(NA,ncol=R,nrow = K)
  tau=vector(length=K)
  alpha_new=matrix(NA,ncol=R,nrow = K)
  delta_new=matrix(NA,ncol=R,nrow = K)
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
    de=vector()
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
      de[r]=(1-mu[r])*term[r]
      if(de[r]<1)
      {
        de[r]=2
      }
      tau <- length(mem)/C
    }
    return(list(al=al,de=de,tau=tau))
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
  delta=matrix(unlist(ini_theta[mid:len_theta]),nrow=K,ncol=R,byrow=T)
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

    z_estimation = function(x,al,de,p,R,N)
    {
      l3=1
      for(r in 1:R) ##no. of samples
      {

        l3=l3*stats::dbeta(x[,(((r-1)*N)+1):(((r-1)*N)+N)],al[r],de[r])
      }

      l4=apply(l3,1,prod)
      return(p*l4)
    }

    ## parallelization of E-step
    z = foreach::foreach(k = 1:K, .combine=cbind) %dopar%
      z_estimation(x,alpha[k,],delta[k,],tau[k],R,N)

    z <- sweep(z, 1, STATS = rowSums(z), FUN = "/")


    ## M-step
    theta_estimation = function(x,z1,R,N)
    {
      y11=vector()
      y22=vector()
      al_new=vector()
      de_new=vector()
      for(r in 1:R)
      {
        y11[r]=(sum(z1*rowSums(log(x[,(((r-1)*N)+1):(((r-1)*N)+N)]))))/(N*sum(z1))
        y22[r]=(sum(z1*rowSums(log(1-x[,(((r-1)*N)+1):(((r-1)*N)+N)]))))/(N*sum(z1))
        term=0
        term=((exp(-y11[r])-1)*(exp(-y22[r])-1))-1
        al_new[r]=0.5+(0.5*exp(-y22[r])/term)
        de_new[r]= (0.5*exp(-y22[r])*(exp(-y11[r])-1))/term
      }
      return(list(al_new=al_new,de_new=de_new))
    }

    ## parallelization of M-step
    theta = foreach::foreach(k = 1:K, .combine=rbind) %dopar%
      theta_estimation(x,z[,k],R,N)

    ## unlisting the values returned from parallelized code
    mid=K+1
    len_theta=K*2
    alpha_new=matrix(unlist(theta[1:K]),nrow=K,ncol=R,byrow=T)
    delta_new=matrix(unlist(theta[mid:len_theta]),nrow=K,ncol=R,byrow=T)


    ### Update z using new alpha and delta

    z_new = foreach::foreach(k = 1:K, .combine=cbind) %dopar%
      z_estimation(x,alpha_new[k,],delta_new[k,],tau[k],R,N)

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
  complete_data<-matrix(NA,C,(N*R+1))
  mem_final<-apply(z_new, 1, which.max)
  complete_data<-cbind(x,mem_final)
 # cluster_count=table(mem_final)


  ### uncertainty
  #cert=apply(z_new,1,max)
  #uc=1-cert

  parallel::stopCluster(cl=my.cluster)

  #### Sorting the clusters as per interest
  mean_delta=as.data.frame(alpha/(alpha+delta))
  mean_delta$diff<-abs(mean_delta$V1-mean_delta$V2)
  mean_sorted<-mean_delta[order(mean_delta$diff,decreasing = T),]
  mean_row<-row.names(mean_sorted)
  data_final=as.data.frame(complete_data)
  data_final$mem_final<-as.numeric(data_final$mem_final)
  data_final$new_mem_final<-NA
  mean_row<-as.numeric(mean_row)
  mean_row<-cbind(mean_row,seq(1,K,by=1))
  alpha_final=matrix(NA,nrow=K,ncol=R)
  delta_final=matrix(NA,nrow=K,ncol=R)
  tau_final=vector(length=K)
  z_final=matrix(NA,nrow=C,ncol=K)
  for(i in 1:K)
  {
    data_final["new_mem_final"][data_final["mem_final"]==mean_row[i,1]]<-mean_row[i,2]
    swapped_row=mean_row[i,1]
    alpha_final[i,]=alpha[swapped_row,]
    delta_final[i,]=delta[swapped_row,]
    tau_final[i]=tau[swapped_row]
    z_final[,i]=z_new[,swapped_row]
  }
  data_final<-data_final[,-(N*R+1)]
  colnames(data_final)[(N*R+1)]<-"mem_final"
  data_final$mem_final<-as.factor(data_final$mem_final)
  cert_final=apply(z_final,1,max)
  uc_final=1-cert_final
  cluster_count=table(data_final$mem_final)

  #### Return data
  #return(list(cluster_count=cluster_count,llk=llk_iter,data=complete_data,alpha=alpha,delta=delta,tau=tau,z=z_new,uncertainty=uc))
  return(list(cluster_size=cluster_count,llk=llk_iter,data=data_final,
              alpha=alpha_final,delta=delta_final,tau=tau_final,
              z=z_final,uncertainty=uc_final))

}
