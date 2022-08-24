#' @title Thresholds for K.. and KN. models
#' @description An empirical cumulative distribution function (ECDF) plot for a betaclust object.
#' @details This function plots the ECDF graphs of the differentially methylated CpG sites identified using the K.R model for all patient samples.
#' The graph visualises the methylation profile changes between the different DNA samples for each patient.
#' @export
#' @param object A dataframe containing methylation values of identified differentially methylated regions related to a gene.
#' @return The thresholds calculated for the clustering solution
#'
#' @seealso \code{\link{betaclust}}
#' @seealso \code{\link{beta_k}}
#' @seealso \code{\link{beta_kn}}


threshold <- function(object,model_name){
  threshold_func_low<-function(data_x,alpha,delta,tau,i)
  {
    num<-tau[i]*dbeta(data_x,alpha[i],delta[i])
    deno<-0
    for(j in 1:length(alpha))
    {
      if(j!=i)
      {

        deno<-deno+(tau[j]*dbeta(data_x,alpha[j],delta[j]))
      }
    }

    r=num/deno
    index=which(r>=1)
    l=data_x[min(index)]
    return(l)
  }
  threshold_func_upper<-function(data_x,alpha,delta,tau,i)
  {
    num<-tau[i]*dbeta(data_x,alpha[i],delta[i])
    deno<-0
    for(j in 1:length(alpha))
    {
      if(j!=i)
      {
        deno<-deno+(tau[j]*dbeta(data_x,alpha[j],delta[j]))
      }
    }

    r=num/deno
    index=which(r>=1)
    u=data_x[max(index)]
    return(u)
  }

  if(model_name == "K..")
  {
    data_x=sort(object$data[,1])
    mode<-(object$alpha-1)/(object$alpha+object$delta-2)
    cluster<-c(1,2,3)
    mode_vec<-cbind(mode,cluster)
    mode_vec<-as.data.frame(mode_vec)
    mode_vec<-mode_vec[order(mode_vec$mode),]
    hypo<-as.numeric(mode_vec[1,2])
    hyper<-as.numeric(mode_vec[3,2])
    th_1<-threshold_func_upper(data_x,object$alpha,object$delta,object$tau,hypo)
    th_2<-threshold_func_low(data_x,object$alpha,object$delta,object$tau,hyper)
    th_vec<-c(th_1,th_2)
    th_new1<-unique(round(th_vec,3))
  }else{
    for(i in 1:object$K)
    {
      th_new<-vector(length = 2)
      data_x=sort(object$data[,i])
      mode<-(object$alpha[,i]-1)/(object$alpha[,i]+object$delta[,i]-2)
      cluster<-c(1,2,3)
      mode_vec<-cbind(mode,cluster)
      mode_vec<-as.data.frame(mode_vec)
      mode_vec<-mode_vec[order(mode_vec$mode),]
      hypo<-as.numeric(mode_vec[1,2])
      hyper<-as.numeric(mode_vec[3,2])
      th_1<-threshold_func_upper(data_x,object$alpha[,i],object$delta[,i],object$tau,hypo)
      th_2<-threshold_func_low(data_x,object$alpha[,i],object$delta[,i],object$tau,hyper)
      th_vec<-c(th_1,th_2)
      th_new1<-unique(round(th_vec,3))
      th_new<-cbind(th_new,th_new1)
    }
    th_new1<-th_new[,-1]
  }
return(th_new1)
}
