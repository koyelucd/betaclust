#' @title Thresholds under the K.. and the KN. models
#' @description An objective method to calculate the threshold points for the clustering solution of the K.. and the KN. models.
#' @details As the K.. model constrains the shape parameters to be equal for all patients, a single pair of threshold points are calculated for all patients. The KN. model allows patient-specific shape parameters which results in a pair of threshold points for each patient based on the shape parameters for that patient.
#' The first threshold point means that any CpG site with beta value less than this value is likely to be hypomethylated.
#' The second threshold point means that any CpG site with beta value greater than this is highly likely to be hypermethylated.
#' A CpG site with beta value lying between the two threshold points is likely to be hemimethylated.
#' @export
#' @param object A \code{\link[betaclust:beta_k]{beta_k}} or \code{\link[betaclust:beta_kn]{beta_kn}} object.
#' @param data A dataframe of dimension \eqn{C \times NR} containing methylation values for \eqn{C} CpG sites from \eqn{R} sample types collected from \eqn{N} patients which was passed as an argument to the \code{\link[betaclust:betaclust]{betaclust}} function.
#' @param model_name The name of the model for which the thresholds need to be calculated.
#' @return thresholds - The threshold points calculated for the selected model. A vector containing two threshold points are returned for the K.. model whereas a matrix containing two threshold points for each patient is returned for the KN. model.
#'
#' @importFrom stats C
#' @seealso \code{\link{beta_k}}
#' @seealso \code{\link{beta_kn}}
#' @seealso \code{\link{betaclust}}


threshold <- function(object,data,model_name){
  threshold_func_low<-function(data_x,alpha,delta,tau,i)
  {
    num<-tau[i]*stats::dbeta(data_x,alpha[i],delta[i])
    deno<-0
    for(j in 1:length(alpha))
    {
      if(j!=i)
      {

        deno<-deno+(tau[j]*stats::dbeta(data_x,alpha[j],delta[j]))
      }
    }

    r=num/deno
    index=which(r>=1)
    l=data_x[min(index)]
    return(l)
  }
  threshold_func_upper<-function(data_x,alpha,delta,tau,i)
  {
    num<-tau[i]*stats::dbeta(data_x,alpha[i],delta[i])
    deno<-0
    for(j in 1:length(alpha))
    {
      if(j!=i)
      {
        deno<-deno+(tau[j]*stats::dbeta(data_x,alpha[j],delta[j]))
      }
    }

    r=num/deno
    index=which(r>=1)
    u=data_x[max(index)]
    return(u)
  }

  if(model_name == "K..")
  {
    data_x=sort(data[,1])
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

    th_new<-vector(length = 2)
    col_th<-vector(length = 1)
    for(i in 1:(ncol(object$alpha)))
    {
      data_x=sort(data[,i])
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
      col_th<-c(col_th,paste0("Patient ",i))
    }
    th_new1<-th_new[,-1]
    th_new1<-as.data.frame(th_new1)
    col_th<-col_th[-1]
    colnames(th_new1)<-col_th
  }
  return(list(threholds=th_new1))
}
