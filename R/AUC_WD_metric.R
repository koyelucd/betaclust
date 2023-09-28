#' @title AUC and WD function
#' @description Function to find the level of similarities between the \eqn{R} cumulative distributions estimated in each of the \eqn{K} clusters.
#' @details Function to find the level of similarities between the \eqn{R} cumulative distributions estimated in each of the \eqn{K} clusters.
#' @export
#' @seealso \code{\link{betaclust}}
#' @param alpha The first shape parameter estimated for the beta mixture model.
#' @param delta The second shape parameter estimated for the beta mixture model.
#' @param K The number of clusters estimated.
#' @param R The number of sample types in the dataset.
#' @return The list with AUC and WD values.
#' @importFrom pROC auc
#' @importFrom stats rbeta pbeta

#AUC_WD_metric<-function(object,K,R){
  AUC_WD_metric<-function(alpha,delta,K,R){

  comb=factorial(R) / (factorial(2) * factorial(R - 2))
  n <- 1000
  auc_mat=matrix(NA,nrow=K,ncol=comb)
  d_mat=matrix(NA,nrow=K,ncol=comb)
  L <- 1001
  d_vec_sim=vector()
  z <- seq(0, 1, length = L)
  text=vector()
  for(j in 1:K)
  {
    shape_1 <- alpha[j, 1:R]
    shape_2 <- delta[j, 1:R]
    vec_count=1
    for(i in 1: (R-1))
    {
      for(k in (i+1):R)
      {
        group_1 <- stats::rbeta(n, shape1 = shape_1[i], shape2 = shape_2[i])
        group_2 <- stats::rbeta(n, shape1 = shape_1[k], shape2 = shape_2[k])
        auc_dat <- data.frame(predictor=c(group_1, group_2),
                              response=factor(c(rep(0,n), rep(1,n))))
        auc_value <- pROC::auc(predictor = auc_dat$predictor,
                         response  = auc_dat$response,quiet=TRUE)
        auc_mat[j,vec_count]=unlist(auc_value)
        d <- sum(abs(stats::pbeta(z, shape_1[i], shape_2[i])-
                       stats::pbeta(z, shape_1[k], shape_2[k]))) / L
        d_mat[j,vec_count]=d
        if(j==1)
        {
          text[vec_count]=paste(i,"->",k)
        }
        vec_count=vec_count+1
      }
    }

  }
  auc_mat=as.data.frame(auc_mat)
  colnames(auc_mat)=text
  d_mat=as.data.frame(d_mat)
  colnames(d_mat)=text

  out=list(AUC=auc_mat, WD=d_mat)
  return(out)
}
