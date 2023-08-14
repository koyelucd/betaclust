library(betaclust)
library(e1071)
library(pROC)

### Simulation study
set.seed(03082023) ## setting seed as a whole before generating random numbers
seed_vec=runif(100,1,999)

## Simulation study for K.. and KN. model
df_threshold_list=list()
beta_threshold_out=list()
ari_threshold=vector()
og_group=vector()
for(i in 1:length(seed_vec))
{
  set.seed(seed_vec[i])
  C=600000
  probs=c(0.35,0.7)
  index=runif(C)
  R=4
  df=matrix(NA,ncol = R,nrow = C)
  og_group=vector(length = C)
  true_dmc=vector(length = C)
  for(c in 1:C)
  {
    if(index[c]<probs[1])
    {
      df[c,]=rbeta(4,2,20) #hypomethylated
      og_group[c]=1
    }else if(index[c]<probs[2])
    {
      df[c,]=rbeta(R,4,3) #hemimethylated
      og_group[c]=2
    }else 
    {
      df[c,]=rbeta(4,20,2) #hypermethylated
      og_group[c]=3
    }
  }
  
  df=cbind(df,og_group)
  df=as.data.frame(df)
  colnames(df)<-c("P_Sam1_1","P_Sam1_2","P_Sam1_3","P_Sam1_4","Original_Group")
  
  df_threshold_list[[i]]=df
  beta_sim<-betaclust(df,3,4,1,c("K..","KN."),parallel_process = TRUE)
  beta_threshold_out[[i]]=beta_sim
  ari_threshold[i]<-unlist(
    classAgreement(table(og_group,
                         beta_sim$optimal_model_results$classification))[4])
}

ari_threshold_mean=mean(ari_threshold)
ari_threshold_sd=sd(ari_threshold)

## Simulation study for just KN. model as K.. model is 
## selected as the optimal model each time.
df_threshold_list_kn=list()
beta_threshold_out_kn=list()
ari_threshold_kn=vector()
og_group_kn=vector()
for(i in 1:length(seed_vec))
{
  set.seed(seed_vec[i])
  C=600000
  probs=c(0.35,0.7)
  index=runif(C)
  R=4
  df=matrix(NA,ncol = R,nrow = C)
  og_group_kn=vector(length = C)
  true_dmc_kn=vector(length = C)
  for(c in 1:C)
  {
    if(index[c]<probs[1])
    {
      df[c,]=rbeta(4,2,20) #hypomethylated
      og_group_kn[c]=1
    }else if(index[c]<probs[2])
    {
      df[c,]=rbeta(R,4,3) #hemimethylated
      og_group_kn[c]=2
    }else 
    {
      df[c,]=rbeta(4,20,2) #hypermethylated
      og_group_kn[c]=3
    }
  }
  
  df=cbind(df,og_group_kn)
  df=as.data.frame(df)
  colnames(df)<-c("P_Sam1_1","P_Sam1_2","P_Sam1_3","P_Sam1_4","Original_Group")
  
  df_threshold_list_kn[[i]]=df
  beta_sim_kn<-betaclust(df,3,4,1,c("KN."),parallel_process = TRUE)
  beta_threshold_out_kn[[i]]=beta_sim_kn
  ari_threshold_kn[i]<-unlist(
    classAgreement(table(og_group_kn,
                         beta_sim_kn$optimal_model_results$classification))[4])
}

ari_threshold_mean_kn=mean(ari_threshold_kn)
ari_threshold_sd_kn=sd(ari_threshold_kn)

ari_threshold_k_kn=vector()
for(i in 1:100)
{
  k=beta_threshold_out[[i]]$optimal_model_results$classification
  kn=beta_threshold_out_kn[[i]]$optimal_model_results$classification
  ari_threshold_k_kn[i]<-unlist(classAgreement(table(k,kn))[4])
}
ari_threshold_mean_k_kn=mean(ari_threshold_k_kn)
ari_threshold_sd_k_kn=sd(ari_threshold_k_kn)


## Simulation study for K.R model
### simulation checking

df_list=list()
auc_df=matrix(NA,nrow=100,ncol=9)
wd_df=matrix(NA,nrow=100,ncol=9)
fdr=vector()
sens=vector()
spec=vector()
ppv=vector()
ari=vector()
beta_out=list()
beta_out_dmc=list()

for(i in 1:length(seed_vec))
{
  #i=1
  set.seed(seed_vec[i])
  C=600000
  probs=c(0.1,0.15,0.3,0.4,0.45,0.6,0.75,0.90)
  index=runif(C)
  R=8
  df=matrix(NA,ncol = R,nrow = C)
  og_group=vector(length = C)
  true_dmc=vector(length = C)
  for(c in 1:C)
  {
    if(index[c]<probs[1])
    {
      df[c,1:4]=rbeta(4,20,2) #hyper
      df[c,5:8]=rbeta(4,4,3) #hemimethylated
      og_group[c]=1
      true_dmc[c]=1
    }else if(index[c]<probs[2])
    {
      df[c,]=rbeta(R,20,2) #hypermethylated-hypermethylated
      og_group[c]=2
      true_dmc[c]=0
    }else if(index[c]<probs[3])
    {
      df[c,1:4]=rbeta(4,20,2) #hypermethylated
      df[c,5:8]=rbeta(4,2,20) #hypomethylated
      og_group[c]=3
      true_dmc[c]=1
    }else if(index[c]<probs[4])
    {
      df[c,1:4]=rbeta(4,4,3) #hemimethylated
      df[c,5:8]=rbeta(4,2,20) #hypomethylated
      og_group[c]=4
      true_dmc[c]=1
    }else if(index[c]<probs[5])
    {
      df[c,1:4]=rbeta(4,2,20) #hypomethylated
      df[c,5:8]=rbeta(4,4,3) #hemimethylated
      og_group[c]=5
      true_dmc[c]=1
    }else if(index[c]<probs[6])
    {
      df[c,]=rbeta(R,4,3) #hemimethylated-hemimethylated
      og_group[c]=6
      true_dmc[c]=0
    }else if(index[c]<probs[7])
    {
      df[c,]=rbeta(R,2,20) #hypomethylated - hypomethylated
      og_group[c]=7
      true_dmc[c]=0
    }else if(index[c]<probs[8])
    {
      df[c,1:4]=rbeta(4,2,20) #hypomethylated
      df[c,5:8]=rbeta(4,20,2) #hypermethylated
      og_group[c]=8
      true_dmc[c]=1
    }else
    {
      df[c,1:4]=rbeta(4,4,3) #hemimethylated
      df[c,5:8]=rbeta(4,20,2) #hypermethylated
      og_group[c]=9
      true_dmc[c]=1
    }
  }
  
  df=cbind(df,og_group,true_dmc)
  df=as.data.frame(df)
  colnames(df)<-c("P_Sam1_1","P_Sam1_2","P_Sam1_3","P_Sam1_4","P_Sam2_1",
                  "P_Sam2_2","P_Sam2_3","P_Sam2_4","Original_Group","True_DMC")
  
  df_list[[i]]=df
  beta_sim<-betaclust(df[,1:8],3,4,2,"K.R",parallel_process = TRUE)
  beta_out[[i]]=beta_sim
  
  ## Using AUC and WD metric for finding similarity in the distributions 
  ## estimated in each cluster
  alpha=beta_sim$optimal_model_results$alpha
  delta=beta_sim$optimal_model_results$delta
  n <- 1000
  auc_vec_sim=vector()
  L <- 1001
  d_vec_sim=vector()
  z <- seq(0, 1, length = L)
  for(j in 1:beta_sim$K)
  {
    shape1_1 = alpha[j,1]
    shape1_2 = delta[j,1];
    shape2_1 = alpha[j,2];
    shape2_2 = delta[j,2];
    group_1 <- rbeta(n, shape1 = shape1_1, shape2 = shape1_2)
    group_2 <- rbeta(n, shape1 = shape2_1, shape2 = shape2_2)
    auc_dat <- data.frame(predictor=c(group_1, group_2), 
                          response=factor(c(rep(0,n), rep(1,n))))
    auc_value <- auc(predictor = auc_dat$predictor, 
                     response  = auc_dat$response)
    auc_vec_sim[j]=unlist(auc_value)
    
    
    alpha1 = alpha[j,1]
    delta1 = delta[j,1];
    alpha2 = alpha[j,2];
    delta2 = delta[j,2];
    
    d <- sum(abs(pbeta(z, alpha1, delta1)- pbeta(z, alpha2, delta2))) / L
    d_vec_sim[j]=d
    
  }
  auc_df[i,]=auc_vec_sim
  wd_df[i,]=d_vec_sim
  auc_th=which(auc_vec_sim>0.85)
  beta_sim_true_dmc=ifelse(beta_sim$optimal_model_results$classification %in% 
                             auc_th,1,0)
  beta_out_dmc[[i]]=beta_sim_true_dmc
  ## FDR
  conf_matrix<-table(true_dmc,beta_sim_true_dmc)
  tn <- conf_matrix[1, 1]
  fp <- conf_matrix[1, 2]
  fn <- conf_matrix[2, 1]
  tp <- conf_matrix[2, 2]
  
  # Calculate the false discovery rate
  fdr[i] <- fp / (fp + tp)
  spec[i] <-tn/(fp+tn)
  sens[i] <-tp/(tp+fn)
  ppv[i] <-tp/(fp + tp)
  ari[i]<-unlist(classAgreement(
    table(og_group,beta_sim$optimal_model_results$classification))[4])

  
}
auc_mean=colMeans(auc_df)
auc_sd=sapply(auc_df[1:100,],sd)
auc_df=rbind(auc_df,auc_mean)
auc_df=as.data.frame(auc_df)
wd_mean=colMeans(wd_df)
wd_sd=sapply(wd_df[1:100,],sd)
wd_df=rbind(wd_df,wd_mean)
wd_df=as.data.frame(wd_df)

## Mean of each performance metric
perf_metric_bmm=vector()
perf_metric_bmm[1]=mean(fdr)
perf_metric_bmm[2]=mean(sens)
perf_metric_bmm[3]=mean(spec)
perf_metric_bmm[4]=mean(ppv)
perf_metric_bmm[5]=mean(ari)

## SD of each performance metric
perf_metric_sd_bmm=vector()
perf_metric_sd_bmm[1]=sd(fdr)
perf_metric_sd_bmm[2]=sd(sens)
perf_metric_sd_bmm[3]=sd(spec)
perf_metric_sd_bmm[4]=sd(ppv)
perf_metric_sd_bmm[5]=sd(ari)

