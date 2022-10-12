library(devtools)
install_github('koyelucd/betaclust',force = TRUE)
library(betaclust)
install.packages("e1071")
library(e1071)


## Hypomethylated generated from Beta(2,20)
## Hemimethylated generated from Beta(4,3)
## Hypermethylated generated from Beta(20,2)


## 10 Simulated dataset 

## simulation study by changing seed
## 190,200,210,2200,1500,3000,100,9156,8237,172
set.seed(190)
C=20000
probs=c(0.1,0.15,0.3,0.4,0.45,0.6,0.75,0.90)
index=runif(C)
R=8
df=matrix(NA,ncol = R,nrow = C)
sd10_og1=vector(length = C)
sd10_og2=vector(length = C)
for(c in 1:C)
{
  if(index[c]<probs[1])
  {
    df[c,1:4]=rbeta(4,20,2) #hyper
    df[c,5:8]=rbeta(4,4,3) #hemimethylated
    sd10_og1[c]=1
    sd10_og2[c]=1
  }else if(index[c]<probs[2])
  {
    df[c,]=rbeta(R,20,2) #hypermethylated-hypermethylated
    sd10_og1[c]=1
    sd10_og2[c]=2
  }else if(index[c]<probs[3])
  {
    df[c,1:4]=rbeta(4,20,2) #hypermethylated
    df[c,5:8]=rbeta(4,2,20) #hypomethylated
    sd10_og1[c]=1
    sd10_og2[c]=3
  }else if(index[c]<probs[4])
  {
    df[c,1:4]=rbeta(4,4,3) #hemimethylated
    df[c,5:8]=rbeta(4,2,20) #hypomethylated
    sd10_og1[c]=2
    sd10_og2[c]=4
  }else if(index[c]<probs[5])
  {
    df[c,1:4]=rbeta(4,2,20) #hypomethylated
    df[c,5:8]=rbeta(4,4,3) #hemimethylated
    sd10_og1[c]=3
    sd10_og2[c]=5
  }else if(index[c]<probs[6])
  {
    df[c,]=rbeta(R,4,3) #hemimethylated-hemimethylated
    sd10_og1[c]=2
    sd10_og2[c]=6
  }else if(index[c]<probs[7])
  {
    df[c,]=rbeta(R,2,20) #hypomethylated - hypomethylated
    sd10_og1[c]=3
    sd10_og2[c]=7
  }else if(index[c]<probs[8])
  {
    df[c,1:4]=rbeta(4,2,20) #hypomethylated
    df[c,5:8]=rbeta(4,20,2) #hypermethylated
    sd10_og1[c]=3
    sd10_og2[c]=8
  }else
  {
    df[c,1:4]=rbeta(4,4,3) #hemimethylated
    df[c,5:8]=rbeta(4,20,2) #hypermethylated
    sd10_og1[c]=2
    sd10_og2[c]=9
  }
}


df=as.data.frame(df)
colnames(df)<-c("P_Sam1_1","P_Sam1_2","P_Sam1_3","P_Sam1_4","P_Sam2_1","P_Sam2_2","P_Sam2_3","P_Sam2_4")

df1=df


M=3
patients=4
samples=1
my.seed=190
out_c1<-betaclust(df1[,1:4],M,patients,samples,model_names = c("K.."),model_selection = "BIC",seed=my.seed)
out_c2<-betaclust(df2[,1:4],M,patients,samples,model_names = c("K.."),model_selection = "BIC",seed=my.seed)
out_c3<-betaclust(df3[,1:4],M,patients,samples,model_names = c("K.."),model_selection = "BIC",seed=my.seed)
out_c4<-betaclust(df4[,1:4],M,patients,samples,model_names = c("K.."),model_selection = "BIC",seed=my.seed)
out_c5<-betaclust(df5[,1:4],M,patients,samples,model_names = c("K.."),model_selection = "BIC",seed=my.seed)
out_c6<-betaclust(df6[,1:4],M,patients,samples,model_names = c("K.."),model_selection = "BIC",seed=my.seed)
out_c7<-betaclust(df7[,1:4],M,patients,samples,model_names = c("K.."),model_selection = "BIC",seed=my.seed)
out_c8<-betaclust(df8[,1:4],M,patients,samples,model_names = c("K.."),model_selection = "BIC",seed=my.seed)
out_c9<-betaclust(df9[,1:4],M,patients,samples,model_names = c("K.."),model_selection = "BIC",seed=my.seed)
out_c10<-betaclust(df10[,1:4],M,patients,samples,model_names = c("K.."),model_selection = "BIC",seed=my.seed)

M=3
patients=4
samples=1
my.seed=190
out_cn1<-betaclust(df1[,1:4],M,patients,samples,model_names = c("KN."),model_selection = "BIC",seed=my.seed)
out_cn2<-betaclust(df2[,1:4],M,patients,samples,model_names = c("KN."),model_selection = "BIC",seed=my.seed)
out_cn3<-betaclust(df3[,1:4],M,patients,samples,model_names = c("KN."),model_selection = "BIC",seed=my.seed)
out_cn4<-betaclust(df4[,1:4],M,patients,samples,model_names = c("KN."),model_selection = "BIC",seed=my.seed)
out_cn5<-betaclust(df5[,1:4],M,patients,samples,model_names = c("KN."),model_selection = "BIC",seed=my.seed)
out_cn6<-betaclust(df6[,1:4],M,patients,samples,model_names = c("KN."),model_selection = "BIC",seed=my.seed)
out_cn7<-betaclust(df7[,1:4],M,patients,samples,model_names = c("KN."),model_selection = "BIC",seed=my.seed)
out_cn8<-betaclust(df8[,1:4],M,patients,samples,model_names = c("KN."),model_selection = "BIC",seed=my.seed)
out_cn9<-betaclust(df9[,1:4],M,patients,samples,model_names = c("KN."),model_selection = "BIC",seed=my.seed)
out_cn10<-betaclust(df10[,1:4],M,patients,samples,model_names = c("KN."),model_selection = "BIC",seed=my.seed)


M=3
patients=4
samples=2
my.seed=190
out_cr1<-betaclust(df1[,1:8],M,patients,samples,model_names = c("K.R"),model_selection = "BIC",seed=my.seed)
out_cr2<-betaclust(df2[,1:8],M,patients,samples,model_names = c("K.R"),model_selection = "BIC",seed=my.seed)
out_cr3<-betaclust(df3[,1:8],M,patients,samples,model_names = c("K.R"),model_selection = "BIC",seed=my.seed)
out_cr4<-betaclust(df4[,1:8],M,patients,samples,model_names = c("K.R"),model_selection = "BIC",seed=my.seed)
out_cr5<-betaclust(df5[,1:8],M,patients,samples,model_names = c("K.R"),model_selection = "BIC",seed=my.seed)
out_cr6<-betaclust(df6[,1:8],M,patients,samples,model_names = c("K.R"),model_selection = "BIC",seed=my.seed)
out_cr7<-betaclust(df7[,1:8],M,patients,samples,model_names = c("K.R"),model_selection = "BIC",seed=my.seed)
out_cr8<-betaclust(df8[,1:8],M,patients,samples,model_names = c("K.R"),model_selection = "BIC",seed=my.seed)
out_cr9<-betaclust(df9[,1:8],M,patients,samples,model_names = c("K.R"),model_selection = "BIC",seed=my.seed)
out_cr10<-betaclust(df10[,1:8],M,patients,samples,model_names = c("K.R"),model_selection = "BIC",seed=my.seed)

## Calculating adjusted Rand index for comparing
# K.. and KN. models
ari_k_kn<-vector()

# K.. model and true membership
ari_k_og<-vector()

# KN. model and true membership
ari_kn_og<-vector()

# K.R and true membership
ari_kr_og<-vector()


tab1<-table(out_c1$optimal_model_results$data[,5],sd1_og1)
ari_k_og[1]=classAgreement(tab1)[[4]]
tab2<-table(out_c2$optimal_model_results$data[,5],sd2_og1)
ari_k_og[2]=classAgreement(tab2)[[4]]
tab3<-table(out_c3$optimal_model_results$data[,5],sd3_og1)
ari_k_og[3]=classAgreement(tab3)[[4]]
tab4<-table(out_c4$optimal_model_results$data[,5],sd4_og1)
ari_k_og[4]=classAgreement(tab4)[[4]]
tab5<-table(out_c5$optimal_model_results$data[,5],sd5_og1)
ari_k_og[5]=classAgreement(tab5)[[4]]
tab6<-table(out_c6$optimal_model_results$data[,5],sd6_og1)
ari_k_og[6]=classAgreement(tab6)[[4]]
tab7<-table(out_c7$optimal_model_results$data[,5],sd7_og1)
ari_k_og[7]=classAgreement(tab7)[[4]]
tab8<-table(out_c8$optimal_model_results$data[,5],sd8_og1)
ari_k_og[8]=classAgreement(tab8)[[4]]
tab9<-table(out_c9$optimal_model_results$data[,5],sd9_og1)
ari_k_og[9]=classAgreement(tab9)[[4]]
tab10<-table(out_c10$optimal_model_results$data[,5],sd10_og1)
ari_k_og[10]=classAgreement(tab10)[[4]]

tab1<-table(out_cr1$optimal_model_results$data[,9],sd1_og2)
ari_kr_og[1]=classAgreement(tab1)[[4]]
tab2<-table(out_cr2$optimal_model_results$data[,9],sd2_og2)
ari_kr_og[2]=classAgreement(tab2)[[4]]
tab3<-table(out_cr3$optimal_model_results$data[,9],sd3_og2)
ari_kr_og[3]=classAgreement(tab3)[[4]]
tab4<-table(out_cr4$optimal_model_results$data[,9],sd4_og2)
ari_kr_og[4]=classAgreement(tab4)[[4]]
tab5<-table(out_cr5$optimal_model_results$data[,9],sd5_og2)
ari_kr_og[5]=classAgreement(tab5)[[4]]
tab6<-table(out_cr6$optimal_model_results$data[,9],sd6_og2)
ari_kr_og[6]=classAgreement(tab6)[[4]]
tab7<-table(out_cr7$optimal_model_results$data[,9],sd7_og2)
ari_kr_og[7]=classAgreement(tab7)[[4]]
tab8<-table(out_cr8$optimal_model_results$data[,9],sd8_og2)
ari_kr_og[8]=classAgreement(tab8)[[4]]
tab9<-table(out_cr9$optimal_model_results$data[,9],sd9_og2)
ari_kr_og[9]=classAgreement(tab9)[[4]]
tab10<-table(out_cr10$optimal_model_results$data[,9],sd10_og2)
ari_kr_og[10]=classAgreement(tab10)[[4]]

tab1<-table(out_cn1$optimal_model_results$data[,5],sd1_og1)
ari_kn_og[1]=classAgreement(tab1)[[4]]
tab2<-table(out_cn2$optimal_model_results$data[,5],sd2_og1)
ari_kn_og[2]=classAgreement(tab2)[[4]]
tab3<-table(out_cn3$optimal_model_results$data[,5],sd3_og1)
ari_kn_og[3]=classAgreement(tab3)[[4]]
tab4<-table(out_cn4$optimal_model_results$data[,5],sd4_og1)
ari_kn_og[4]=classAgreement(tab4)[[4]]
tab5<-table(out_cn5$optimal_model_results$data[,5],sd5_og1)
ari_kn_og[5]=classAgreement(tab5)[[4]]
tab6<-table(out_cn6$optimal_model_results$data[,5],sd6_og1)
ari_kn_og[6]=classAgreement(tab6)[[4]]
tab7<-table(out_cn7$optimal_model_results$data[,5],sd7_og1)
ari_kn_og[7]=classAgreement(tab7)[[4]]
tab8<-table(out_cn8$optimal_model_results$data[,5],sd8_og1)
ari_kn_og[8]=classAgreement(tab8)[[4]]
tab9<-table(out_cn9$optimal_model_results$data[,5],sd9_og1)
ari_kn_og[9]=classAgreement(tab9)[[4]]
tab10<-table(out_cn10$optimal_model_results$data[,5],sd10_og1)
ari_kn_og[10]=classAgreement(tab10)[[4]]

tab1<-table(out_c1$optimal_model_results$data[,5],out_cn1$optimal_model_results$data[,5])
ari_k_kn[1]=classAgreement(tab1)[[4]]
tab2<-table(out_c2$optimal_model_results$data[,5],out_cn2$optimal_model_results$data[,5])
ari_k_kn[2]=classAgreement(tab2)[[4]]
tab3<-table(out_c3$optimal_model_results$data[,5],out_cn3$optimal_model_results$data[,5])
ari_k_kn[3]=classAgreement(tab3)[[4]]
tab4<-table(out_c4$optimal_model_results$data[,5],out_cn4$optimal_model_results$data[,5])
ari_k_kn[4]=classAgreement(tab4)[[4]]
tab5<-table(out_c5$optimal_model_results$data[,5],out_cn5$optimal_model_results$data[,5])
ari_k_kn[5]=classAgreement(tab5)[[4]]
tab6<-table(out_c6$optimal_model_results$data[,5],out_cn6$optimal_model_results$data[,5])
ari_k_kn[6]=classAgreement(tab6)[[4]]
tab7<-table(out_c7$optimal_model_results$data[,5],out_cn7$optimal_model_results$data[,5])
ari_k_kn[7]=classAgreement(tab7)[[4]]
tab8<-table(out_c8$optimal_model_results$data[,5],out_cn8$optimal_model_results$data[,5])
ari_k_kn[8]=classAgreement(tab8)[[4]]
tab9<-table(out_c9$optimal_model_results$data[,5],out_cn9$optimal_model_results$data[,5])
ari_k_kn[9]=classAgreement(tab9)[[4]]
tab10<-table(out_c10$optimal_model_results$data[,5],out_cn10$optimal_model_results$data[,5])
ari_k_kn[10]=classAgreement(tab10)[[4]]



## 10 Simulated dataset

## simulation study by changing seed
# ##  140, 3987,6543,8632,3456,4567,5678,1076,5555,4444
set.seed(4444)
C=20000
probs=c(0.1,0.15,0.3,0.4,0.45,0.6,0.75,0.90)
index=runif(C)
R=8
df=matrix(NA,ncol = R,nrow = C)
set2_sd10_og1=vector(length = C)
set2_sd10_og2=vector(length = C)
for(c in 1:C)
{
  if(index[c]<probs[1])
  {
    df[c,1:4]=rbeta(4,20,2) #hypermethylated
    df[c,5:8]=rbeta(4,4,3) #hemimethylated
    set2_sd10_og1[c]=1
    set2_sd10_og2[c]=1
  }else if(index[c]<probs[2])
  {
    df[c,]=rbeta(R,20,2) #hypermethylated-hypermethylated
    set2_sd10_og1[c]=1
    set2_sd10_og2[c]=2
  }else if(index[c]<probs[3])
  {
    df[c,1:4]=rbeta(4,20,2) #hypermethylated
    df[c,5:8]=rbeta(4,2,20) #hypomethylated
    set2_sd10_og1[c]=1
    set2_sd10_og2[c]=3
  }else if(index[c]<probs[4])
  {
    df[c,1:4]=rbeta(4,4,3) #hemimethylated
    df[c,5:8]=rbeta(4,2,20) #hypomethylated
    set2_sd10_og1[c]=2
    set2_sd10_og2[c]=4
  }else if(index[c]<probs[5])
  {
    df[c,1:4]=rbeta(4,2,20) #hypomethylated
    df[c,5:8]=rbeta(4,4,3) #hemimethylated
    set2_sd10_og1[c]=3
    set2_sd10_og2[c]=5
  }else if(index[c]<probs[6])
  {
    df[c,]=rbeta(R,4,3) #hemimethylated-hemimethylated
    set2_sd10_og1[c]=2
    set2_sd10_og2[c]=6
  }else if(index[c]<probs[7])
  {
    df[c,]=rbeta(R,2,20) #hypomethylated - hypomethylated
    set2_sd10_og1[c]=3
    set2_sd10_og2[c]=7
  }else if(index[c]<probs[8])
  {
    df[c,1:4]=rbeta(4,2,20) #hypomethylated
    df[c,5:8]=rbeta(4,20,2) #hypermethylated
    set2_sd10_og1[c]=3
    set2_sd10_og2[c]=8
  }else
  {
    df[c,1:4]=rbeta(4,4,3) #hemimethylated
    df[c,5:8]=rbeta(4,20,2) #hypermethylated
    set2_sd10_og1[c]=2
    set2_sd10_og2[c]=9
  }
}


df=as.data.frame(df)
colnames(df)<-c("P_Sam1_1","P_Sam1_2","P_Sam1_3","P_Sam1_4","P_Sam2_1","P_Sam2_2","P_Sam2_3","P_Sam2_4")

df10_set2=df

M=3
patients=4
samples=1
my.seed=190
out1_set2<-betaclust(df1_set2[,1:4],M,patients,samples,model_names = c("K..","KN."),model_selection = "BIC",seed=my.seed)
out2_set2<-betaclust(df2_set2[,1:4],M,patients,samples,model_names = c("K..","KN."),model_selection = "BIC",seed=my.seed)
out3_set2<-betaclust(df3_set2[,1:4],M,patients,samples,model_names = c("K..","KN."),model_selection = "BIC",seed=my.seed)
out4_set2<-betaclust(df4_set2[,1:4],M,patients,samples,model_names = c("K..","KN."),model_selection = "BIC",seed=my.seed)
out5_set2<-betaclust(df5_set2[,1:4],M,patients,samples,model_names = c("K..","KN."),model_selection = "BIC",seed=my.seed)
out6_set2<-betaclust(df6_set2[,1:4],M,patients,samples,model_names = c("K..","KN."),model_selection = "BIC",seed=my.seed)
out7_set2<-betaclust(df7_set2[,1:4],M,patients,samples,model_names = c("K..","KN."),model_selection = "BIC",seed=my.seed)
out8_set2<-betaclust(df8_set2[,1:4],M,patients,samples,model_names = c("K..","KN."),model_selection = "BIC",seed=my.seed)
out9_set2<-betaclust(df9_set2[,1:4],M,patients,samples,model_names = c("K..","KN."),model_selection = "BIC",seed=my.seed)
out10_set2<-betaclust(df10_set2[,1:4],M,patients,samples,model_names = c("K..","KN."),model_selection = "BIC",seed=my.seed)


M=3
patients=4
samples=1
my.seed=190
out_set2_c1<-betaclust(df1_set2[,1:4],M,patients,samples,model_names = c("K.."),model_selection = "BIC",seed=my.seed)
out_set2_c2<-betaclust(df2_set2[,1:4],M,patients,samples,model_names = c("K.."),model_selection = "BIC",seed=my.seed)
out_set2_c3<-betaclust(df3_set2[,1:4],M,patients,samples,model_names = c("K.."),model_selection = "BIC",seed=my.seed)
out_set2_c4<-betaclust(df4_set2[,1:4],M,patients,samples,model_names = c("K.."),model_selection = "BIC",seed=my.seed)
out_set2_c5<-betaclust(df5_set2[,1:4],M,patients,samples,model_names = c("K.."),model_selection = "BIC",seed=my.seed)
out_set2_c6<-betaclust(df6_set2[,1:4],M,patients,samples,model_names = c("K.."),model_selection = "BIC",seed=my.seed)
out_set2_c7<-betaclust(df7_set2[,1:4],M,patients,samples,model_names = c("K.."),model_selection = "BIC",seed=my.seed)
out_set2_c8<-betaclust(df8_set2[,1:4],M,patients,samples,model_names = c("K.."),model_selection = "BIC",seed=my.seed)
out_set2_c9<-betaclust(df9_set2[,1:4],M,patients,samples,model_names = c("K.."),model_selection = "BIC",seed=my.seed)
out_set2_c10<-betaclust(df10_set2[,1:4],M,patients,samples,model_names = c("K.."),model_selection = "BIC",seed=my.seed)

M=3
patients=4
samples=1
my.seed=190
out_set2_cn1<-betaclust(df1_set2[,1:4],M,patients,samples,model_names = c("KN."),model_selection = "BIC",seed=my.seed)
out_set2_cn2<-betaclust(df2_set2[,1:4],M,patients,samples,model_names = c("KN."),model_selection = "BIC",seed=my.seed)
out_set2_cn3<-betaclust(df3_set2[,1:4],M,patients,samples,model_names = c("KN."),model_selection = "BIC",seed=my.seed)
out_set2_cn4<-betaclust(df4_set2[,1:4],M,patients,samples,model_names = c("KN."),model_selection = "BIC",seed=my.seed)
out_set2_cn5<-betaclust(df5_set2[,1:4],M,patients,samples,model_names = c("KN."),model_selection = "BIC",seed=my.seed)
out_set2_cn6<-betaclust(df6_set2[,1:4],M,patients,samples,model_names = c("KN."),model_selection = "BIC",seed=my.seed)
out_set2_cn7<-betaclust(df7_set2[,1:4],M,patients,samples,model_names = c("KN."),model_selection = "BIC",seed=my.seed)
out_set2_cn8<-betaclust(df8_set2[,1:4],M,patients,samples,model_names = c("KN."),model_selection = "BIC",seed=my.seed)
out_set2_cn9<-betaclust(df9_set2[,1:4],M,patients,samples,model_names = c("KN."),model_selection = "BIC",seed=my.seed)
out_set2_cn10<-betaclust(df10_set2[,1:4],M,patients,samples,model_names = c("KN."),model_selection = "BIC",seed=my.seed)


M=3
patients=4
samples=2
my.seed=190
out_set2_cr1<-betaclust(df1_set2[,1:8],M,patients,samples,model_names = c("K.R"),model_selection = "BIC",seed=my.seed)
out_set2_cr2<-betaclust(df2_set2[,1:8],M,patients,samples,model_names = c("K.R"),model_selection = "BIC",seed=my.seed)
out_set2_cr3<-betaclust(df3_set2[,1:8],M,patients,samples,model_names = c("K.R"),model_selection = "BIC",seed=my.seed)
out_set2_cr4<-betaclust(df4_set2[,1:8],M,patients,samples,model_names = c("K.R"),model_selection = "BIC",seed=my.seed)
out_set2_cr5<-betaclust(df5_set2[,1:8],M,patients,samples,model_names = c("K.R"),model_selection = "BIC",seed=my.seed)
out_set2_cr6<-betaclust(df6_set2[,1:8],M,patients,samples,model_names = c("K.R"),model_selection = "BIC",seed=my.seed)
out_set2_cr7<-betaclust(df7_set2[,1:8],M,patients,samples,model_names = c("K.R"),model_selection = "BIC",seed=my.seed)
out_set2_cr8<-betaclust(df8_set2[,1:8],M,patients,samples,model_names = c("K.R"),model_selection = "BIC",seed=my.seed)
out_set2_cr9<-betaclust(df9_set2[,1:8],M,patients,samples,model_names = c("K.R"),model_selection = "BIC",seed=my.seed)
out_set2_cr10<-betaclust(df10_set2[,1:8],M,patients,samples,model_names = c("K.R"),model_selection = "BIC",seed=my.seed)

## Calculating adjusted Rand index for comparing
# K.. and KN. models
ari_set2_k_kn<-vector()

# K.. model and true membership
ari_set2_k_og<-vector()

# KN. model and true membership
ari_set2_kn_og<-vector()

# K.R and true membership
ari_set2_kr_og<-vector()

tab1<-table(out1_set2$optimal_model_results$data[,5],set2_sd1_og1)
ari_set2_k_og[1]=classAgreement(tab1)[[4]]
tab2<-table(out2_set2$optimal_model_results$data[,5],set2_sd2_og1)
ari_set2_k_og[2]=classAgreement(tab2)[[4]]
tab3<-table(out3_set2$optimal_model_results$data[,5],set2_sd3_og1)
ari_set2_k_og[3]=classAgreement(tab3)[[4]]
tab4<-table(out4_set2$optimal_model_results$data[,5],set2_sd4_og1)
ari_set2_k_og[4]=classAgreement(tab4)[[4]]
tab5<-table(out5_set2$optimal_model_results$data[,5],set2_sd5_og1)
ari_set2_k_og[5]=classAgreement(tab5)[[4]]
tab6<-table(out6_set2$optimal_model_results$data[,5],set2_sd6_og1)
ari_set2_k_og[6]=classAgreement(tab6)[[4]]
tab7<-table(out7_set2$optimal_model_results$data[,5],set2_sd7_og1)
ari_set2_k_og[7]=classAgreement(tab7)[[4]]
tab8<-table(out8_set2$optimal_model_results$data[,5],set2_sd8_og1)
ari_set2_k_og[8]=classAgreement(tab8)[[4]]
tab9<-table(out9_set2$optimal_model_results$data[,5],set2_sd9_og1)
ari_set2_k_og[9]=classAgreement(tab9)[[4]]
tab10<-table(out10_set2$optimal_model_results$data[,5],set2_sd10_og1)
ari_set2_k_og[10]=classAgreement(tab10)[[4]]

tab1<-table(out_set2_cr1$optimal_model_results$data[,9],set2_sd1_og2)
ari_set2_kr_og[1]=classAgreement(tab1)[[4]]
tab2<-table(out_set2_cr2$optimal_model_results$data[,9],set2_sd2_og2)
ari_set2_kr_og[2]=classAgreement(tab2)[[4]]
tab3<-table(out_set2_cr3$optimal_model_results$data[,9],set2_sd3_og2)
ari_set2_kr_og[3]=classAgreement(tab3)[[4]]
tab4<-table(out_set2_cr4$optimal_model_results$data[,9],set2_sd4_og2)
ari_set2_kr_og[4]=classAgreement(tab4)[[4]]
tab5<-table(out_set2_cr5$optimal_model_results$data[,9],set2_sd5_og2)
ari_set2_kr_og[5]=classAgreement(tab5)[[4]]
tab6<-table(out_set2_cr6$optimal_model_results$data[,9],set2_sd6_og2)
ari_set2_kr_og[6]=classAgreement(tab6)[[4]]
tab7<-table(out_set2_cr7$optimal_model_results$data[,9],set2_sd7_og2)
ari_set2_kr_og[7]=classAgreement(tab7)[[4]]
tab8<-table(out_set2_cr8$optimal_model_results$data[,9],set2_sd8_og2)
ari_set2_kr_og[8]=classAgreement(tab8)[[4]]
tab9<-table(out_set2_cr9$optimal_model_results$data[,9],set2_sd9_og2)
ari_set2_kr_og[9]=classAgreement(tab9)[[4]]
tab10<-table(out_set2_cr10$optimal_model_results$data[,9],set2_sd10_og2)
ari_set2_kr_og[10]=classAgreement(tab10)[[4]]

tab1<-table(out_set2_cn1$optimal_model_results$data[,5],set2_sd1_og1)
ari_set2_kn_og[1]=classAgreement(tab1)[[4]]
tab2<-table(out_set2_cn2$optimal_model_results$data[,5],set2_sd2_og1)
ari_set2_kn_og[2]=classAgreement(tab2)[[4]]
tab3<-table(out_set2_cn3$optimal_model_results$data[,5],set2_sd3_og1)
ari_set2_kn_og[3]=classAgreement(tab3)[[4]]
tab4<-table(out_set2_cn4$optimal_model_results$data[,5],set2_sd4_og1)
ari_set2_kn_og[4]=classAgreement(tab4)[[4]]
tab5<-table(out_set2_cn5$optimal_model_results$data[,5],set2_sd5_og1)
ari_set2_kn_og[5]=classAgreement(tab5)[[4]]
tab6<-table(out_set2_cn6$optimal_model_results$data[,5],set2_sd6_og1)
ari_set2_kn_og[6]=classAgreement(tab6)[[4]]
tab7<-table(out_set2_cn7$optimal_model_results$data[,5],set2_sd7_og1)
ari_set2_kn_og[7]=classAgreement(tab7)[[4]]
tab8<-table(out_set2_cn8$optimal_model_results$data[,5],set2_sd8_og1)
ari_set2_kn_og[8]=classAgreement(tab8)[[4]]
tab9<-table(out_set2_cn9$optimal_model_results$data[,5],set2_sd9_og1)
ari_set2_kn_og[9]=classAgreement(tab9)[[4]]
tab10<-table(out_set2_cn10$optimal_model_results$data[,5],set2_sd10_og1)
ari_set2_kn_og[10]=classAgreement(tab10)[[4]]

tab1<-table(out_set2_c1$optimal_model_results$data[,5],out_set2_cn1$optimal_model_results$data[,5])
ari_set2_k_kn[1]=classAgreement(tab1)[[4]]
tab2<-table(out_set2_c2$optimal_model_results$data[,5],out_set2_cn2$optimal_model_results$data[,5])
ari_set2_k_kn[2]=classAgreement(tab2)[[4]]
tab3<-table(out_set2_c3$optimal_model_results$data[,5],out_set2_cn3$optimal_model_results$data[,5])
ari_set2_k_kn[3]=classAgreement(tab3)[[4]]
tab4<-table(out_set2_c4$optimal_model_results$data[,5],out_set2_cn4$optimal_model_results$data[,5])
ari_set2_k_kn[4]=classAgreement(tab4)[[4]]
tab5<-table(out_set2_c5$optimal_model_results$data[,5],out_set2_cn5$optimal_model_results$data[,5])
ari_set2_k_kn[5]=classAgreement(tab5)[[4]]
tab6<-table(out_set2_c6$optimal_model_results$data[,5],out_set2_cn6$optimal_model_results$data[,5])
ari_set2_k_kn[6]=classAgreement(tab6)[[4]]
tab7<-table(out_set2_c7$optimal_model_results$data[,5],out_set2_cn7$optimal_model_results$data[,5])
ari_set2_k_kn[7]=classAgreement(tab7)[[4]]
tab8<-table(out_set2_c8$optimal_model_results$data[,5],out_set2_cn8$optimal_model_results$data[,5])
ari_set2_k_kn[8]=classAgreement(tab8)[[4]]
tab9<-table(out_set2_c9$optimal_model_results$data[,5],out_set2_cn9$optimal_model_results$data[,5])
ari_set2_k_kn[9]=classAgreement(tab9)[[4]]
tab10<-table(out_set2_c10$optimal_model_results$data[,5],out_set2_cn10$optimal_model_results$data[,5])
ari_set2_k_kn[10]=classAgreement(tab10)[[4]]



# Mean ARI and its standard deviation for comparing K.. and KN. models
m1=c(ari_k_kn,ari_set2_k_kn) 
mean(m1)
sd(m1)

# Mean ARI and its standard deviation for comparing K.. model solution and true membership
m2=c(ari_k_og,ari_set2_k_og)
mean(m2)
sd(m2)

# Mean ARI and its standard deviation for comparing KN. model solution and true membership
m3=c(ari_kn_og,ari_set2_kn_og)
mean(m3)
sd(m3)

# Mean ARI and its standard deviation for comparing K.R model solution and true membership
m4=c(ari_kr_og,ari_set2_kr_og)
mean(m4)
sd(m4)


