## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(eval = TRUE)
#knitr::opts_chunk$set(dev = 'png')
#knitr::opts_chunk$set(dpi=100)

## ----setup2, echo=FALSE,include=FALSE-----------------------------------------
library(plotly)
library(ggplot2)
library(foreach)
library(doParallel)

 # library(devtools)
 # install_github('koyelucd/betaclust',force = TRUE)
 # library(betaclust)

## ----package, include=TRUE, echo=TRUE, message=FALSE,warning=FALSE------------
library(betaclust)

## ----data,include=TRUE, echo=TRUE---------------------------------------------
data(pca.methylation.data)
head(pca.methylation.data)

## ----thresholds,include=TRUE, echo=TRUE---------------------------------------
M <- 3 ## No. of methylation states in a DNA sample
N <- 4 ## No. of patients
R <- 1 ## No. of DNA samples
my.seed <- 190 ## set seed for reproducibility

threshold_out <- betaclust(pca.methylation.data[,2:5], M, N, R,
                         model_names = c("K..","KN."),
                         model_selection = "BIC", parallel_process = FALSE,
                         seed = my.seed)


## ----output0,include=TRUE, echo=TRUE------------------------------------------
summary(threshold_out)

## ----wrapperoutput3,include=TRUE, echo=TRUE,fig.width = 4, fig.height = 4,dev = 'png'----
plot(threshold_out, what = "information criterion", plot_type = "ggplot")

## ----output1,include=TRUE, echo=TRUE------------------------------------------
threshold_points <- threshold_out$optimal_model_results$thresholds
threshold_points$threholds

## ----output2,include=TRUE, echo=TRUE,fig.width = 5, fig.height = 4,dev = 'png'----
plot(threshold_out, what = "fitted density", threshold = TRUE, data = pca.methylation.data[,2:5], patient_number = 1, plot_type = "ggplot")

## ----output3,include=TRUE, echo=TRUE,fig.width = 5, fig.height = 4,dev = 'png',warning=FALSE----
plot(threshold_out, what = "kernel density", threshold = TRUE, data = pca.methylation.data[,2:5], plot_type = "ggplot")

## ----output4,include=TRUE, echo=TRUE,fig.width = 6, fig.height = 5, dev = 'png'----
plot(threshold_out, what = "uncertainty")

## ----dmc,include=TRUE, echo=TRUE----------------------------------------------
M <- 3  ## No. of methylation states in a DNA sample
N <- 4  ## No. of patients
R <- 2  ## No. of DNA samples
my.seed <- 190 ## set seed for reproducibility

dmc_output <- betaclust(pca.methylation.data[,2:9], M, N, R,
                      model_names = "K.R", parallel_process = FALSE,
                      seed = my.seed)

## ----dmcoutput4b, include=TRUE, echo=TRUE-------------------------------------
print(dmc_output$optimal_model_results$DM$AUC)
print(dmc_output$optimal_model_results$DM$WD)

## ----dmcoutput1,include=TRUE, echo=TRUE---------------------------------------
summary(dmc_output)

## ----dmcoutput2,include=TRUE, echo=TRUE,fig.width = 6.5, fig.height = 5,dev = 'png'----
plot(dmc_output, what = "fitted density", plot_type = "ggplot", data = pca.methylation.data[,2:9], sample_name = c("Benign","Tumour"))

## ----dmcoutput3,include=TRUE, echo=TRUE,fig.width = 6.5, fig.height = 5,dev = 'png'----
plot(dmc_output, what = "kernel density", plot_type = "ggplot", data = pca.methylation.data[,2:9])

## ----dmcoutput4,include=TRUE, echo=TRUE,fig.width = 6, fig.height = 5, dev = 'png'----
plot(dmc_output, what = "uncertainty", plot_type = "ggplot")

## ----dmcoutput6,include=TRUE, echo=TRUE,fig.width = 6, fig.height = 5, dev = 'png'----

## the dataframe containing methylation data for the CpG sites identified as the most differentially methylated ones
dmc_df <- DMC_identification(dmc_output, data = pca.methylation.data[,2:9],
                             pca.methylation.data[,1], threshold = 0.06,
                             metric="WD")

## ----dmcoutput5,include=TRUE, echo=TRUE,fig.width = 6, fig.height = 5,dev = 'png'----

##read the legacy data
data(legacy.data)
head(legacy.data)

## create empty dataframes and matrices to store the DMCs related to the RARB genes
ecdf_rarb <- data.frame()
cpg_rarb <- data.frame(matrix(NA, nrow = 0, ncol = 6))
colnames(cpg_rarb) <- colnames(legacy.data)

## split the UCSC_RefGene_name column to read the gene name related to that CpG site
## select the CpG sites related to the RARB genes
a <- 1
for(i in 1:nrow(legacy.data))
{
  str_vec = unlist(strsplit(as.character(legacy.data[i,"UCSC_RefGene_Name"]),"[;]"))
  if(length(str_vec) != 0)
  {
    if(str_vec[[1]] == "RARB")
    {
      cpg_rarb[a,] <- legacy.data[i,]
      a <- a+1
    }
  }
}

## Read the DMCs related to the RARB genes
ecdf_rarb <- dmc_df[dmc_df$IlmnID %in% cpg_rarb$IlmnID,]

## Plot the ecdf of the selected DMCs
ecdf.betaclust(ecdf_rarb[,2:9], R = 2, sample_name = c("Benign","Tumour"))

