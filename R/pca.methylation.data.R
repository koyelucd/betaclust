#' @title DNA methylation dataset of patients suffering from prostate cancer disease.
#'
#' @description The dataset contains pre-processed beta methylation values of R=2 samples which are collected from N=4 patients suffering from prostate cancer disease.
#'
#' @details The raw methylation array data was first quality controlled and preprocessed using RnBeads package. The array data was then normalized and and probes located outside of CpG sites and on the sex chromosome were filtered out. The CpG sites with missing values were removed from the resulting dataset.
#' @seealso \code{\link{legacy.data}}
#' @format A data frame with 694820 rows and 9 columns. The data contains no missing values.
#' \describe{
#'   \item{IlmnID}{This contains the Unique identifier from the Illumina CG database. (The probe ID).}
#'   \item{Patient_benign_1}{This contains the methylation values from benign prostate tissue collected from patient 1.}
#'   \item{Patient_benign_2}{This contains the methylation values from benign prostate tissue collected from patient 2.}
#'   \item{Patient_benign_3}{This contains the methylation values from benign prostate tissue collected from patient 3.}
#'   \item{Patient_benign_4}{This contains the methylation values from benign prostate tissue collected from patient 4.}
#'   \item{Patient_benign_1}{This contains the methylation values from tumor prostate tissue collected from patient 1.}
#'   \item{Patient_benign_2}{This contains the methylation values from tumor prostate tissue collected from patient 2.}
#'   \item{Patient_benign_3}{This contains the methylation values from tumor prostate tissue collected from patient 3.}
#'   \item{Patient_benign_4}{This contains the methylation values from tumor prostate tissue collected from patient 4.}
#'    }
#' @usage data(pca.methylation.data)
"pca.methylation.data"
