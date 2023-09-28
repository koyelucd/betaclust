#' @title DNA methylation data from patients with prostate cancer
#'
#' @description A dataset containing pre-processed beta methylation values from \eqn{R=2} sample types, collected from \eqn{N=4} patients with prostate cancer.
#'
#' @details The raw methylation array data was first quality controlled and preprocessed using the \href{https://rnbeads.org/}{RnBeads} package. The array data were then normalized and and probes located outside of CpG sites and on the sex chromosome were filtered out. The CpG sites with missing values were removed from the resulting dataset. A subset of the complete dataset has been uploaded in the package for testing purposes. The complete dataset is available on \href{https://github.com/koyelucd/betaclust}{GitHub}.
#' @seealso \code{\link{legacy.data}}
#' @format A data frame with 5,067 rows and 9 columns. The data contain no missing values.
#' \itemize{
#'   \item{IlmnID: The unique identifier from the Illumina CG database, i.e. the probe ID.}
#'   \item{Benign_Patient_1: Methylation values from benign prostate tissue from patient 1.}
#'   \item{Benign_Patient_2: Methylation values from benign prostate tissue from patient 2.}
#'   \item{Benign_Patient_3: Methylation values from benign prostate tissue from patient 3.}
#'   \item{Benign_Patient_4: Methylation values from benign prostate tissue from patient 4.}
#'   \item{Tumour_Patient_1: Methylation values from tumor prostate tissue from patient 1.}
#'   \item{Tumour_Patient_2: Methylation values from tumor prostate tissue from patient 2.}
#'   \item{Tumour_Patient_3: Methylation values from tumor prostate tissue from patient 3.}
#'   \item{Tumour_Patient_4: Methylation values from tumor prostate tissue from patient 4.}
#'    }
#' @references {Mueller F, Scherer M, Assenov Y, Lutsik P, Walter J, Lengauer T, Bock C (2019). “RnBeads 2.0: comprehensive analysis of DNA methylation data.” Genome Biology, 20(55).}
#' @references {Assenov Y, Mueller F, Lutsik P, Walter J, Lengauer T, Bock C (2014). “Comprehensive Analysis of DNA Methylation Data with RnBeads.” Nature Methods, 11(11), 1138–1140.}
#' @usage data(pca.methylation.data)
"pca.methylation.data"
