#' @title DNA methylation dataset of patients suffering from prostate cancer disease.
#'
#' @description The dataset contains pre-processed beta methylation values from \eqn{R=2} sample, collected from \eqn{N=4} patients suffering from prostate cancer disease.
#'
#' @details The raw methylation array data was first quality controlled and preprocessed using the RnBeads package. The array data was then normalized and and probes located outside of CpG sites and on the sex chromosome were filtered out. The CpG sites with missing values were removed from the resulting dataset.
#' @seealso \code{\link{legacy.data}}
#' @format A data frame with 694820 rows and 9 columns. The data contains no missing values.
#' \describe{
#'   \item{IlmnID}{This contains the unique identifier from the Illumina CG database, i.e. the probe ID.}
#'   \item{Patient_benign_1}{This contains the methylation values from benign prostate tissue collected from patient 1.}
#'   \item{Patient_benign_2}{This contains the methylation values from benign prostate tissue collected from patient 2.}
#'   \item{Patient_benign_3}{This contains the methylation values from benign prostate tissue collected from patient 3.}
#'   \item{Patient_benign_4}{This contains the methylation values from benign prostate tissue collected from patient 4.}
#'   \item{Patient_benign_1}{This contains the methylation values from tumor prostate tissue collected from patient 1.}
#'   \item{Patient_benign_2}{This contains the methylation values from tumor prostate tissue collected from patient 2.}
#'   \item{Patient_benign_3}{This contains the methylation values from tumor prostate tissue collected from patient 3.}
#'   \item{Patient_benign_4}{This contains the methylation values from tumor prostate tissue collected from patient 4.}
#'    }
#' @references {Mueller F, Scherer M, Assenov Y, Lutsik P, Walter J, Lengauer T, Bock C (2019). “RnBeads 2.0: comprehensive analysis of DNA methylation data.” Genome Biology, 20(55). doi: 10.1186/s13059-019-1664-9, https://rnbeads.org.}
#' @references {Assenov Y, Mueller F, Lutsik P, Walter J, Lengauer T, Bock C (2014). “Compehensive Analysis of DNA Methylation Data with RnBeads.” Nature Methods, 11(11), 1138–1140. doi: 10.1038/nmeth.3115, https://rnbeads.org.}
#' @usage data(pca.methylation.data)
"pca.methylation.data"
