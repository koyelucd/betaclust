#' MethylationEPIC manifest data.
#'
#' A dataset containing a subset of the manifest data from the Illumina MethylationEPIC beadchip array. A subset of the complete dataset has been uploaded in the package for testing purpose. The complete dataset is available on \href{https://github.com/koyelucd/betaclust}{GitHub}.
#'
#' @seealso \code{\link{pca.methylation.data}}
#' @format A data frame with 5,080 rows and 6 columns.
#' \itemize{
#'   \item{IlmnID: The unique identifier from the Illumina CG database, i.e. the probe ID.}
#'   \item{Genome_Build: The genome build referenced by the Infinium MethylationEPIC manifest.}
#'   \item{CHR: The chromosome containing the CpG (Genome_Build = 37).}
#'   \item{MAPINFO: The chromosomal coordinates of the CpG sites.}
#'   \item{UCSC_RefGene_Name: The target gene name(s), from the UCSC database.
#'   Note: multiple listings of the same gene name indicate splice variants.}
#'   \item{UCSC_CpG_Islands_Name: The chromosomal coordinates of the CpG Island from UCSC.}
#'    }
#' @usage data(legacy.data)
"legacy.data"
