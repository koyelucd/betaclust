#' MethylationEPIC manifest data.
#'
#' The dataset contains the manifest data from the Illumina MethylationEPIC beadchip array.
#'
#' @seealso \code{\link{pca.methylation.data}}
#' @format A data frame with 867,525 rows and 6 columns.
#' \itemize{
#'   \item{IlmnID: the unique identifier from the Illumina CG database, i.e. the probe ID.}
#'   \item{Genome_Build: the genome build referenced by the Infinium MethylationEPIC manifest.}
#'   \item{CHR: the chromosome containing the CpG (Genome_Build = 37).}
#'   \item{MAPINFO: the chromosomal coordinates of the CpG.}
#'   \item{UCSC_RefGene_Name: the target gene name(s), from the UCSC database.
#'   Note: multiple listings of the same gene name indicate splice variants.}
#'   \item{UCSC_CpG_Islands_Name: the chromosomal coordinates of the CpG Island from UCSC.}
#'    }
#' @usage data(legacy.data)
"legacy.data"
