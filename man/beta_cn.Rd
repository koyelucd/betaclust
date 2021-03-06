% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/beta_cn.R
\name{beta_cn}
\alias{beta_cn}
\title{The CN. model}
\usage{
beta_cn(data, M = 3, seed, register = NULL)
}
\arguments{
\item{data}{Methylation values for \eqn{C} CpG sites from \eqn{R=1} samples collected from \eqn{N} patients.}

\item{M}{Number of methylation profiles to be identified in a DNA sample.}

\item{seed}{Seed to allow for reproducibility.}

\item{register}{Setting for registering the parallel backend with the 'foreach' package. To start parallel execution of R code on machine with multiple cores, 'NULL' value needs to be assigned to this parameter.}
}
\value{
A list containing:
\itemize{
   \item cluster_size - the total number of CpG sites identified in each of the K clusters.
   \item llk - a vector containing the log-likelihood value at each step of the EM algorithm.
   \item data - this contains the methylation dataset along with the cluster label for each CpG site.
   \item alpha - this contains the shape parameter 1 for the beta mixture model.
   \item delta - this contains the shape parameter 2 for the mixture model.
   \item tau - the proportion of CpG sites in each cluster.
   \item z - a matrix containing the probability for each CpG site of belonging to each of the \eqn{K} clusters.
   \item uncertainty - the uncertainty of each CpG site's clustering.    }
}
\description{
Fit the CN. model from the family of beta mixture models for DNA methylation data.
             The CN. model analyses a single DNA sample and identifies the thresholds for the different methylation profiles.
}
\details{
This model clusters each of the \eqn{C} CpG sites into one of \eqn{K} methylation profiles, based on data from \eqn{N} patients for one DNA sample (i.e. \eqn{R=1}).
As each CpG site can belong to either of the \eqn{M=3} methylation profiles (hypomethylated, hemimethylated or hypermethylated), the default value of \eqn{K=M=3}.
The CN. model differs from the C.. model as it is less parsimonious, allowing cluster and patient-specific shape parameters.
}
\examples{
\dontrun{
data(pca.methylation.data)
my.seed=190
M=3
data_output=beta_cn(pca.methylation.data[,2:5],M,seed=my.seed)
}
}
\references{
{Microsoft, Weston, S. (2022): foreach: Provides Foreach Looping Construct. R package version 1.5.2. https://CRAN.R-project.org/package=foreach.}
}
\seealso{
\code{\link{beta_c}}

\code{\link{betaclust}}
}
