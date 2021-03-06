% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.betaclust.R
\name{plot.betaclust}
\alias{plot.betaclust}
\title{Plots for visualizing the betaclust class object}
\usage{
\method{plot}{betaclust}(
  object,
  what = "density",
  plot_type = "ggplot",
  title = NULL,
  scale_param = "free_y"
)
}
\arguments{
\item{object}{A betaclust object.}

\item{what}{The different plots that can be obtained are either "density","uncertainty" or "InformationCriterion". (default="density").}

\item{plot_type}{The plot type to be displayed are either "ggplot" or "plotly". (default="ggplot").}

\item{title}{The title that the user wants to display on the graph. If no title is to be displayed the default is "NULL" value.}

\item{scale_param}{The axis that needs to be fixed for density estimates plot for visualizing the C.R clustering solution are either "free_y","free_x" or "free". (default = "free_y").}
}
\description{
This function helps visualise the clustering solution by plotting the density estimates, the uncertainty and the information criterion.
}
\details{
The density estimates under the optimal clustering solution by specifying what = "density" in the function. Interactive plots can also be produced using
 plot_type = "plotly". The uncertainty in the clustering solution can be plotted using what="uncertainty".
The information criterion values for all fitted models can be plotted using what = "InformationCriterion".
}
\seealso{
\code{\link{betaclust}}
}
