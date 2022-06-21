#' @keywords  summary of betaclust object

print.summary.betaclust <- function(x, digits = getOption("digits"), ...)
{
  txt <- paste(rep("-", min(nchar(x$title), getOption("width"))), collapse = "")
  cat(txt)
  cat("\n")
  cat(x$title)
  cat("\n")
  cat(txt)
  #
  cat("\n")
  if(x$K == 0)
  {
    cat("betaclust model with only a noise component:")
  } else
  {
    cat(paste0("betaclust ", x$modelName, " model with ",
               x$K, ifelse(x$K > 1, " components", " component"),
               ":"))
  }
  cat("\n")



  tab <- data.frame("log-likelihood" = x$loglik,  "Information-criterion"=x$information_criterion,"IC-value" = x$ic_output,
                    "CpG-sites" = x$C, "Patients" = x$N, "Samples" = x$R,
                    row.names = "", check.names = FALSE)
  print(tab)
  #
  cat("\nClustering table:")
  #print(table(x$classification))
  print(x$classification)

  cat("\nProportion of CpG sites in each cluster: \n")
  cat(x$prop_data)
  #
  invisible(x)
}
