#' Reinke reference data time series
#'
#' Reconstructed data from a dataset published by Reinke et al.
#' Time is in hours, scaled on the Oudenaarden 20C time series.
#' 
#' To obtain this reference, an ICA was performed joining the Oudenaarden and 
#' Reinke datasets to get clear development dynamics, and interpolation done on the 
#' 28-87h time period. This was necessary as the data is quite old and noisy.
#'
#' @docType data
#'
#' @usage data(reinke_ref)
#'
#' @format a list with \code{X} being the gene expression matrix and \code{time.series} being the timepoints (in minutes)
#'
#' @keywords datasets
#'
#' @references Reinke, et al.  (2004) Development 131.2 : 311-323.
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/14668411}{PubMed})
#'
#' @source \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE737}{GEO dataset}
#'
#' @examples
#' data(reinke_ref)
#' times <- reinke_ref$time.series
#' gene1 <- reinke_ref$X[500,]
#' \donttest{plot(times, gene1, type='l')}
"reinke_ref"