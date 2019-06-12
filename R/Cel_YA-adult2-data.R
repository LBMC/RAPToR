#' C. elegans Young Adult/Adult reference data time series (2)
#'
#' Dataset published by Reinke et al. in 2004
#' Time is in hours, scaled on the Cel_larval 20C time series.
#' Note that this reference series is quite old and very noisy.
#' 
#'
#' @docType data
#'
#' @usage data(Cel_YA_adult2)
#'
#' @format a list with \code{X} being the gene expression matrix and \code{time.series} being the timepoints (in hours post-hatching)
#'
#' @keywords datasets
#'
#' @references Reinke, et al.  (2004) Development 131.2 : 311-323.
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/14668411}{PubMed})
#'
#' @source \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE737}{GEO dataset}
#'
#' @examples
#' data(Cel_YA_adult2)
#' times <- Cel_YA_adult2$time.series
#' gene1 <- Cel_YA_adult2$X[500,]
#' \donttest{plot(times, gene1, type='l')}
"Cel_YA_adult2"