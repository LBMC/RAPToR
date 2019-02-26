#' Oudenaarden reference data time series
#'
#' Data from a paper published by Oudenaarden et al.
#' with C elegans gene expression levels (log+1 from rpkm) at 
#' 26 time points.
#'
#' @docType data
#'
#' @usage data(oud_ref)
#'
#' @format a list with \code{X} being the gene expression matrix and \code{time.series} being the timepoints (in hours)
#'
#' @keywords datasets
#'
#' @references Oudenaarden et al. (2013) Nat Genet. 2013 Nov; 45(11): 1337–1344. 
#' (\href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3812263/}{PubMed})
#'
#' @source \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49043}{GEO dataset}
#'
#' @examples
#' data(oud_ref)
#' times <- attr(oud_ref, "time.series")
#' gene1 <- oud_ref$X[,1]
#' \donttest{plot(times, gene1, type='l')}
"oud_ref"