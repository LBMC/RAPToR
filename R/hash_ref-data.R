#' Hashimshony reference data time series
#'
#' Data from a paper published by Hashimshony et al.
#' with C elegans gene expression levels (log+1 from rpkm) at 
#' 56 time points. The time is in minutes after the 4-cell stage. 
#' Two timepoints thus have negative times (the 1 and 2-cell stage embryos).
#' This dataset is a reference time series for C. elegans embryo development at 20C 
#' from the 1-cell stage to hatching (830min).
#'
#' @docType data
#'
#' @usage data(hash_ref)
#'
#' @format a list with \code{X} being the gene expression matrix and \code{time.series} being the timepoints (in minutes)
#'
#' @keywords datasets
#'
#' @references Hashimshony, et al.  (2015) Nature 519.7542 : 219.
#' (\href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4359913/}{PubMed})
#'
#' @source \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50548}{GEO dataset}
#'
#' @examples
#' data(hash_ref)
#' times <- attr(hash_ref, "time.series")
#' gene1 <- hash_ref$X[500,]
#' \donttest{plot(times, gene1, type='l')}
"hash_ref"