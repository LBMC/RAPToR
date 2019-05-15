#' Sterken reference data time series
#'
#' Data from the ArrayExpress submitted by Mark G. Sterken in 2014, 
#' recently published by Van Sluijs et al. (2019).
#' Only samples of Standard N2 strain were kept here.
#' 
#' Time is in hours from hatching, scaled on the Oudenaarden 20C time series, 
#' re-estimated on the Reinke and Oudenaarden references.
#' 
#' 
#'
#' @docType data
#'
#' @usage data(sterken_ref)
#'
#' @format a list with \code{X} being the gene expression matrix and \code{time.series} being the timepoints (in hours post-hatching)
#'
#' @keywords datasets
#'
#' @references Van Sluijs, et al.  (2019) BioRxiv : 579151
#' (\href{https://www.biorxiv.org/content/10.1101/579151v3}{bioRxiv})
#'
#' @source \href{https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-7574/}{ArrayExpress dataset}
#'
#' @examples
#' data(sterken_ref)
#' times <- sterken_ref$time.series
#' gene1 <- sterken_ref$X[700,]
#' \donttest{plot(times, gene1, type='b')}
"sterken_ref"