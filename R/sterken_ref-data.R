#' C. elegans Young Adult/Adult reference data time series (1)
#'
#' Data from the ArrayExpress submitted by Mark G. Sterken in 2014
#' (recently published by Van Sluijs et al. (2019)), with C elegans gene expression
#' levels (log+1 from intensities).
#' 
#' Only samples of Standard N2 strain were kept here.
#' 
#' Time is in hours from hatching, scaled on the Oudenaarden 20C time series, 
#' re-estimated on the Reinke and Oudenaarden references.
#' 
#' 
#'
#' @docType data
#'
#' @usage data(Cel_YA_adult1)
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
#' data(Cel_YA_adult1)
#' times <- Cel_YA_adult1$time.series
#' gene1 <- Cel_YA_adult1$X[700,]
#' \donttest{plot(times, gene1, type='b')}
"Cel_YA_adult1"