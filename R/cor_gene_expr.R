#' Compute correlation between sample and reference data
#' 
#' Computes the correlation between the gene expressions of the sample 
#' and the reference data using the \code{\link{cor}} function.
#' 
#' @param samp the sample matrix, gene as rows, individuals as columns
#' @param refdata the reference time series matrix, same format as \code{samp}
#' @param method \code{method} parameter passed on to the \code{\link{cor}} function
#' 
#' @return a \code{corg} object, being a matrix of correlation scores samples as columns, refdata as rows
#' 
#' @export
#' 
#' @examples 
#' data(oud_ref)
#' 
#' samp <- oud_ref$X[,c(3,8,12,20)]
#' cc <- cor.gene_expr(samp, oud_ref$X)
#' \donttest{
#' plot(cc, margins=c(10,5))
#' }
#' 
cor.gene_expr <- function(samp, refdata, method="pearson")
{
  requireNamespace("stats", quietly = T)
  if(all(rownames(samp)!=rownames(refdata))){
    stop("Sample and reference matrices do not contain the same gene set.")
  }
  if(any(is.na(samp))){
    na.rows <- nrow(samp)-nrow(na.omit(samp))
    samp <- na.omit(samp)
    refdata <- refdata[rownames(samp),]
    warning(paste("NA values in samp, removed", na.rows, "rows with NA."))
  }
  cors <- apply(samp,2, function(cb){
    apply(refdata, 2, stats::cor, y=cb, method=method)
  })
  class(cors) <- 'corg'
  return(cors)
}



#' Plot a corg object
#' 
#' Plots the correlation matrix returned by \code{\link{cor.gene_expr}} 
#' as a heatmap.
#' 
#' @param cors a \code{corg} object, as returned by \code{\link{cor.gene_expr}} 
#' @param ... additional arguments passed on to \code{\link{heatmap}}
#' 
#' @export
#' 
#' @examples
#' data(oud_ref)
#' 
#' samp <- oud_ref$X[,c(3,8,12,20)]
#' cc <- cor.gene_expr(samp, oud_ref$X)
#' \donttest{
#' plot(cc, margins=c(10,5))
#' }
#' 
plot.corg <- function(cors, col=cm.colors(256), ...)
{ 
  heatmap(cors, Rowv=NA, Colv=NA, 
          col = col, scale="column", ...)
}
