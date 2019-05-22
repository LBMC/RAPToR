#' Compute correlation between sample and reference data
#' 
#' Computes the correlation between the gene expressions of the sample 
#' and the reference data using the \code{\link{cor}} function ;
#' default correlation index is spearman.
#' 
#' @param samp the sample matrix, gene as rows, individuals as columns
#' @param refdata the reference time series matrix, same format as \code{samp}
#' @param cor.method \code{method} parameter passed on to the \code{\link{cor}} function
#' 
#' @return a matrix of correlation scores with samples as columns, refdata as rows
#' 
#' @export
#' 
#' @examples 
#' data(Cel_larval)
#' 
#' samp <- Cel_larval$X[,c(3,8,12,20)]
#' cc <- cor.gene_expr(samp, Cel_larval$X)
#' \donttest{
#' plot(cc, margins=c(10,5))
#' }
#' 
#' @importFrom stats cor na.omit
#' @importFrom data.table frank
#' 
cor.gene_expr <- function(samp, refdata, cor.method="spearman")
{
  if(all(rownames(samp)!=rownames(refdata))){
    stop("Sample and reference matrices do not contain the same gene set.")
  }
  if(any(is.na(samp))){
    na.rows <- nrow(samp)-nrow(stats::na.omit(samp))
    samp <- na.omit(samp)
    refdata <- refdata[rownames(samp),]
    warning(paste("NA values in samp, removed", na.rows, "rows with NA."))
  }
  if(cor.method=="spearman"){
    # using a combination of data.table's frank
    # and R's cor is much faster than the default 
    # spearman method.
    cors <- stats::cor(apply(refdata, 2, data.table::frank), 
                       apply(samp, 2, data.table::frank), 
                       method="pearson")
  }
  else{
    cors <- stats::cor(refdata, samp, method=cor.method)
  }
  if(!is.matrix(cors)){
    cors <- as.matrix(cors)
  }
  
  return(cors)
}





