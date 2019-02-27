#' Compute correlation between sample and reference data
#' 
#' Computes the correlation between the gene expressions of the sample 
#' and the reference data using the \code{\link{cor}} function.
#' 
#' @param samp the sample matrix, gene as rows, individuals as columns
#' @param refdata the reference time series matrix, same format as \code{samp}
#' @param cor.method \code{method} parameter passed on to the \code{\link{cor}} function
#' 
#' @return a \code{corg} object, being a matrix of correlation scores with samples as columns, refdata as rows
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
cor.gene_expr <- function(samp, refdata, cor.method="pearson")
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
    apply(refdata, 2, stats::cor, y=cb, method=cor.method)
  })
  class(cors) <- 'corg'
  return(cors)
}



#' Estimate the developmental age of individuals
#' 
#' This function estimates the developmental age of individuals based on given
#' reference data as well as an approximate time estimate of the age.
#' This is achieved using a two-step process.
#' First, simple correlation between the gene expression profiles of
#' the sample and reference data are computed using \code{\link{cor.gene_expr}}; 
#' this gives correlation profiles according to reference time per individual. 
#' Then, the maxima of the correlation profiles are evaluated and scored according to 
#' a gaussian distribution around the time estimate given as input. 
#' The maxima scores range from 0 to 1, 1 being the case where a correlation peak is 
#' exactly at the given approximate time. 
#' 
#' @param samp the sample matrix, gene as rows, individuals as columns
#' @param refdata the reference time series matrix, same format as \code{samp}
#' @param ref.time_series the reference time series (\textit{e.g.} \code{colnames(refdata)} if using interpolated reference data)
#' @param est.time the approximate development time of the samples, must be in the same units than \code{ref.time_series}
#' @param time.sd the std. deviation of the gaussian scoring distribution. \textit{Note that setting this value too low can cause a significant bias in the age estimation.}
#' @param cor.method correlation method argument passed on to \code{\link{cor.gene_expr}}
#' @param all.peaks logical; if TRUE, returns all correlation peaks (potential age estimates) and their respective scores for every individuals, as a list. If FALSE, only returns the best estimate for each individual, as a dataframe.
#' 
#' 
#' @return an '\code{ae}' object, which is a list of the correlation matrix between sample and reference, and the age estimates (either as list of individuals or dataframe, depending on \code{all.peaks}).
#' 
#' @export
#' 
#' @examples 
#' data(oud_ref)
#' 
#' samp <- oud_ref$X[,13:15]
#' age.est <- estimate.worm_age(samp, oud_ref$X, oud_ref$time.series, 26)
#' age.est$age.estimates
#' \donttest{
#' }
#' 
#' 
estimate.worm_age <- function(samp, refdata, ref.time_series, est.time,
                              time.sd=5, cor.method="pearson", all.peaks=F)
{
  if(length(ref.time_series)!=ncol(refdata)){
    stop("Reference data and time series don't match")
  }
  if(est.time>max(ref.time_series)|est.time<min(ref.time_series)){
    stop("Estimated time is outside reference time series' range")
  }
  ref.time_series <-  as.numeric(ref.time_series)
  
  # compute correlations
  cors <- cor.gene_expr(samp, refdata, cor.method = cor.method)
  
  # compute gaussian around estimated time
  ref.gauss <- dnorm(ref.time_series, mean = est.time, sd = time.sd)
  m.gauss <- max(ref.gauss)
  
  age.estimates <- lapply(1:ncol(samp), function(i){
    # get correlation maxima (peaks) positions
    cor.maxs.i <- which(diff(sign(diff(cors[,i])))==-2)+1
    cor.maxs <- cors[cor.maxs.i, i]
    cor.maxs.times <- ref.time_series[cor.maxs.i]
    
    #compute scores based on gaussian of reference time
    cor.maxs.scores <- round(ref.gauss[cor.maxs.i]/m.gauss, 4)
    
    age.estimate <- cbind(time=cor.maxs.times, 
                            cor.score=cor.maxs,
                            proba.score=cor.maxs.scores)
    
    age.estimate <- age.estimate[order(age.estimate[, "proba.score"], decreasing = T),]
    
    return(age.estimate)
  })
  
  names(age.estimates) <- colnames(samp)
  if(!all.peaks){
    # only keep best estimate
    age.estimates <- sapply(age.estimates, function(a.e){a.e[1,]})
  }
  res <- list(cors=cors, age.estimates=age.estimates)
  class(res) <- "ae"
  return(res)
  
  
}










#' Plot a corg object
#' 
#' Plots the correlation matrix returned by \code{\link{cor.gene_expr}} 
#' as a heatmap.
#' 
#' @param cors a \code{corg} object, as returned by \code{\link{cor.gene_expr}} 
#' @param col a color scheme, passed on to \code{\link{heatmap}}
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
  requireNamespace("graphics", quietly = T)
  requireNamespace("grDevices", quietly = T)
  heatmap(cors, Rowv=NA, Colv=NA, 
          col = col, scale="column", ...)
}
