#' Interpolation of gene expression on time series from reference data
#' 
#' This function computes the interpolated gene expression data from a reference time series.
#' This is done with a multi-target Partial Least Square Rregression (PLSR) model, using a spline of time as descriptive variables.
#' 
#' @param X gene expression matrix of reference time series, genes as rows, (ordered) individuals as columns.
#' @param time.series timepoints of the reference (X).
#' @param df the df parameter passed on to the \code{\link[splines]{ns}()} function.
#' @param n.inter number of timepoints to return in interpolated data, defaults to 200.
#' @param tmin,tmax defaults to min and max of \code{time.series} ; start and end times of interpolated time series.
#' @param scale defaults to TRUE, passed on to the \code{\link[pls]{plsr}()} function.
#' @param knots defaults to NULL, passed on to the \code{\link[splines]{ns}()} function.
#' @param return.model if TRUE, returns the PLSR model object and the input df value.
#' 
#' @export
#' 
#' @examples 
#' 
#' data("Cel_embryo")
#' 
#' iGE <- plsr_interpol(Cel_embryo$X, Cel_embryo$time.series, df = 9, n.inter = 100)
#' 
#' \donttest{
#' # plot random genes with their interpolations
#' par(mfrow=c(2,2))
#' invisible(sapply(sample(1:nrow(Cel_embryo$X), 4), function(i){
#'   plot(Cel_embryo$time.series, Cel_embryo$X[i,], lwd = 2, 
#'        xlab = "time", ylab = "gene expression", 
#'        main = rownames(Cel_embryo$X)[i])
#'   points(iGE$time.series, iGE$interpGE[i,], lwd = 2,
#'          type = "l", col = "royalblue")
#' }))
#' }
#'
#' 
#' @importFrom stats predict
#' @importFrom splines ns
#' @import pls 

plsr_interpol <- function(X, time.series, df, 
                          n.inter = 200, 
                          tmin = min(time.series), tmax = max(time.series),
                          scale = T, knots = NULL,
                          return.model = FALSE)
{
  if(length(time.series)!=ncol(X)){
    stop("time series must be of length ncol(X)")
  }
  if(tmin < min(time.series) | tmax > max(time.series)){
    stop("tmin and tmax must be within time.series")
  }
  if(n.inter<ncol(X)){
    warning("n.inter should be larger than ncol(X)")
  }

  
  if(df>=(length(time.series)))
    stop("Spline degree must be less than number of unique samples")
  
  t.mat <- splines::ns(time.series, df=df, knots=knots) # Build polynomial from time series
  m.X <- pls::plsr(t(X)~t.mat, scale=scale, validation='CV') # apply plsr model
  
  # get appropriate ncomp for plsr with built-in CV
  cv <- pls::RMSEP(m.X)$val['CV',,]
  nc <- which.min(colMeans(cv)) - 1
  if(nc==0)
    nc <- which.min(colMeans(cv)[-1])
  
  # Build interpolation
  inter <- seq(tmin, tmax, length.out = n.inter)
  
  inter.mat <- stats::predict(t.mat, inter) # make interpolated ns time series 
  pred.mat <- stats::predict(m.X, newdata = inter.mat, comps=1:nc) # use model to predict gene expr
  
  
  res <- list(
    time.series = inter,   # time series (interpolated)
    interpGE = t(pred.mat) # model predictions (interpolated)
    )
  if(isTRUE(return.model)){
    res$plsr.model <- m.X  # PLSR model
    res$df <- df           # df param
  }
  return(res)            
} 
