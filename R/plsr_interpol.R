#' Interpolation of gene expression on time series from reference data
#' 
#' This function computes the interpolated gene expression data from a reference time series.
#' This is done with a multi-target PLSR model
#' 
#' @param X gene expression matrix of reference time series, genes as rows, (ordered) individuals as columns
#' @param n.inter number of timepoints to return in interpolated data
#' @param time.series timepoints of the reference (X).
#' @param t.min,t.max defaults to min and max of \code{time.series} ; start and end times of new time series ; ignored if new.timepoints is given.
#' @param new.timepoints vector of length n.inter with the new time series (overrides t.min and t.max)
#' @param return.model if TRUE, returns the PLSR model object
#' 
#' @export
#' 
#' @examples 
#' 
#' \donttest{
#' data(Cel_larval)
#' par(mfrow=c(2,2))
#'
#' interpold <- interpol_refdata(X = Cel_larval$X, n.inter = 200, 
#'                               time.series = Cel_larval$time.series, ,
#'                               ica.nc = 10, keep.c = 1:10,
#'                               plot = TRUE)
#'
#' pb <- sapply(c(2,5,13,50), function(i){
#'    plot(Cel_larval$time.series, Cel_larval$X[i,],
#'         type = 'l', lwd=2, 
#'         main=rownames(Cel_larval$X)[i], xlab = 'time')
#'    points(interpold$time.series, interpold$interpol.gene_expr[i,],
#'           type = 'l', lwd=2, col='firebrick')
#'   if(i==5){
#'     legend('topright', legend = c("Initial ref. data", "Interpolated data"),
#'            lwd=3, col=c('black', 'firebrick'), bty = "n")
#'   }
#' })
#' }
#'
#' 
#' @importFrom stats predict
#' @importFrom splines ns
#' @importFrom pls plsr RMSEP

plsr_interpol <- function(X, time.series, df, 
                          n.inter=100, scale=T, 
                          knots = NULL)
{
  if(n.inter<ncol(X)){
    stop("n.inter must be larger than ncol(X)")
  }
  if(length(time.series)!=ncol(X)){
    stop("time series must be of length ncol(X)")
  }

  
  if(df>=(length(time.series)))
    stop("Spline degree must be less than number of unique samples")
  
  t.mat <- splines::ns(time.series, df=df, knots=knots) # Build polynomial from time series
  m.X <- pls::plsr(t(X)~t.mat, scale=scale, validation='CV') # apply plsr model
  
  # get appropriate ncomp for plsr with built-in CV
  cv <- RMSEP(m.X)$val['CV',,]
  nc <- which.min(colMeans(cv)) -1
  if(nc==0)
    nc <- which.min(colMeans(cv)[-1])
  
  # Build interpolation
  inter <- seq(min(time.series), max(time.series), length.out = n.inter)
  inter.mat <- stats::predict(t.mat, inter) # make interpolated ns time series 
  pred.mat <- stats::predict(m.X, newdata = inter.mat, comps=1:nc)
  
  return(list(plsr.model=m.X,    # PLSR model
              ts=inter,          # time series (interpolated)
              pred=t(pred.mat),  # model predictions (interpolated)
              df=df))            # spline df
  
} 
