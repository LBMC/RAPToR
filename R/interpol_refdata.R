#' function that fits a loess curve on an ICA component and returns predictions on given x series
#' 
#' @param ICA an ICA object (returned by icafast)
#' @param time.series the initial time series
#' @param comp the index of the ICA component on which to fit the curve
#' @param pred.x the new time series to predict component values on
#' @param span span argument given to loess
#' @param plot if TRUE, plot the component and its fitted curve
#' 
#' @export
#' 
#' @examples
get.spline <- function(ICA, time.series, comp, pred.x,
                       span=.25, plot=T)
{
  r <- loess(ICA$M[,comp]~time.series, span = span)
  pred.y <- predict(r, newdata=pred.x)
  
  if(plot==T){
    plot(ICA$M[,comp]~time.series, main=paste("Comp.", comp))
    lines(pred.x, pred.y, lwd=2)
  }
  return(list(reg=r, curve=data.frame(x=pred.x, y=pred.y)))
}



#' function that computes interpolated data from a reference time series
#' 
#' @param X gene expression matrix of reference time series, genes as rows, (ordered) individuals as columns
#' @param n.inter number of timepoints to return in interpolated data
#' @param time.series timepoints of the reference (X). If none given, 1:ncol(X) is used.
#' @param ica.nc number of components to keep in icafast
#' @param keep.c indices of components to keep for interpolation
#' @param t.min start time of new time series (ignored if new.timepoints is given)
#' @param t.max end time of new time series (ignored if new.timepoints is given)
#' @param new.timepoints vector of length n.inter with the new time series (overrides t.min and t.max)
#' @param span value(s) given to loess for curve fitting on the ica components
#' @param plot if TRUE, plots the components and their fitted curves
#' @param return.fits if TRUE, returns the list of fitted loess objects and curves
#' 
#' @export
#' 
#' @examples
interpol_refdata <- function(X, n.inter,
                             time.series=NULL,
                             ica.nc=16, keep.c=1:10,
                             t.min=NULL, t.max=NULL, new.timepoints=NULL,
                             span=0.25, plot=F, return.fits=F)
{
  if(n.inter<ncol(X)){
    stop("n.inter must be larger than ncol(X)")
  }
  if(!is.null(time.series)&length(time.series)!=ncol(X)){
    stop("time series must be of length ncol(X)")
  }
  if(is.null(time.series)){
    warning("no time.series given, using 1:ncol(X) as reference")
    time.series <- 1:ncol(X)
  }
  if(length(keep.c)>ica.nc|any(!keep.c%in%(1:ica.nc))){
    stop("keep.c must be included in 1:ica.nc")
  }
  if(t.min==t.max&is.null(t.min)&is.null(new.timepoints)){
    warning("No time scale given, setting t.min as 1 and t.max as ncol(X)")
    t.min <- 1
    t.max <- ncol(X)
  }
  if(!is.null(new.timepoints)&length(new.timepoints)!=n.inter){
    stop("new.timepoints must be of length n.inter")
  }
  if(!is.null(new.timepoints)){
    if(!is.null(t.min)|!is.null(t.max)){
      warning("ignoring t.min and t.max, using new.timepoints as reference")
    }
  }
  if(!is.null(t.min)&!is.null(t.min)){
    new.timepoints <- seq(t.min, t.max, length.out = n.inter)
  }
  
  if(length(span)<length(keep.c)){
    span <- rep(span, length.out=length(keep.c))
    names(span) <- keep.c
  }
  
  require(ica)
  # compute ICA
  ICA <- icafast(X, ica.nc)
  
  # get splines and predictions
  rs <- lapply(keep.c, function(i){
    get.spline(ICA, time.series, i, new.timepoints,
               span = span[i], plot=plot)
  })
  
  # extract interpolated individual data
  interpol.indivs <- as.matrix(sapply(1:length(keep.c), function(i){
    rs[[i]]$curve$y
  }))
  
  # compute interpolated gene expression data
  interpol.gene_expr <- tcrossprod(ICA$S[,keep.c], interpol.indivs)
  colnames(interpol.gene_expr) <- new.timepoints
  
  
  if(return.fits){
    return(list(interpol.gene_expr = interpol.gene_expr, 
                time.series = new.timepoints,
                fits = rs))
  }
  
  return(list(interpol.gene_expr = interpol.gene_expr, 
              time.series = new.timepoints))
  
  
}