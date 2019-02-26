library(ica)

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