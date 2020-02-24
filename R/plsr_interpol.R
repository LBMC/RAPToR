#' Interpolation of gene expression on time series from reference data
#' 
#' *NOTE : this function is deprecated, it is recommended to use ge_im functions to interpolate on gene expression data*
#' 
#' This function computes the interpolated gene expression data from a reference time series.
#' This is done with a multi-target Partial Least Square Regression (PLSR) model, using a spline of time as descriptive variables.
#' 
#' @param X gene expression matrix of reference time series, genes as rows, (ordered) individuals as columns.
#' @param time.series timepoints of the reference (`X`).
#' @param df the df parameter passed on to the \code{\link[splines]{ns}} function.
#' @param covar a covariate to include in the model (*e.g* batch).  
#' @param topred a level of `covar` to use for model predictions ; defaults to the first level.
#' @param n.inter number of timepoints to return in interpolated data, defaults to 200.
#' @param tmin,tmax defaults to min and max of `time.series` ; start and end times of interpolated time series.
#' @param scale defaults to TRUE, passed on to the \code{\link[pls]{plsr}} function.
#' @param knots defaults to NULL, passed on to the \code{\link[splines]{ns}} function.
#' @param plsr.nc the number of components to use for PLSR prediction. If NULL, determined by CV.
#' @param return.model if TRUE, returns the PLSR model object and the input df value.
#' 
#' @export
#' 
#' @eval interpol_example()
#' 
#' @importFrom stats predict formula
#' @importFrom splines ns
#' @import pls 

plsr_interpol <- function(X, time.series, df, 
                          covar = NULL, topred = NULL,
                          n.inter = 200, 
                          tmin = min(time.series), tmax = max(time.series),
                          scale = T, knots = NULL, plsr.nc = NULL,
                          return.model = FALSE)
{
  if (length(time.series) != ncol(X)) {
    stop("time series must be of length ncol(X)")
  }
  if (tmin < min(time.series) | tmax > max(time.series)) {
    stop("tmin and tmax must be within time.series")
  }
  if (n.inter < ncol(X)) {
    warning("n.inter should be larger than ncol(X)")
  }
  if (df >= (length(time.series))) 
    stop("Spline degree must be less than number of unique samples")
  
  
  # for building the model
  dat <- data.frame(
    X = I(t(X)),
    Y = I(splines::ns(time.series, df = df, knots = knots))
  )
  
  m_formula <- stats::formula(X ~ Y)
  
  # for predictions
  inter <- seq(tmin, tmax, length.out = n.inter)
  ndat <- data.frame( 
    Y = I(stats::predict(dat$Y, inter))
  )
  
  
  if(!is.null(covar)){
    if (length(covar) != length(time.series)){
      stop("covar should be of length ncol(X)")
    }
    
    covar <- factor(covar)
    
    if(is.null(topred)){
      # interpolate as the first dataset/condition by default
      topred <- covar[which(covar == levels(covar)[1])[1]]
    }
    else{
      # topred must be specified as a level of covar
      if(! topred %in% levels(covar)){
        stop("topred must be a level of covar")
      }
      topred <- covar[which(covar == topred)[1]]
    }
    
    dat$covar <- covar
    m_formula <- stats::formula(X ~ Y + covar)
    ndat$covar <- rep(topred, n.inter)
    
  } 
  
  m.X <- pls::plsr(m_formula, data = dat, scale = scale, validation = "CV")
  
  # choose nc for plsr interpol
  if(is.null(plsr.nc)){
    cv <- pls::RMSEP(m.X)$val["CV", , ]
    nc <- which.min(colMeans(cv)) - 1
    if (nc == 0){
      nc <- which.min(colMeans(cv)[-1])
    }
  } else {
    nc <- plsr.nc
  }
  
  pred <- stats::predict(m.X, newdata = ndat, comps = 1:nc)
  
  # format results
  res <- list(time.series = inter, interpGE = t(pred))
  
  if (isTRUE(return.model)) {
    res$plsr.model <- m.X
    res$plsr.nc <- nc
    res$df <- df
  }
  return(res)
} 
