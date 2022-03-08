#' Create a reference object
#' 
#' Make a reference object from a geim model
#' 
#' 
#' @param m a geim model object (as returned by \code{\link{ge_im}}).
#' @param from,to start/end of the interpolation (defaults to first and last time points) 
#' @param cov.levels a list with potential model covariate levels (e.g batch, strain) to predict as (defaults to first level). If multiple, give in the formula order.
#' @param n.inter interpolation resolution, as in seq(start, end, length.out = n.inter). One of \code{n.inter} or \code{by.inter} must be specified.
#' @param by.inter interpolation resolution, as in seq(start, end, by = by.inter). One of \code{n.inter} or \code{by.inter} must be specified.
#' @param t.unit an optional string specifying the time unit and t-zero, e.g "h past egg-laying".
#' @param ... extra arguments passed on to model functions.
#'
#' @return a '\code{ref}' object to use with ae() for age estimation.
#' 
#' @export
#' 
#' @eval interpol_example()
#' 
#' @importFrom mgcv interpret.gam
#' @importFrom Rdpack reprompt
#' 
make_ref <- function(m,
                     from = NULL, to = NULL,
                     cov.levels = NULL,
                     n.inter = 500, by.inter = NULL,
                     t.unit = "",
                     ...)
{
  p <- attr(m, "pdata")
  f <- as.formula(attr(m, "formula"))
  vars <- mgcv::interpret.gam(f)$pred.names
  # extract time variable
  t.var <- p[,vars[which(is.numeric(p[,vars])[1])]]
  
  if(length(vars) > 1){ # there is a covariate
    if(is.null(cov.levels)){ # no specified level
      warning("No covariate level specified: using first level fro interpolation.")
      cov.levels <- lapply(vars[-1], function(i) levels(p[, i])[1]) # get 1st level of each covariate
    } else{
      if(length(cov.levels) != (length(vars) -1)){
        stop("Levels for all (or none of) the covariates must be specified, in the formula order.")
      }
    }
  }
  
  # from & to param handling 
  if(!is.null(from)){
    if(!is.numeric(from) | from < min(t.var)){
      stop("'from' must be a numeric value within the reference time span.")
    }
  } else {
    from <- min(t.var)
  }
  
  if(!is.null(to)){
    if(!is.numeric(to) | to > max(t.var)){
      stop("'to' must be a numeric value within the reference time span.")
    }
  } else {
    to <- max(t.var)
  }
  
  
  ndat <- data.frame(time = seq(min(tvar), max(tvar), l = n.inter))
  
}