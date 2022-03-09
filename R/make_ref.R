#' Create a reference object
#' 
#' Make a reference object from a geim model
#' 
#' 
#' @param m a geim model object (as returned by \code{\link{ge_im}}).
#' @param from,to start/end of the interpolation (defaults to first and last time points) 
#' @param cov.levels a named list with potential model covariate levels (e.g batch, strain) to predict as (defaults to first level).
#' @param n.inter interpolation resolution, as in seq(start, end, length.out = n.inter). One of \code{n.inter} or \code{by.inter} must be specified.
#' @param by.inter interpolation resolution, as in seq(start, end, by = by.inter). One of \code{n.inter} or \code{by.inter} must be specified.
#' @param t.unit an optional string specifying the time unit and t-zero, e.g "h past egg-laying".
#' @param metadata an optional named list with reference metadata (e.g. organism, tissue).
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
                     metadata = list())
{
  p <- attr(m, "pdata")
  f <- attr(m, "formula")
  vars <- mgcv::interpret.gam(as.formula(f))$pred.names
  
  # extract time variable
  ti <- which(sapply(vars, function(v) is.numeric(p[,v])))[1]
  t.var <- p[,vars[ti]]
  cvars <- vars[-ti]
  
  # cov.levels param handling
  if(length(cvars) > 0){ # there is a covariate
    if(is.null(cov.levels)){ # no specified level
      cov.levels <- lapply(cvars, function(i) levels(p[, i])[1]) # get 1st level of each covariate
      names(cov.levels) <- cvars
      wl <- paste(sapply(seq_along(cvars), 
                          function(i) paste0("\n\t", cvars[i], ": ", cov.levels[[i]])), 
                   collapse = "")
      wl <- paste0("No covariate level specified: using first level for interpolation.\n",
                   "  Model formula:\n\t", as.character(f), "\n  Covariate level(s) used:", wl)
      warning(wl)
    } else{
      if(!is.list(cov.levels) | length(cov.levels) != length(cvars)){ # not all cov levels specified
        stop("Levels to predict for all the covariates must be specified as a named list (with variable names).")
      }
      if(!all(sapply(seq_along(cvars), function(i) cov.levels[i] %in% levels(p[,vars[i]])))){
        stop("A level of cov. levels does not match with factor levels in the data.")
      }
    }
  } else {
    cov.levels <- list() # no covariate, no cov.levels.
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
  
  # n/by.inter param handling
  if(!is.null(by.inter)){
    if(!is.numeric(by.inter)){
      stop("by.inter must be a single numeric value.")
    }
    ts <- seq(min(t.var), max(t.var), by = by.inter)
    l <- length(ts)
  } else {
    ts <- seq(min(t.var), max(t.var), by = n.inter)
    l <- n.inter
  }
  
  # t.unit param handling
  
  
  # make the new predictor dataframe
  ndat <- data.frame(time = ts)
  colnames(ndat)[1] <- vars[ti]
  
  for(v in cvars){
    ndat[, v] <- factor(rep(cov.levels[[v]], l=l))
  }
  
  # build ref. object
  ref <- list(interpGE = predict(m, ndat),
              time = ts)
  
  # assign object class and attributes
  class(ref, 'ref')
  attr(ref, "t.unit") <- t.unit
  attr(ref, "metadata") <- metadata
  attr(ref, "formula") <- f
  attr(ref, "cov.level") <- cov.levels
  
  return(ref)
  
}