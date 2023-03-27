#' Create a reference object
#' 
#' Make a reference object from a geim model
#' 
#' 
#' @param m a geim model object (as returned by \code{\link{ge_im}}).
#' @param from,to start/end of the interpolation (defaults to first and last time points) 
#' @param cov.levels a named list with potential model covariate levels (e.g batch, strain) to predict as (defaults to first level).
#' @param n.inter interpolation resolution, as in \code{seq(start, end, length.out = n.inter)}. One of \code{n.inter} or \code{by.inter} must be specified.
#' @param by.inter interpolation resolution, as in \code{seq(start, end, by = by.inter)}. One of \code{n.inter} or \code{by.inter} must be specified.
#' @param t.unit an optional string specifying the time unit and t-zero, e.g "h past egg-laying".
#' @param metadata an optional named list with reference metadata (e.g. organism, tissue).
#'
#' @return a '\code{ref}' object to use with \code{\link{ae}} for age estimation.
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
                     t.unit = "no unit specified",
                     metadata = list())
{
  cov.levels <- .cov_check(m, cov.levels) # param check 
  
  p <- attr(m, "pdata")
  t.var <- attr(m, "vars")$t.var # time variable
  tv <- p[, t.var]
  cvars <- attr(m, "vars")$cvars # covariates
  
  
  # from & to param handling 
  if(!is.null(from)){
    if(!is.numeric(from) | from < min(tv)){
      stop("'from' must be a numeric value within the reference time span.")
    }
  } else {
    from <- min(tv)
  }
  if(!is.null(to)){
    if(!is.numeric(to) | to > max(tv)){
      stop("'to' must be a numeric value within the reference time span.")
    }
  } else {
    to <- max(tv)
  }
  
  # n/by.inter param handling
  if(is.null(by.inter) & is.null(n.inter)){
    stop("one of n.inter or by.inter must be specified")
  }
  if(!is.null(by.inter)){
    if(!is.numeric(by.inter)){
      stop("by.inter must be a single numeric value.")
    }
    ts <- seq(from, to, by = by.inter)
    l <- length(ts)
  } else {
    ts <- seq(from, to, l = n.inter)
    l <- n.inter
  }
  
  # t.unit param handling
  if(!is.character(t.unit)){
    stop("t.unit must be a string.")
  }
  # metadata param handling
  if(0 != length(metadata)){
    if(!is.list(metadata)){
      stop("metadata must be given as a named list.")
    }
    if(is.null(names(metadata))){
      names(metadata) <- rep("unnamed", length(metadata))
    }
  } else {
    metadata <- list("no metatdata"="")
  }
  
  # make the new predictor dataframe
  ndat <- data.frame(time = ts)
  colnames(ndat)[1] <- t.var
  
  for(v in cvars){
    ndat[, v] <- factor(rep(cov.levels[[v]], l=l))
  }
  
  # build ref. object
  ref <- list(interpGE = predict.geim(m, ndat),
              time = ts)
  
  # get geim params
  gp <- c("formula", "method", "dim_red", "nc")
  geim.params <- lapply(gp, attr, x=m)
  names(geim.params) <- gp
  
  # assign object class and attributes
  class(ref) <- 'ref'
  attr(ref, "t.unit") <- t.unit
  attr(ref, "metadata") <- metadata
  attr(ref, "geim.params") <- geim.params
  attr(ref, "cov.levels") <- cov.levels
  
  return(ref)
}




#' Print a ref object
#' 
#' Prints a \code{ref} object
#' 
#' @param x a \code{ref} object, as returned by \code{\link{make_ref}}.
#' @param ... arguments passed on to \code{\link{print}}
#' 
#' @export
#' 
print.ref <- function(x, ...){
  d <- dim(x$interpGE)
  ts <- diff(x$time)
  cat("RAPToR reference object\n---")
  cat("\n interpGE:\t( ",d[1], " x ", d[2]," )",
      "\n time:\t [ ", min(x$time), " - ", max(x$time), " ] ", attr(x, "t.unit"),", by ",
      ifelse(all(sapply(ts, all.equal, current=ts[1])), ts[1], "varying time"), " steps",
      sep = "")
  
  md <- attr(x, "metadata")
  if(1 == length(md) & md[[1]]==""){
    cat("\n\n (no metadata)")
  } else {
  cat('\n\n Metadata:\n\t', 
      paste(sapply(seq_along(md), function(i) paste0(names(md)[i], ": ", md[[i]])),
            collapse = "\n\t"),
      sep = "")
  }
  
  cat("\n\n GEIM:")
  ats <- attr(x, "geim.params")
  if(ats$method != "limma"){
    cat("\n\t", casefold(ats$method, upper = T), "fit on", ats$nc,
        casefold(ats$dim_red, upper = T), "components:")
  } else {
    cat("\n\t", casefold(ats$method, upper = T), " model fit on genes with:")
  }
  cat("\n\t", ats$formula)
  # cov. levels
  cl <- attr(x, "cov.levels")
  if(length(cl)>0){
    
    cat("\n\t with covariate levels", 
        paste(sapply(seq_along(cl), 
                     function(i) paste0(names(cl)[i], ": ", cl[[i]])),
              collapse = ", "))
  }
 
  cat("\n---\n")
}


.cov_check <- function(m, cov.levels){
  p <- attr(m, "pdata")
  f <- attr(m, "formula")
  cvars <- attr(m, "vars")$cvars # covariates

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
      if(!all(sapply(seq_along(cvars), function(i) cov.levels[[i]] %in% levels(p[, cvars[i]])))){
        stop("A level of 'cov.levels' does not match with factor levels in the data.")
      }
    }
  } else {
    cov.levels <- list() # no covariate, no cov.levels.
  }
  return(cov.levels)
}
