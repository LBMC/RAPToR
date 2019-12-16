#' Find the optimal spline df through cross-valdation
#' 
#' This function performs cross-validation (CV) with the aim of finding the optimal spline df parameter 
#' for a multi-target PLSR regression.
#' By default, CV is done using the sample loadings on ICA components rather than the full gene expression matrix as target variable.
#' The reason for this is a non-negligible reduction of computing costs (>100x faster) for equal results. 
#' Using components as "eigen genes" is not uncommon to find model parameters fitting the whole gene set.
#' 
#' @param X gene expression matrix of reference time series, genes as rows, (ordered) individuals as columns.
#' @param time.series timepoints of the reference (`X`).
#' @param covar a covariate to include in the model (*e.g* batch).  
#' @param dfs vector of spline df parameters to scan (passed on to \code{\link[splines]{ns}} function).
#' @param cv.n number of cross-validation repeats.
#' @param cv.s ratio of samples to use for training set. If `cv.s > 1`, then `cv.s` samples are used for the training set.
#' @param err.func the error function to use to compute the CV error. Defaults to sum square of differences (the \code{\link{ef}} function).
#' @param ica.use boolean ; if TRUE (default), sample loadings of ICA components are used for CV rather than the whole gene expression matrix.
#' @param ica.nc number of components to use for the ica. Defaults to 16 or `ncol(X)` if there are less than 16 samples.
#' @param nb.cores number of cores to use for parallelization.
#' 
#' @return a '`dfCV`' object, which is a list of the CV error table (`cv.n` x `dfs`), the dfs, 
#' the CV parameters and PLSR model components used for prediction
#' There are `plot`, `print` and `summary` methods for this object.
#' 
#' @export
#' 
#' @eval interpol_example()
#' 
#' @importFrom parallel parLapply stopCluster makeForkCluster
#' @importFrom stats quantile dnorm predict formula
#' @importFrom splines ns
#' @importFrom ica icafast
#' @import pls 
#' 
#' 
df_CV <- function(X, time.series, 
                  dfs = 1:(ncol(X)-1), 
                  covar = NULL,
                  cv.n = 50, cv.s = 0.8, 
                  err.func = ef,
                  ica.use = TRUE, ica.nc = min(ncol(X), 16),
                  nb.cores = 2){

  y <- time.series
  n <- length(y)
  if (isTRUE(ica.use)) {
    X <- ica::icafast(X, nc = ica.nc)$M
  }
  else {
    X <- t(X)
  }
  if (cv.s < 1) 
    n_samp <- round(cv.s * n, 0)
  else n_samp <- cv.s
  
  # for building the model
  dat <- data.frame(
    X = I(X)
  )
  m_formula <- stats::formula(X ~ Y)
  
  if(!is.null(covar)){
    if (length(covar) != length(time.series)){
      stop("covar should be of length ncol(X)")
    }
    covar <-  factor(covar)
    
    dat$covar <- covar
    m_formula <- stats::formula(X ~ Y + covar)
    
    # need to keep at least first and last points from each dataset/condition in training set
    df <- data.frame(t = time.series, c = covar, id = seq_along(covar))
    dontsample <- c(
      by(df, df$c, function(x) x$id[which.min(x$t)]),
      by(df, df$c, function(x) x$id[which.max(x$t)])
    )
  }
  else{
    dontsample <- c(which.min(time.series), which.max(time.series))
  }
  
  # find the nc for plsr interpol with scanned dfs.
  ncomps <- sapply(dfs, function(df) {
    dat$Y = I(splines::ns(y, df = df))
    cv <- suppressWarnings(pls::RMSEP(
      pls::plsr(m_formula, data = dat,
                scale = F, validation = "CV"))$val["CV", , ])
    nc <- which.min(colMeans(cv)) - 1
    if (nc == 0) 
      nc <- which.min(colMeans(cv)[-1])
    return(nc)
  })
  names(ncomps) <- paste0("df.", dfs)
  
  # randomly split training/validation sets 
  sels <- lapply(1:cv.n, function(i) {
    sample(seq_along(time.series)[- dontsample], n - n_samp, replace = F)
  })
  
  # Build splines for scanned dfs
  Ys <- lapply(dfs, function(df) {
    splines::ns(y, df = df)
  })
  
  # build cluster for parallel computing & export needed variables
  cl <- parallel::makeCluster(nb.cores, type = ifelse(.Platform$OS.type == "windows", "PSOCK", "FORK"))
  parallel::clusterExport(cl, varlist = c('m_formula', 'Ys', 'dat', 
                                          'sels', 'ncomps', 'dfs'),
                          envir = environment())
  
  cv.errors <- do.call("rbind", 
                       parallel::parLapply(cl, 1:cv.n, 
                                           function(i) {
                                             sel <- sels[[i]]
                                             
                                             errs <- sapply(dfs, function(df) {
                                               dat$Y <- I(Ys[[which(dfs == df)]])
                                               m <- pls::plsr(m_formula, data = dat[-sel,], scale = F)
                                               xv <- stats::predict(m, newdata = dat[sel,],
                                                                    comps = 1:ncomps[which(dfs == df)])
                                               err <- err.func(dat$X[sel, ], xv)
                                               return(err)
                                             })
                                             return(errs)
                                           }))
  colnames(cv.errors) <- paste0("df.", dfs)
  
  # stop cluster and free space 
  parallel::stopCluster(cl)
  gc(verbose = F)
  
  # format results
  res <- list(cv.errors = cv.errors, dfs = dfs, cv.n = cv.n, 
              cv.s = cv.s, plsr.nc = ncomps)
  if (isTRUE(ica.use)) {
    res$ica.nc <- ica.nc
  }
  
  class(res) <- "dfCV"
  return(res)
}



#' CV error function
#' 
#' Computes sum((xt-xe)^2)
#' 
#' @param xt,xe vectors of validation data and predictions respectively.
#' 
#' @return The sum squared prediction error.
ef <- function(xt, xe){
  return(sum((xt-xe)^2))
}






#' dfCV object summary
#' 
#' Prints a summary of the `dfCV` object
#' 
#' @param object a `dfCV` object, as returned by \code{\link{df_CV}}.
#' @param digits the number of digits passed on to \code{\link{round}}
#' @param ... ignored (needed to match the S3 standard)
#' 
#' @return invisibly returns the string output
#' 
#' @export
#' 
#' @eval interpol_example()
#' 
#' @importFrom stats median
summary.dfCV <- function(object, digits = 3, ...){

  out_str <- paste0('\nCross-Validation on spline df :\n',
                    '\nPLSR models fitted on ', 
                    ifelse(is.null(object$ica.nc), 
                           'whole gene expression matrix.',
                           paste0('ICA sample loadings (',object$ica.nc,' comps).')),
                    '\nPerformed ', object$cv.n, ' CV repeats with a ', 
                    ifelse(object$cv.s<1, paste0(round(object$cv.s*100),'/', round(100-object$cv.s*100), ' training/validation ratio.'),
                           paste0(object$cv.s, ' sample training set.')),
                    '\nScanned ', length(object$dfs),' df values  [', min(object$dfs), ' , ', max(object$dfs), ']',
                    '\n')
  cat(out_str)
  cat('\nMedian CVError per df :\n')
  print(apply(object$cv.errors, 2, stats::median), digits = digits)

  invisible(out_str)
}


#' dfCV object print
#' 
#' Prints a summary of the `dfCV` object  (same as summary)
#' 
#' @param x a `dfCV` object, as returned by \code{\link{df_CV}}.
#' @param digits the number of digits passed on to \code{\link{round}}
#' @param ... ignored (needed to match the S3 standard)
#' 
#' @return invisibly returns the string output
#' 
#' @export
#' 
#' @eval interpol_example()
#' 
print.dfCV <- function(x, digits = 3, ...){
  invisible(summary.dfCV(x, digits = digits))
}
