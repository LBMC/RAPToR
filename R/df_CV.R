#' Find the optimal spline df through cross-valdation
#' 
#' Some text about how cross-validation is performed here. And using ICA components.
#' 
#' @param X gene expression matrix of reference time series, genes as rows, (ordered) individuals as columns.
#' @param time.series timepoints of the reference (X).
#' @param dfs vector of spline df parameters to scan (passed on to \code{\link[splines]{ns}()} function).
#' @param cv.n number of cross-validation repeats.
#' @param cv.s ratio of samples to use for training set. If cv.s > 1, then cv.s samples are used for the training set.
#' @param err.func the error function to use to compute the CV error. Defaults to sum square of differences (the \code{\link{ef}()} function).
#' @param ica.use boolean ; if TRUE, sample loadings of ICA components are used for CV rather than the whole gene expression matrix.
#' @param ica.nc number of components to use for the ica. Defaults to 16 or ncol(X) if there are less than 16 samples.
#' @param nb.cores number of cores to use for parallelization.
#' 
#' @return a '\code{dfCV}' object, which is a list of ...
#' There are plot, print and summary methods for this object.
#' 
#' @export
#' 
#' @examples 
#' data(Cel_embryo)
#' 
#' dfCVembryo <- df_CV(Cel_embryo$X, Cel_embryo$time.series, 
#'                     dfs = 3:17, cv.n = 100, cv.s = 0.8, 
#'                     nb.cores = 4)
#' 
#' \donttest{
#' 
#' }
#' 
#' @importFrom parallel parLapply stopCluster makeForkCluster
#' @importFrom stats quantile dnorm predict
#' @importFrom splines ns
#' @importFrom ica icafast
#' @import pls 
#' 
#' 
df_CV <- function(X, time.series, 
                  dfs = 1:(ncol(X)-1), 
                  cv.n = 50, cv.s = 0.8, 
                  err.func = ef,
                  ica.use = TRUE, ica.nc = min(ncol(X), 16),
                  nb.cores = 2){

  y <- time.series
  n <- length(y)
  
  if(isTRUE(ica.use)){
    # Use ICA sample loadings for CV
    X <- ica::icafast(X, nc = ica.nc)$M 
  }
  else{
    X <- t(X)
  }
  
  if(cv.s<1.0)
    n_samp <- round(cv.s*n, 0)
  else
    n_samp <- cv.s
  
  # Determine the correct ncomps to use in plsr for scanned df values
  ncomps <- sapply(dfs, function(df){
    cv <- suppressWarnings(RMSEP(plsr(X~ns(y, df = df), scale=F, validation='CV'))$val['CV',,])
    nc <- which.min(colMeans(cv))-1
    if(nc==0)
      nc <- which.min(colMeans(cv)[-1])
    nc
  })
  names(ncomps) <- paste0('df.', dfs)
  
  # Make cv.n random subsets
  sels <- lapply(1:cv.n, function(i){
    sample(2:(n-1), n-n_samp, replace = F)
  })
  
  # Build splines for scanned dfs
  Ys <- lapply(dfs, function(df){
    ns(y, df=df)
  })
  
  # Set up cluster for parallelization
  cl <- parallel::makeCluster(nb.cores, 
                              type = ifelse(.Platform$OS.type=="windows", 
                                            "PSOCK", "FORK"))
  
  
  # Perform cross-validation on df parameter
  cv.errors <- do.call('rbind', 
                       parallel::parLapply(cl, 1:cv.n, function(i){
                         sel <- sels[[i]]
                         
                         errs <- sapply(dfs, function(df){
                           Y <- Ys[[which(dfs==df)]]
                           yi <- Y[-sel,, drop=F] # used for model
                           yv <- Y[sel,, drop=F] # used for cv
                           
                           m <- plsr(X[-sel,]~yi, scale=F)
                           
                           xv <- predict(m, newdata=yv, comps=1:ncomps[which(dfs==df)])
                           err <- err.func(X[sel,], xv)
                           return(err)
                           
                         })
                         return(errs)
                       }))
  
  
  # Stop cluster and free space
  parallel::stopCluster(cl)
  gc(verbose = F)
  
  res <- list(
    cv.errors = cv.errors,
    dfs = dfs,
    cv.n = cv.n,
    cv.s = cv.s,
    plsr.nc = ncomps)
  
  if(isTRUE(ica.use))
    res$ica.nc <- ica.nc
  
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

