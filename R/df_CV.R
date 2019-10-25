#' Find the optimal spline df through cross-valdation
#' 
#' Some text about how crosss-validation is performed here.
#' 
#' @param X the sample matrix, gene as rows, individuals as columns
#' 
#' @return a '\code{dfCV}' object, which is a list of the age estimates, the correlation matrix between sample and reference, 
#' the reference time series as well as the bootstrap correlation matrices and age estimates.
#' There are plot, print and summary methods for this object.
#' 
#' @export
#' 
#' @examples 
#' data(Cel_embryo)
#' 
#' \donttest{
#' 
#' }
#' 
#' @importFrom parallel parApply parSapply parLapply stopCluster makeForkCluster
#' @importFrom stats quantile dnorm
#' @importFrom pryr standardise_call
#' 
df_CV <- function(X, y, dfs=1:length(y), n.cv=50, s.cv=0.8, err.func=ef, ncores=2, ...){
  require(parallel)
  require(splines)
  require(pls)
  
  n <- length(y)
  if(s.cv<1.0)
    n_samp <- round(s.cv*n, 0)
  else
    n_samp <- s.cv
  
  # determine the correct ncomps to use in plsr for scanned df values
  ncomps <- sapply(dfs, function(df){
    cv <- suppressWarnings(RMSEP(plsr(X~ns(y, df = df), scale=F, validation='CV'))$val['CV',,])
    nc <- which.min(colMeans(cv))-1
    if(nc==0)
      nc <- which.min(colMeans(cv)[-1])
    nc
  })
  names(ncomps) <- paste0('df.', dfs)
  
  # Make n.cv random subsets
  sels <- lapply(1:n.cv, function(i){
    sample(2:(n-1), n-n_samp, replace = F)
  })
  
  # Build splines for scanned dfs
  Ys <- lapply(dfs, function(df){
    ns(y, df=df)
  })
  
  # Perform cross-validation on df parameter
  cv_errors <- do.call('rbind', 
                       parallel::mclapply(1:n.cv, function(i){
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
                       }
                       , mc.cores = ncores
                       ))
  
  res <- list(dfs = dfs,
              cv_errors = cv_errors,
              n.cv = n.cv,
              s.cv = s.cv,
              plsr.nc = ncomps)
}

# CV error function between Xtrue and Xestimated (sum squared)
ef <- function(xt, xe){
  sum((xt-xe)^2)
}

