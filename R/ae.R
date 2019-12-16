#' Age Estimate
#' 
#' This function estimates the developmental age of sample individuals based on 
#' correlation with given reference data. 
#' 
#' The implemented bootstrap procedure re-estimates the age on random gene subsets
#' of fixed size to evaluate the robustness of the estimate, given in the form of
#' the Median Absolute Deviation of the bootstrap age estimates to the global estimate. 
#' 
#' Using interpolated reference data gives more precise results.
#' 
#' A prior can be given to help with the estimate, in which case the peaks 
#' of the correlation profiles will be scored according to a gaussian
#' of the specified parameters.
#' 
#' @param samp the sample matrix, genes as rows, individuals as columns
#' @param refdata the reference time series matrix, same format as \code{samp}
#' @param ref.time_series the reference time series (\emph{e.g.} \code{interpol$time.series} if using interpolated reference data)
#' @param nb.cores the number of cores on which to parallelize the process, defaults to 2.
#' @param cor.method correlation method argument passed on to \code{\link{cor.gene_expr}}. Note that the Spearman coefficient performs much better than Pearson (while a bit slower).
#' @param bootstrap.n the number of bootstrap steps. Should ideally be >5
#' @param bootstrap.set_size the size of the random sub genesets for the bootstrap, defaults to ngenes/3 (ngenes being the number of \emph{overlapping} genes between sample and reference).
#' @param prior a vector with an approximate development time of the samples, must be in the same units than \code{ref.time_series}. Vector is recycled if its length is smaller than the number of samples
#' @param prior.params the std. deviation of the prior scoring distribution. \emph{Note that setting this value too low can cause a significant bias in the age estimation.}
#' @param verbose boolean ; if TRUE, displays messages of the various steps of the method.
#' 
#' @return an `ae` object, which is a list of the age estimates, the correlation matrix between sample and reference, 
#' the reference time series as well as the bootstrap correlation matrices and age estimates.
#' There are `plot`, `print` and `summary` methods for this object.
#' 
#' @export
#' 
#' @eval ae_example()
#' 
#' 
#' @importFrom parallel parApply parSapply parLapply stopCluster makeForkCluster
#' @importFrom stats quantile dnorm
#' @importFrom pryr standardise_call
#' 
ae <- function(samp, refdata, ref.time_series,
               cor.method = "spearman", nb.cores = 2,
               bootstrap.n = 30, bootstrap.set_size = NULL,
               prior = NULL, prior.params = NULL,
               verbose = T)
{
  if(length(ref.time_series)!=ncol(refdata)){
    stop("Reference data and time series don't match")
  }
  ref.time_series <- as.numeric(ref.time_series)
  
  
  ncs <- ncol(samp)
  dup <- FALSE
  if(ncs<=1){
    # if there is only one sample, double it to avoid 
    # dimension drop problems with R
    samp <- cbind(samp,dup=samp)
    ncs <- 2
    dup <- TRUE
  }
  
  if(bootstrap.n<=5){
    warning("bootstrap.n should ideally be > 5")
  }
  
  # get matching geneset
  if(!identical(rownames(refdata),rownames(samp))){
    overlap <- format_to_ref(samp, refdata, verbose = verbose)
    samp <- overlap$samp
    refdata <- overlap$refdata
    rm(overlap); gc(verbose = F)
  }
  if(is.null(bootstrap.set_size)){
    # default set size
    bootstrap.set_size <- round(nrow(samp)/3)
    if(verbose){
      message(paste("Bootstrap set size is", bootstrap.set_size))
    }
  }
  if(bootstrap.set_size>nrow(samp)){
    stop("bootstrap.set_size must be smaller than the overlapping geneset between sample and reference")
  }
  
  
  if(!is.null(prior.params)){
    if(is.null(prior)){
      stop("prior value must be specified if prior.params are defined")
    }
    if(length(prior)!=ncs){
      prior <- rep(prior, ncs)
    }
    prior.params <- rep(prior.params, ncs)
    
    if(any(prior>max(ref.time_series)|prior<min(ref.time_series))){
      stop("Some priors are outside reference time series' range")
    }
    # build prior densities (normed)
    range01 <- function(x){(x-min(x))/(max(x)-min(x))}
    priors <- lapply(1:ncs, function(i){
      range01(stats::dnorm(ref.time_series, mean = prior[i], sd = prior.params[i]))
    })
    
    # function to get the cor peak with prior
    get.cor_peak.prior <- function(cor.s, i){
      cor.maxs.i <- unique(c(which(diff(sign(diff(cor.s))) == -2) + 1, which.max(cor.s)))
      cor.maxs <- cor.s[cor.maxs.i]
      cor.max <- max(cor.maxs)
      cor.min <- min(cor.s)
      cor.maxs.times <- ref.time_series[cor.maxs.i]
      cor.maxs.scores <- (priors[[i]][cor.maxs.i] + ((cor.maxs-cor.min)/(cor.max-cor.min)))/2
      
      chosen <- which.max(cor.maxs.scores)
      
      return(cbind(time = cor.maxs.times[chosen], cor.score = cor.maxs[chosen]))
    }
    
  }
  else {
    
    # function to get the cor peak
    get.cor_peak <- function(cor.s){
      
      cor.max.i <- which.max(cor.s)
      cor.max <- cor.s[cor.max.i]
      cor.max.time <- ref.time_series[cor.max.i]
      
      return(cbind(time = cor.max.time, cor.score = cor.max))
    }
  }

  # build clusters for parallel comp. (detecting OS for cluster type)
  cl <- parallel::makeCluster(nb.cores, 
                              type = ifelse(.Platform$OS.type=="windows", 
                                            "PSOCK", "FORK"))
  
  samp.seq <- 1:ncol(samp)
  boot.seq <- 1:bootstrap.n
  
  if(verbose){
    message("Performing age estimation...")
  }
  
  # do age estimation on whole dataset
  cors <- cor.gene_expr(samp, refdata, cor.method = cor.method)
  if(is.null(prior)){
    #no prior
    age.estimates <- parallel::parApply(cl, cors, 2, get.cor_peak)
  } else {
    #with prior
    age.estimates <- parallel::parSapply(cl, samp.seq, function(i){
      get.cor_peak.prior(cors[,i], i)
    })
  }
  
  if(verbose){
    message("Bootstrapping...")
  }
  ### Bootstrap
  # generate gene subsets
  if(verbose){
    message("\tBuilding gene subsets...")
  }
  totalset <- 1:nrow(samp)
  gene_subsets <- lapply(boot.seq, function(i){
    sample(totalset, size = bootstrap.set_size, replace = F)
  })
  
  # get bootstrap correlations 
  if(verbose){
    message("\tComputing correlations...")
  }
  boot.cors <- simplify2array(parallel::parLapply(cl, boot.seq, function(j){
    cor.gene_expr(samp[gene_subsets[[j]], ,drop=F], refdata[gene_subsets[[j]], ,drop=F],
                  cor.method=cor.method)
  }))
  
  # get bootstrap age estimates
  if(verbose){
    message("\tPerforming age estimation...")
  }
  if(is.null(prior)){
    #no prior
    boots <- simplify2array(parallel::parLapply(cl, boot.seq,function(j){
      bcors <- boot.cors[,,j]
      age.estimate <- apply(bcors, 2, get.cor_peak)
    }))
  } else {
    #with prior
    boots <- simplify2array(parallel::parLapply(cl, boot.seq,function(j){
      bcors <- boot.cors[,,j]
      age.estimate <- sapply(samp.seq, function(i){
        get.cor_peak.prior(bcors[,i], i)
      })
      age.estimate
    }))
  }
  
  
  if(verbose){
    message("Computing summary statistics...")
  }
  # compute MAD Confidence Interval from boostrap
  resolution <- mean(diff(ref.time_series))/2
  conf.inter <- t(parallel::parSapply(cl, 1:dim(boots)[2], function(i) {
    m <- stats::mad(boots[1, i, ], center=age.estimates[1, i], na.rm = T) + resolution
    return(age.estimates[1, i]+c(-m, m))
  }))
  
  # get IC95 and median on the bootstrap correlation curves
  bc95 <- parallel::parApply(cl, boot.cors, c(1,2), quantile, probs=c(0.025, 0.5, 0.975), na.rm = T)
  
  # check for edge-of-reference estimates
  qref <- stats::quantile(ref.time_series, probs=c(.05,.95))
  if(any(  (age.estimates[1,]<qref[1])|(age.estimates[1,]>qref[2])  ))
    warning("Some estimates come near the edges of the reference.\nIf possible, stage those on a different reference for confirmation.")
  
  # data formatting
  colnames(conf.inter) <- c('lb', 'ub')
  age.estimates <- cbind(age.estimate=age.estimates[1,], conf.inter,
                         cor.score=age.estimates[2,])
  rownames(age.estimates) <- colnames(samp)
  
  
  # stop cluster and free space
  parallel::stopCluster(cl)
  rm(boot.cors, samp, refdata)
  if(is.null(prior))
    rm(get.cor_peak)
  else
    rm(get.cor_peak.prior)
  gc(verbose = F)
  
  res <- list(age.estimates = age.estimates, 
              ref.time_series = ref.time_series, 
              cors = cors,  cors.95 = bc95,
              boots = boots,
              prior = prior)
  
  if(dup){
    # if single sample was doubled, get only one result
    res <- list(age.estimates = age.estimates[1,,drop=F], 
                ref.time_series = ref.time_series, 
                cors = cors[,1,drop=F],  cors.95 = bc95[,,1,drop=F],
                boots = boots[,1,,drop=F],
                prior = prior[1])
  }
  
  res$call <- deparse(pryr::standardise_call(sys.call()))
  
  class(res) <- "ae"
  return(res)
  
}











#' Print an ae object
#' 
#' Prints the \code{age.estimates} dataframe of an \code{ae} object
#' 
#' @param x an \code{ae} object, as returned by \code{\link{estimate.worm_age}}.
#' @param digits the number of digits passed on to \code{\link{round}}
#' @param ... arguments passed on to \code{\link{print}}
#' 
#' @export
#' 
#' @eval ae_example()
#' 
print.ae <- function(x, digits=3, ...){
  cat("ae object\n\n")
  cat("Call : \n")
  sapply(x$call, function(s){cat(paste0('\t',s,'\n'))})
  cat('\n')
  print(round(x$age.estimates, digits = digits), ...)
}



#' Print an ae object summary
#' 
#' Prints a summary of the \code{age.estimates} dataframe of an \code{ae} object
#' 
#' @param object an \code{ae} object, as returned by \code{\link{estimate.worm_age}}.
#' @param digits the number of digits passed on to \code{\link{round}}
#' @param ... ignored (needed to match the S3 standard)
#' 
#' @return a list holding a data frame of ordered age estimates and confidence intervals, the span of estimates and the range.
#' 
#' @export
#' 
#' @eval ae_example()
#' 
summary.ae <- function(object, digits=3, ...){
  # rank the samples by age
  ord <- order(object$age.estimates[,1]) 
  ae_tab <-  cbind(round(object$age.estimates[ord,1:3], digits = digits)
                   # ,w=sapply(object$age.estimates[ord,'IC.imbalance'],
                   #          function(im){ifelse(im>5, '*','')})
                   )
  # colnames(ae_tab)[4] <- ' '
  
  bar <- paste0(rep('-', 50+digits*3), collapse = '')
  
  # Display timespan of samples
  mn <- as.numeric(min(object$age.estimates[,1], na.rm = T))
  mx <- as.numeric(max(object$age.estimates[,1], na.rm = T))
  
  cat(paste0('\nSpan of samples : ',
             round(mx-mn, digits = digits),
             '\nRange of samples :  [ ', 
             round(mn, digits = digits),' , ',round(mx, digits = digits),
             ' ]',
             '\n', bar, '\n'))
  # Print the sample table
  print(ae_tab, quote = F, right = T)
  cat(bar)
  # if(any(ae_tab[,4]=='*'))
  #   cat('\n\t* : Warning, estimate jumps around on bootstrap')
  # 
  invisible(list(age.estimates=object$age.estimates[ord,1:3], span=c(mn,mx), range=mx-mn))
}


#' Age Estimate - **DEPRECATED**
#' 
#' **DEPRECATED**. Please use \code{\link{ae}} instead.
#' 
#' @param samp passed on to ae
#' @param refdata passed on to ae
#' @param ref.time_series passed on to ae
#' @param ... passed on to ae
#' 
#' @export
#'
estimate.worm_age <- function(samp, refdata, ref.time_series, ...)
{
  warning("DEPRECATED. Please use ae() instead.")
  return(ae(samp, refdata, ref.time_series, ...))
}
