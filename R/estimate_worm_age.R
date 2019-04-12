#' Estimate the developmental age of individuals
#' 
#' This function estimates the developmental age of individuals based on given
#' reference data as well as an approximate time estimate of the age.
#' This is achieved using a two-step process.
#' First, simple correlation between the gene expression profiles of
#' the sample and reference data are computed using \code{\link{cor.gene_expr}}; 
#' this gives correlation profiles according to reference time per individual. 
#' Then, the maxima of the correlation profiles are evaluated and scored according to 
#' a gaussian distribution around the time estimate given as input and their correlation score. 
#' The maxima scores range from 0 to 1, 1 being the case where the highest correlation peak is 
#' exactly at the given approximate time. 
#' The implemented bootstrap procedure re-estimates the age from the given times with a 
#' random uniform noise and returns the average best time, as well as the 95% quantile 
#' interval of the values.
#' 
#' Do note that using interpolated reference data (from \code{\link{interpol_refdata}})
#' gives the best results.
#' 
#' @param samp the sample matrix, gene as rows, individuals as columns
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
#' @return an '\code{ae}' object, which is a list of the correlation matrix between sample and reference, the age estimates, the reference time series as well as the bootstrap correlation matrices and age estimates.
#' 
#' @export
#' 
#' @examples 
#' data(oud_ref)
#' 
#' samp <- oud_ref$X[,13:15]
#' age.est <- estimate.worm_age(samp, oud_ref$X, oud_ref$time.series)
#' age.est$age.estimates
#' \donttest{
#' plot(age.est)
#' }
#' 
#' 
estimate.worm_age <- function(samp, refdata, ref.time_series,
                              cor.method="spearman", nb.cores=2,
                              bootstrap.n = 30, bootstrap.set_size = NULL,
                              prior=NULL, prior.params=NULL,
                              verbose=T)
{
  requireNamespace('parallel', quietly = T)
  if(length(ref.time_series)!=ncol(refdata)){
    stop("Reference data and time series don't match")
  }
  ref.time_series <- as.numeric(ref.time_series)
  
  ncs <- ncol(samp)
  dup <- FALSE
  if(ncs<=1){
    # if there is only one sample, double it to avoid 
    # dimension problems with R
    samp <- cbind(samp,dup=samp)
    ncs <- 2
    dup <- TRUE
  }
  
  if(bootstrap.n<=5){
    message("Note : bootstrap.n should ideally be > 5")
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
    message(paste("Bootstrap set size is", bootstrap.set_size))
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
      stop("Some estimated times are outside reference time series' range")
    }
    # build prior densities (normed)
    range01 <- function(x){(x-min(x))/(max(x)-min(x))}
    priors <- lapply(1:ncs, function(i){
      range01(dnorm(ref.time_series, mean = prior[i], sd = prior.params[i]))
    })
    
    # function to get the cor peak with prior
    get.cor_peak <- function(cor.s, i){
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
  
  # build clusters for parallel comp.
  cl <- parallel::makeForkCluster(nb.cores)
  
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
      get.cor_peak(cors[,i], i)
    })
  }
  
  print(age.estimates)
  
  if(verbose){
    message("Bootstrapping...")
  }
  ### Bootstrap
  # generate gene subsets
  if(verbose){
    message("\tBuilding gene subsets...")
  }
  totalset <- 1:nrow(samp)
  gene_subsets <- parallel::parLapply(cl, boot.seq, function(i){
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
    boots <- simplify2array(parallel::parLapply(cl, boot.seq,function(j){
      bcors <- boot.cors[,,j]
      age.estimate <- sapply(samp.seq, function(i){
        get.cor_peak(bcors[,i], i)
      })
      age.estimate
    }))
  }
  
  
  if(verbose){
    message("Computing summary statistics...")
  }
  # get average & IC95 over boostrap
  resolution <- mean(diff(ref.time_series))/2
  age.est95 <- t(parSapply(cl, 1:dim(boots)[2], function(i){
    quantile(boots[1,i,], probs=c(0.025,0.975), na.rm = T)+c(-1,1)*resolution
  }))
  
  # get IC95 and median on the bootstrap correlation curves
  bc95 <- parallel::parApply(cl, boot.cors, c(1,2), quantile, probs=c(0.025, 0.5, 0.975), na.rm = T)
  
  
  # data formatting
  age.estimates <- cbind(age.estimate=age.estimates[1,], age.est95,
                         cor.score=age.estimates[2,])
  rownames(age.estimates) <- colnames(samp)
  
  
  # stop cluster and free space
  parallel::stopCluster(cl)
  rm(boot.cors, samp, refdata, get.cor_peak) ;  gc(verbose = F)
  
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
  
  class(res) <- "ae"
  return(res)
  
}





#' Plot an ae object
#' 
#' Plots the correlation score curves from samples against a reference series \code{\link{cor.gene_expr}} 
#' 
#' @param age.est an \code{ae} object, as returned by \code{\link{estimate.worm_age}} 
#' @param errbar.width the width of the error bars
#' @param show.init_estimate logical ; if TRUE, shows the initial time estimate(s) on the plot
#' @param col.i the color of the initial estimate marker.
#' @param groups a factor with sample categories, as passed on to \code{\link{dotchart}}.
#' @param pch the pch parameter passed on to \code{\link{dotchart}}.
#' @param cex sizing parameter applied to various elements of the plot.
#' @param ... additional arguments passed on to \code{\link{dotchart}}
#' 
#' @export
#' 
#' @examples
#' data(oud_ref)
#' 
#' samp <- oud_ref$X[,13:15]
#' age.est <- estimate.worm_age(samp, oud_ref$X, oud_ref$time.series)
#' \donttest{
#' plot(age.est)
#' }
#' 
plot.ae <- function(age_est, errbar.width=0.1, 
                    show.init_estimate=F, col.i=1,
                    groups=NULL, 
                    pch=16, cex=1, 
                    xlab="Estimated ages", ...)
{
  err.inf <- age_est$age.estimates[,2]
  err.sup <- age_est$age.estimates[,3]
  n <- nrow(age_est$age.estimates)
  
  dc <- dotchart(age_est$age.estimates[,1], labels = rownames(age_est$age.estimates),
                 xlab=xlab, groups = groups,
                 xlim=range(c(err.inf, err.sup, age_est$prior)),
                 pch=pch, cex=cex,
                 ...)
  
  # Adjusting Y positions of error bars to dotchart layout
  y <- 1L:n
  o <- y
  
  if(!is.null(groups)){
    o <- sort.list(as.numeric(groups), decreasing = TRUE)
    err.inf <- err.inf[o]
    err.sup <- err.sup[o]
    groups <- groups[o]
    offset <- cumsum(c(0, diff(as.numeric(groups)) != 0))
    y <- 1L:n + 2 * offset
  }
  # plot error bars
  arrows(err.sup, y,
         err.inf, y,
         angle=90, code=3, length=errbar.width)
  
  # adding initial estimate to plot
  if(show.init_estimate){
    inis <- age_est$prior[o]
    col.i <- rep(col.i, n)
    col.i <- col.i[o]
    points(inis, y, lwd=2, cex=cex*1.1, col=col.i)
    #text(inis[n], y[n], labels = "initial estimate", pos = 3, col=col.i[n])
    legend('bottomleft', legend = ' Initial estimate', col = col.i[n], inset = .02,
           pt.lwd=2, pch=1, bty = 'n', text.col = col.i[n])
  }
}




#' Plot an ae object
#' 
#' Plots the correlation score curves from samples against the reference series.
#' 
#' @param age.est an \code{ae} object, as returned by \code{\link{estimate.worm_age}} 
#' @param subset an index vector of the samples to plot (defaults to all)
#' @param show.init_estimate logical ; if TRUE, shows the initial time estimate on the plot
#' @param c.lwd line width for the correlation score curve
#' @param bar.size size of the estimate 95IC bars
#' @param mx.col color of the age estimate bars
#' @param in.col color of the initial estimate bar
#' @param ... additional arguments passed on to \code{\link{plot}}
#' 
#' @export
#' 
#' @examples
#' data(oud_ref)
#' 
#' samp <- oud_ref$X[,13:15]
#' \donttest{
#' age.est <- estimate.worm_age(samp, oud_ref$X, oud_ref$time.series)
#' plot_cor.ae(age.est)
#' }
#' 
plot_cor.ae <- function (age.est, subset = 1:ncol(age.est$cors), 
                         show.init_estimate = F,
                         c.lwd = 2, bar.size = 1, 
                         mx.col = "firebrick", in.col = "royalblue", 
                         ...) 
{
  pb <- sapply(subset, function(i) {
    # set ylim values
    if(!is.null(age.est$cors.95))
      yl <- range(age.est$cors.95[,,i])
    else yl <- range(age.est$cors[,i])
    
    # plot corr. curve
    plot(age.est$ref.time_series, age.est$cors[,i], type = "l", 
         lwd = c.lwd, main = colnames(age.est$cors)[i], xlab = "reference time",
         ylim = yl*c(1,1.025), ylab = "corr.score", ...)
    
    # if bootstrap cor 95 IC was returned, plot cor curve 95 IC & median
    if(!is.null(age.est$cors.95))
      sapply(1:3, function(j){
        points(age.est$ref.time_series, age.est$cors.95[j,,i], 
               type = 'l', lwd=c(1,2,1)[j], lty=2)
      })
    
    # get values for current sample
    ae <- age.est$age.estimates[i, c("age.estimate", "cor.score", "2.5%", "97.5%")]
    
    # plot 95IC bars & band
    seg.h <- ((yl[2]-yl[1])/10)*bar.size
    xs <- c(ae[3], ae[4])
    y0s <- rep(ae[2]-seg.h, 2)
    y1s <- rep(ae[2]+seg.h, 2)
    segments(xs, y0s, y1 = y1s, lwd=2.5, col = mx.col)
    
    yp <- c(y0s[1]+seg.h/2, y1s[1]-seg.h/2)
    polygon(rep(xs, each=2), y=c(yp[1], yp[2], yp[2], yp[1]), 
            col=makeTransparent(mx.col, alpha = 150), border = NA)
    
    
    # add estimate as text
    text(ae[1],ae[2]-3*seg.h, labels = paste(round(ae[1], 2), sep = ""))
    
    
    if (show.init_estimate&!is.null(age.est$prior)) {
      # show initial estimate 
      init.est <- age.est$prior[i]
      points(init.est, min(age.est$cors[, i]), pch = "|", 
             col = in.col, cex = bar.size)
      text(init.est, min(age.est$cors[, i]), pos = 3, offset = 1, 
           labels = paste(round(init.est, 2), "\n(initial estimate)", 
                          sep = ""))
    }
  })
}

#' Make a color transparent
#' 
#' Makes any given color(s) transparent
#' 
#' @param color any color.
#' @param alpha the alpha channel value, (0:255) - from fully transparent to opaque
#' 
#' @export
#' 
#' @examples
#' \donttest{
#' plot(1:10, col=makeTransparent('firebrick', 120), pch=16, cex=2)
#' }
#' 
makeTransparent<-function(color, alpha=100)
{
  newColor<-col2rgb(color)
  apply(newColor, 2, function(curcoldata){
    rgb(red=curcoldata[1], green=curcoldata[2],
        blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}
