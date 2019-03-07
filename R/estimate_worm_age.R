#' Estimate the developmental age of individuals
#' 
#' This function estimates the developmental age of individuals based on given
#' reference data as well as an approximate time estimate of the age.
#' This is achieved using a two-step process.
#' First, simple correlation between the gene expression profiles of
#' the sample and reference data are computed using \code{\link{cor.gene_expr}}; 
#' this gives correlation profiles according to reference time per individual. 
#' Then, the maxima of the correlation profiles are evaluated and scored according to 
#' a gaussian distribution around the time estimate given as input. 
#' The maxima scores range from 0 to 1, 1 being the case where a correlation peak is 
#' exactly at the given approximate time. 
#' Do note that using interpolated reference data (from \code{\link{interpol_refdata}})
#' gives the best results.
#' 
#' @param samp the sample matrix, gene as rows, individuals as columns
#' @param refdata the reference time series matrix, same format as \code{samp}
#' @param ref.time_series the reference time series (\emph{e.g.} \code{colnames(refdata)} if using interpolated reference data)
#' @param est.time a vector with the approximate development time of the samples, must be in the same units than \code{ref.time_series}. The vector is recycled if its length is smaller than the number of samples
#' @param time.sd the std. deviation of the gaussian scoring distribution. \emph{Note that setting this value too low can cause a significant bias in the age estimation.}
#' @param cor.method correlation method argument passed on to \code{\link{cor.gene_expr}}
#' @param all.peaks logical; if TRUE, returns all correlation peaks (potential age estimates) and their respective scores for every individuals, as a list. If FALSE, only returns the best estimate for each individual, as a dataframe.
#' 
#' 
#' @return an '\code{ae}' object, which is a list of the correlation matrix between sample and reference, the age estimates (either as list of individuals or dataframe, depending on \code{all.peaks}), the initial time estimates and the reference time series.
#' 
#' @export
#' 
#' @examples 
#' data(oud_ref)
#' 
#' samp <- oud_ref$X[,13:15]
#' age.est <- estimate.worm_age(samp, oud_ref$X, oud_ref$time.series, 26)
#' age.est$age.estimates
#' \donttest{
#' }
#' 
#' 
estimate.worm_age <- function(samp, refdata, ref.time_series, est.time,
                              time.sd=5, cor.method="pearson", all.peaks=F)
{
  if(length(ref.time_series)!=ncol(refdata)){
    stop("Reference data and time series don't match")
  }
  if(length(est.time)!=ncol(samp)){
    est.time <- rep(est.time, ncol(samp))
  }
  if(any(est.time>max(ref.time_series)|est.time<min(ref.time_series))){
    stop("Some estimated times are outside reference time series' range")
  }
  ref.time_series <-  as.numeric(ref.time_series)
  
  # compute correlations
  cors <- cor.gene_expr(samp, refdata, cor.method = cor.method)
  
  # compute gaussian around estimated time
  ref.gauss <- lapply(est.time, 
                      function(et){
                        dnorm(ref.time_series, mean = et, sd = time.sd)
                      })
  m.gauss <- lapply(ref.gauss, max)
  
  age.estimates <- lapply(1:ncol(samp), function(i){
    # get correlation maxima (peaks) positions
    cor.maxs.i <- unique(c(which(diff(sign(diff(cors[,i])))==-2)+1, 
                           which.max(cors[,i])))
    cor.maxs <- cors[cor.maxs.i, i]
    cor.maxs.times <- ref.time_series[cor.maxs.i]
    
    # compute scores based on gaussian of reference time
    cor.maxs.scores <- round(ref.gauss[[i]][cor.maxs.i]/m.gauss[[i]], 4)
    
    age.estimate <- cbind(time=cor.maxs.times, 
                          cor.score=cor.maxs,
                          proba.score=cor.maxs.scores)
    
    # order by score & formatting
    age.estimate <- age.estimate[order(age.estimate[, "proba.score"], decreasing = T),]
    if(length(cor.maxs)<2){age.estimate <- as.matrix(t(age.estimate))}
    
    return(age.estimate)
  })
  names(age.estimates) <- colnames(samp)
  if(!all.peaks){
    # only keep best estimate
    age.estimates <- simplify2array(lapply(age.estimates, 
                                           function(a.e){return(a.e[1,])}))
  }
  res <- list(cors=cors, 
              age.estimates=age.estimates, 
              ref.time_series=ref.time_series,
              init.est.times=est.time)
  class(res) <- "ae"
  return(res)
  
}






#' Plot an ae object
#' 
#' Plots the correlation score curves from samples against a reference series \code{\link{cor.gene_expr}} 
#' as a heatmap.
#' 
#' @param age.est an \code{ae} object, as returned by \code{\link{estimate.worm_age}} 
#' @param subset an index vector of the samples to plot (defaults to all)
#' @param show.init_estimate logical ; if TRUE, shows the initial time estimate on the plot
#' @param show.other_estimates logical ; if TRUE and if \code{\link{estimate.worm_age}} was called with \code{all.peaks=TRUE}, displays the other maxima of the curve
#' @param c.lwd line width for the correlation score curve
#' @param bar.size cex of the maxima bars
#' @param mx.col color of the best age estimate bar
#' @param in.col color of the initial estimate bar
#' @param ot.col color of the other maxima bars
#' @param ... additional arguments passed on to \code{\link{plot}}
#' 
#' @export
#' 
#' @examples
#' data(oud_ref)
#' 
#' samp <- oud_ref$X[,13:15]
#' age.est <- estimate.worm_age(samp, oud_ref$X, oud_ref$time.series, 26)
#' \donttest{
#' plot(age.est)
#' }
#' 
plot.ae <- function(age.est, subset=1:ncol(age.est$cors),
                    show.init_estimate=F, show.other_estimates=T,
                    c.lwd=2, bar.size=2,
                    mx.col='firebrick', in.col='royalblue', ot.col='grey50',
                    ...){
  
  pb <- sapply(subset, function(i){
    # plot corr.score curve
    plot(age.est$ref.time_series, age.est$cors[,i], type = 'l', lwd=c.lwd,
         main=colnames(age.est$cors)[i], 
         xlab = 'reference time', ylab='corr.score', ...)
    
    if(is.list(age.est$age.estimates)){
      # estimate.worm_ages called with all.peaks=T
      
      # get age.est df  
      aes <- age.est$age.estimates[[colnames(age.est$cors)[i]]]
      
      # plot best estimate
      ae.max <- t(aes[1, c(1,2)])
      points(ae.max, pch='|', cex=bar.size, col=mx.col)
      text(ae.max, pos=1, 
           labels = paste(round(ae.max[1], 2), sep=''),
           offset = 1)
      
      if(show.other_estimates){
        # plot other estimates
        points(aes[-1,1:2], pch='|', cex=bar.size/1.5, col=ot.col)
      }
      
    }
    else{
      # get age estimation and plot
      ae <- t(age.est$age.estimates[c(1,2),i])
      points(ae, pch='|', cex=bar.size, col=mx.col)
      text(ae, pos=1, 
           labels = paste(round(ae[1], 2), sep=''),
           offset = 1)
    }
    
    if(show.init_estimate){
      # show initial estimate
      init.est <- age.est$init.est.times[i]
      points(init.est, min(age.est$cors[,i]), pch='|', col=in.col, cex=bar.size)
      text(init.est, min(age.est$cors[,i]), pos=3, offset = 1,
           labels = paste(round(init.est, 2), '\n(initial estimate)', sep=''))
    }
  })
}
