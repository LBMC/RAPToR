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
#' @param est.time a vector with the approximate development time of the samples, must be in the same units than \code{ref.time_series}. Vector is recycled if its length is smaller than the number of samples
#' @param time.sd the std. deviation of the gaussian scoring distribution. \emph{Note that setting this value too low can cause a significant bias in the age estimation.}
#' @param cor.method correlation method argument passed on to \code{\link{cor.gene_expr}}. Note that the spearman coefficient performs much better than the others.
#' @param bootstrap.n the number of re-estimates done by the bootstrap. If set to 0, the 95\% interval is computed from the reference time series' resolution
#' @param bootstrap.time_window the width of the window in which bootstrap re-estimates occur.
#' @param cors the correlation matrix between sample and reference data. \bold{This must be exactly} \code{cor.gene_expr(samp, refdata, cor.method = cor.method)}
#' 
#' @return an '\code{ae}' object, which is a list of the correlation matrix between sample and reference, the age estimates, the initial time estimates and the reference time series.
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
                              time.sd=5, cor.method="spearman",
                              bootstrap.n = 50, bootstrap.time_window = 2,
                              cors = NULL)
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
  ref.time_series <- as.numeric(ref.time_series)
  
  if(!identical(rownames(refdata),rownames(samp))){
    overlap <- format_to_ref(samp, refdata)
    samp <- overlap$samp
    refdata <- overlap$refdata
    rm(overlap); gc()
  }
  
  
  if(is.null(cors)){
    # compute correlations
    cors <- cor.gene_expr(samp, refdata, cor.method = cor.method)
  }
  
  if(bootstrap.n>0){
    b.mod <- c(0, runif(bootstrap.n, -.5*bootstrap.time_window,
                        .5*bootstrap.time_window))
  }
  else{
    b.mod <- 0
  }
  boots <- sapply(1:length(b.mod), function(j){
    b.est_time <- est.time + b.mod[j]
    
    ref.gauss <- lapply(b.est_time, function(b.et){
      dnorm(ref.time_series, mean = b.et, sd = time.sd)
    })
    m.gauss <- lapply(ref.gauss, max)
    
    
    age.estimates <- lapply(1:ncol(samp), function(i) {
      
      cor.maxs.i <- which(diff(sign(diff(cors[, i]))) == -2) + 1
      if(length(cor.maxs.i)==0){
        # No maxima found 
        age.estimate <- cbind(time = NA, cor.score = NA,
                              proba.score = NA)
      }
      else{
        cor.maxs <- cors[cor.maxs.i, i]
        cor.max <- max(cor.maxs)
        cor.min <- min(cors[,i])
        cor.maxs.times <- ref.time_series[cor.maxs.i]
        cor.maxs.scores <- round(((ref.gauss[[i]][cor.maxs.i]/m.gauss[[i]]) + 
                                    ((cor.maxs-cor.min)/(cor.max-cor.min)))/2, 4)
        age.estimate <- cbind(time = cor.maxs.times, cor.score = cor.maxs, 
                              proba.score = cor.maxs.scores)
        age.estimate <- age.estimate[order(age.estimate[, "proba.score"], 
                                           decreasing = T), , drop=F]
      }
      return(age.estimate)
    })
    # get best estimate
    age.estimates <- do.call('rbind',lapply(age.estimates, 
                                            function(a.e) {
                                              return(a.e[1, ,drop=F])
                                            }))
    
    return(age.estimates)
  }, simplify = 'array')
  
  
  # get average & IC95 over boostrap
  age.estimates <- sapply(1:dim(boots)[1], function(i){apply(boots[i,,, drop=F],c(1,2),mean)})
  
  resolution <- mean(diff(ref.time_series))/2
  age.est95 <- t(sapply(1:dim(boots)[1], function(i){
    quantile(boots[i,1,], probs=c(0.025,0.975), na.rm = T)+c(-1,1)*resolution
  }))
  
  
  age.estimates <- cbind(age.estimate=age.estimates[1,], age.est95, 
                         cor.score=age.estimates[2,], proba.score=age.estimates[3,])
  rownames(age.estimates) <- colnames(samp)
  
  res <- list(cors = cors, age.estimates = age.estimates, ref.time_series = ref.time_series, 
              init.est.times = est.time)
  
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
#' age.est <- estimate.worm_age(samp, oud_ref$X, oud_ref$time.series, 26)
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
                 xlim=range(c(err.inf, err.sup, age_est$init.est.times)),
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
    inis <- age_est$init.est.times[o]
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
#' age.est <- estimate.worm_age(samp, oud_ref$X, oud_ref$time.series, verbose=F)
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
    
    
    if (show.init_estimate&!is.null(age.est$init.est.times)) {
      # show initial estimate 
      init.est <- age.est$init.est.times[i]
      points(init.est, min(age.est$cors[, i]), pch = "|", 
             col = in.col, cex = bar.size)
      text(init.est, min(age.est$cors[, i]), pos = 3, offset = 1, 
           labels = paste(round(init.est, 2), "\n(initial estimate)", 
                          sep = ""))
    }
  })
}
