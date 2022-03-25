#' Plot an ae object
#' 
#' Plots age estimates of sample with bootstrap confidence interval as a dotchart.
#' 
#' @param x an `ae` object, as returned by \code{\link{ae}}.
#' @param groups a factor with sample categories, as passed on to \code{\link{dotchart}}.
#' @param subset an index vector of the samples to plot (defaults to all).
#' @param show.boot_estimates logical ; if TRUE, shows the individual bootstrapped estimates on the plot as swarms.
#' @param show.prior logical ; if TRUE, shows the input prior(s) on the plot.
#' @param col, color color parameter passed on to \code{\link{dotchart}}.
#' @param errbar.width the width of the error bars.
#' @param col.b,col.p the color of the bootstrapped estimates and prior respectively.
#' @param glob.above logical ; if TRUE, the global estimate is plotted above all else.
#' @param pch passed on to \code{\link{dotchart}}.
#' @param cex sizing parameter applied to various elements of the plot.
#' @param xlim horizontal range for the plot, see \code{\link[graphics]{plot.window}}, for example
#' @param xlab the x axis label, passed on to \code{\link{dotchart}}.
#' @param l.pos the position of the legend when show.prior is \code{TRUE}, passed on to \code{\link[graphics]{legend}}
#' @param ... additional arguments passed on to \code{\link{dotchart}}.
#' 
#' @export
#' 
#' @eval ae_example()
#' 
#' @importFrom graphics plot dotchart points arrows legend
#' @importFrom beeswarm swarmy
#' 
plot.ae <- function(x, groups=NULL, subset=NULL,
                    show.boot_estimates=F, show.prior=F, 
                    color = par("fg"), col = color,
                    col.b=2, col.p=1,
                    errbar.width=0.1, glob.above = F,
                    pch=16, cex=1, xlim=NULL,
                    xlab=NULL, 
                    l.pos='bottomright', ...)
{
  if(!is.null(subset)){
    # subset the data to plot
    x$age.estimates <- x$age.estimates[subset,,drop=F]
    x$prior <- x$prior[subset]
    x$boots <- x$boots[,subset, ,drop=F]
    if(!is.null(groups)){
      groups <- droplevels(groups[subset])
    }
  }
  err.inf <- x$age.estimates[,2]
  err.sup <- x$age.estimates[,3]
  n <- nrow(x$age.estimates)
  
  if(is.null(xlim)){
    xlim <- range(c(err.inf, err.sup, x$prior))
  }
  
  dc <- graphics::dotchart(x$age.estimates[,1], labels = rownames(x$age.estimates),
                           xlab=xlab, groups = groups,
                           xlim=xlim, color = color,
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
  color <- rep(color, n)[o]
  
  # plot error bars
  arrows(err.sup, y,
         err.inf, y,
         angle=90, code=3, length=errbar.width, 
         col = color)
  
  # adding individual bootstrap estimates as swarms
  if(show.boot_estimates){
    nboot <- dim(x$boots)[3]
    xs <- x$boots[1,o,]
    col.b <- rep(col.b, n)
    col.b <- col.b[o]
    invisible(
      sapply(1:n, function(i){
        
        yi <- rep(y[i], nboot)
        sw <- beeswarm::swarmy(xs[i,], yi, cex=.08*cex)
        graphics::points(sw, pch=16, cex=.3*cex, col=col.b[i])
      })
    )
  }
  
  # adding initial estimate to plot
  if(show.prior){
    inis <- x$prior[o]
    col.p <- rep(col.p, n)
    col.p <- col.p[o]
    graphics::points(inis, y, lwd=2, cex=cex*1.1, col=col.p)
    graphics::legend(l.pos, legend = ' Prior', col = col.p[n], inset = .02,
                     pt.lwd=2, pch=1, bty = 'n', text.col = col.p[n])
  }
  
  if(glob.above){
    graphics::points(x$age.estimates[o,1], y, cex=cex, pch=16, col = color, ...)
  }
}




#' Plot an ae object
#' 
#' Plots the correlation score curves from samples against the reference series.
#' 
#' @param ae_obj an \code{ae} object, as returned by \code{\link{ae}} 
#' @param subset an index vector of the samples to plot (defaults to all)
#' @param show.prior logical ; if TRUE, shows the input prior(s) on the plot.
#' @param c.lwd line width for the correlation score curve
#' @param bar.size size of the estimate confidence interval bars
#' @param mx.col color of the age estimate bars
#' @param col.p color of the prior bar
#' @param xlab,ylab the x and y axis labels, passed on to \code{\link{plot}}
#' @param ... additional arguments passed on to \code{\link{plot}}
#' 
#' @export
#' 
#' @eval ae_example()
#' 
#' @importFrom graphics plot points segments text polygon
#' 
plot_cor <- function (ae_obj, subset = 1:ncol(ae_obj$cors), 
                         show.prior = F,
                         c.lwd = 2, bar.size = 1, 
                         mx.col = "firebrick", col.p = "royalblue", 
                         xlab = NULL, ylab = NULL,
                         ...) 
{
  pb <- sapply(subset, function(i) {
    # set ylim values
    if(!is.null(ae_obj$cors.95))
      yl <- range(ae_obj$cors.95[,,i])
    else yl <- range(ae_obj$cors[,i])
    
    if(is.null(xlab){
      xlab <- paste0("Reference time, ", attr(ae_obj, "t.unit"))
    }  
    if(is.null(ylab)){
      ylab <- "Corr. score"
    }
    # plot corr. curve
    graphics::plot(ae_obj$ref.time_series, ae_obj$cors[,i], type = "l", 
                   lwd = c.lwd, main = colnames(ae_obj$cors)[i], xlab = xlab,
                   ylim = yl*c(1,1.025), ylab = ylab, ...)
    
    # if bootstrap cor 95 IC was returned, plot cor curve 95 IC & median
    if(!is.null(ae_obj$cors.95))
      sapply(1:3, function(j){
        graphics::points(ae_obj$ref.time_series, ae_obj$cors.95[j,,i], 
                         type = 'l', lwd=c(1,2,1)[j], lty=2)
      })
    
    # get values for current sample
    ae <- ae_obj$age.estimates[i, c("age.estimate", "cor.score", "lb", "ub")]
    
    # plot IC bars & band
    seg.h <- ((yl[2]-yl[1])/10)*bar.size
    xs <- c(ae[3], ae[4])
    y0s <- rep(ae[2]-seg.h, 2)
    y1s <- rep(ae[2]+seg.h, 2)
    graphics::segments(xs, y0s, y1 = y1s, lwd=2.5, col = mx.col)
    
    yp <- c(y0s[1]+seg.h/2, y1s[1]-seg.h/2)
    graphics::polygon(rep(xs, each=2), y=c(yp[1], yp[2], yp[2], yp[1]), 
                      col=makeTransparent(mx.col, alpha = 150), border = NA)
    
    
    # add estimate as text
    graphics::text(ae[1],ae[2]-3*seg.h, labels = paste(round(ae[1], 2), sep = ""))
    
    
    if (show.prior&!is.null(ae_obj$prior)) {
      # show initial estimate 
      init.est <- ae_obj$prior[i]
      graphics::points(init.est, min(ae_obj$cors[, i]), pch = "|", 
                       col = col.p, cex = bar.size)
      graphics::text(init.est, min(ae_obj$cors[, i]), pos = 3, offset = 1, 
                     labels = paste(round(init.est, 2), 
                                    "\n(prior)", sep = ""))
    }
  })
}

#' (DEPRECATED) Plot an ae object
#' 
#' This function has been renamed. Please use \link{\code{plot_cor()}} instead. 
#' 
#' @param age.est an \code{ae} object, as returned by \code{\link{ae}} 
#' @param subset an index vector of the samples to plot (defaults to all)
#' @param show.prior logical ; if TRUE, shows the input prior(s) on the plot.
#' @param c.lwd line width for the correlation score curve
#' @param bar.size size of the estimate confidence interval bars
#' @param mx.col color of the age estimate bars
#' @param col.p color of the prior bar
#' @param xlab the x axis label, passed on to \code{\link{plot}}
#' @param ... additional arguments passed on to \code{\link{plot}}
#' 
#' @export
#' 
#' @eval ae_example()
#' 
#' @importFrom graphics plot points segments text polygon
#' 
plot_cor.ae <- function (age.est, subset = 1:ncol(age.est$cors), 
                         show.prior = F,
                         c.lwd = 2, bar.size = 1, 
                         mx.col = "firebrick", col.p = "royalblue", 
                         xlab = NULL,
                         ...) 
{
  warning("This function is now deprecated. Please use plot_cor() instead.")
  RAPToR::plot_cor(ae_obj = age.est, subset = subset, 
           show.prior = show.prior,
           c.lwd = c.lwd, bar.size = bar.size, 
           mx.col = mx.col, col.p = col.p, 
           xlab = xlab, ...)
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
#' plot(1:255, col=makeTransparent('firebrick', 1:255), pch=16, cex=2)
#' }
#' 
#' @importFrom grDevices rgb col2rgb
makeTransparent<-function(color, alpha=100)
{
  newColor<-grDevices::col2rgb(color)
  apply(newColor, 2, function(curcoldata){
    grDevices::rgb(red=curcoldata[1], green=curcoldata[2],
                   blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}
