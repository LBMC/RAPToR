#' Plot an ae object
#' 
#' Plots age estimates of samples with bootstrap confidence intervals.
#' 
#' @param x an `ae` object, as returned by \code{\link{ae}}.
#' @param groups a factor with sample categories (e.g. treatment groups).
#' @param subset an index vector of the samples to plot (defaults to all).
#' @param show.boot_estimates logical ; if TRUE (default), shows individual bootstrap estimates as swarms.
#' @param show.prior logical ; if TRUE, shows input prior(s) on the plot.
#' @param col,color the estimate, confidence interval, and label color.
#' @param col.b,col.p,col.l color of bootstrap estimates, priors, and background lines respectively.
#' @param CIbar.width the width of the confidence interval bars.
#' @param truncate_name whether to truncate displayed sample names from start, end, or not (default).
#' @param sn_len number of characters to keep when truncating sample names.
#' @param lmar left margin value, increase to fit sample names.
#' @param g.line position of the group names (margin line).
#' @param glob.above logical ; if TRUE, the global estimate is overlayed above all else.
#' @param pch,cex,xlim graphical parameters.
#' @param xlab x axis label.
#' @param main plot title.
#' @param l.pos position of the legend when show.prior is \code{TRUE}, passed on to \code{\link[graphics]{legend}}
#' @param ... additional arguments passed on to \code{\link{points}}.
#' 
#' @export
#' 
#' @eval ae_example()
#' 
#' @importFrom graphics arrows axis box legend plot points title   
#' @importFrom beeswarm swarmy
#' 
plot.ae <- function(x, groups=NULL, subset=NULL,
         show.boot_estimates=T, show.prior=F, 
         col = par("fg"), color = col,
         col.b=2, col.p=1, col.l='gray',
         pch=16, cex=1,
         truncate_name=c("none", "end", "start"), 
         sn_len = 10, lmar=10, g.line=lmar*.75,
         xlim=NULL, xlab=NULL, main=NULL,
         CIbar.width=0.1, glob.above = F,
         l.pos='bottomright', ...)
{
  op <- par("mar", "fg")
  on.exit(par(op))
  
  # subset the data to plot
  subset <- c(na.omit(subset))
  if(is.null(subset)){
    subset <- 1L:nrow(x$age.estimates)
  } else if(!identical(subset, unique(subset))| any(subset==0) | length(subset)==0){
    stop("subset must be a vector of unique and valid indices.")
  }
  
  a.e <- x$age.estimates[subset,,drop=F]
  prior <- x$prior[subset]
  
  if(!is.null(groups)){
    groups <- droplevels(groups[subset])
  }
  
  ci.inf <- a.e[,2]
  ci.sup <- a.e[,3]
  n <- nrow(a.e)
  truncate_name <- match.arg(truncate_name)
  
  if(is.null(xlim)){
    xlim <- range(c(ci.inf, ci.sup, prior))
  }
  if(is.null(xlab)){
    xlab <- paste0("Reference time, ", attr(x, "t.unit"))
  }
  
  # adjust Y positions to group layout
  y <- 1L:n
  o <- y
  
  if(!is.null(groups)){
    o <- sort.list(as.numeric(groups), decreasing = TRUE)
    ci.inf <- ci.inf[o]
    ci.sup <- ci.sup[o]
    groups <- groups[o]
    offset <- cumsum(c(0, diff(as.numeric(groups)) != 0))
    y <- 1L:n + 2 * offset
  }
  color <- rep(color, n)[o]
  
  # prepare plot window
  par(mar=c(op$mar[1], lmar, op$mar[3:4]))
  graphics::plot.new()
  graphics::plot.window(xlim = xlim, ylim = range(0.5, y+1))
  graphics::axis(1)
  graphics::abline(h=y, lty = "dotted", col = 'gray')
  
  # plot ae
  graphics::points(a.e[o,1], y, 
                   col=color, pch=pch, cex=cex, ...)
  
  # manage sample/group labels
  nms <- rownames(a.e)
  if(truncate_name == "start"){
    nms <- strtrim(nms, sn_len)
  }
  if(truncate_name == "end"){
    nms <- substr(nms, nchar(nms) - (sn_len-1), nchar(nms))
  }
  mtext(nms[o], side = 2, line=1, at=y, adj = 1, col = color, las=1)
  
  # group labels
  if(!is.null(groups)){
    lvg <- levels(groups)
    ylvg <-  rev(cumsum(rev(tapply(groups, groups, length)) + 2) - 1)
    mtext(text = lvg, side = 2, at = ylvg, line = g.line, las=1, font=2, adj = 1)
  }
  
  
  # plot confidence intervals
  graphics::arrows(ci.sup, y,
                   ci.inf, y,
                   angle=90, code=3, length=CIbar.width, 
                   col = color)
  
  # adding individual bootstrap estimates as swarms
  if(show.boot_estimates){
    boots <- x$boots[,subset,,drop=F]
    nboot <- dim(boots)[3]
    xs <- rbind(boots[1,o,])
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
  
  # adding prior to plot
  if(show.prior){
    inis <- prior[o]
    col.p <- rep(col.p, n)
    col.p <- col.p[o]
    graphics::points(inis, y, lwd=2, cex=cex*1.1, col=col.p)
    graphics::legend(l.pos, legend = ' Prior', col = col.p[n], inset = .02,
                     pt.lwd=2, pch=1, bty = 'n', text.col = col.p[n])
  }
  
  if(glob.above){
    graphics::points(a.e[o,1], y, cex=cex, pch=pch, col = color, ...)
  }
  
  graphics::box()
  graphics::title(xlab=xlab, main=main)
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
    
    if(is.null(xlab)){
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
      graphics::points(init.est, yl[1], pch = "|", 
                       col = col.p, cex = bar.size)
      graphics::text(init.est, yl[1], pos = 3, offset = 1, 
                     labels = paste(round(init.est, 2), 
                                    "\n(prior)", sep = ""))
    }
  })
}

#' (DEPRECATED) Plot an ae object
#' 
#' This function has been renamed. Please use \code{\link{plot_cor}} instead. 
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
