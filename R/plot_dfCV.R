#' Plot an ae object
#' 
#' Plots the age estimates along with bootstrap error bars using a dotchart.
#' 
#' @param x a \code{dfCV} object, as returned by \code{\link{df_CV}}.
#' @param col.p the color of the prior estimate marker.
#' @param col.b the color of the bootstrapped estimates.
#' @param pch the pch parameter passed on to \code{\link{dotchart}}.
#' @param cex sizing parameter applied to various elements of the plot.
#' @param xlim horizontal range for the plot, see \code{\link[graphics]{plot.window}}, for example
#' @param xlab the x axis label, passed on to \code{\link{plot}}.
#' @param l.pos the position of the legend when show.prior is \code{TRUE}, as passed on to \code{\link[graphics]{legend}}
#' @param ... additional arguments passed on to \code{\link{plot}}.
#' 
#' @export
#' 
#' @examples
#' data(Cel_larval)
#' 
#' samp <- Cel_larval$X[,13:15]
#' age.est <- estimate.worm_age(samp, Cel_larval$X, Cel_larval$time.series)
#' \donttest{
#' plot(age.est)
#' }
#' 
#' @importFrom graphics plot dotchart points arrows legend
#' @importFrom beeswarm swarmy
#' 

plot.dfCV <- function(x, 
                      stat = c('median', 'mean'),
                      subset = 1:length(x$dfs),
                      t.plot = TRUE, t.col = makeTransparent('black', 30),
                      main = "CV Error", l.pos='left', ...)
{
  
  stat <- match.arg(stat)
  M_errs <- x$cv_errors[,sel]
  
  y <- apply(M_errs, 2, median)
  if(!median)
    y <- colMeans(M_errs)
  plot(x$dfs[sel], y, lwd=2, cex=.5, type='b', 
       xlab = "df", ylab = "CV Error",
       main = main)
  invisible(sapply(1:nrow(M_errs), function(i){
    points(x$dfs[sel],M_errs[i,], type='b', col=t.col, lty=3, cex=.2)
  }))
  
  text((x$dfs[sel])[-1], y[-1], labels = c('-','+')[(diff(y)>0)+1], pos = 3, font = 2, 
       col=c('chartreuse4','firebrick')[(diff(y)>0)+1], cex = 2)
  
  abline(h=min(y), lty=2, col='royalblue')
  text(max(x$dfs), min(y), labels="min", font=2, col='royalblue', pos=3, offset = .25)
  legend(l.pos, bty='n', legend = ifelse(median, "median", "mean"), inset = .01,
         lwd=2, pt.cex = .5, text.font = 2, pch = 1, lty='65', cex=1.2)
}