#' Plot a dfCV object
#' 
#' Plots the `dfCV` object to help determine the optimal df parameter value for PLSR regression.
#' 
#' @param x a `dfCV` object, as returned by \code{\link{df_CV}}.
#' @param stat either "median" or "mean".
#' @param subset indices of a subset of df values to plot
#' @param t.plot boolean ; if TRUE (default), plots the individual CV Error trajectories.
#' @param t.col color for the individual CV Error trajectories.
#' @param signchange boolean ; if TRUE, displays '+'/'-' above the curve to indicate diff from previous point.
#' @param main title of the plot.
#' @param type the line type.
#' @param cex sizing parameter applied to various elements of the plot.
#' @param lwd the line width.
#' @param leg boolean ; if TRUE (default), displays a legend.
#' @param l.pos the position of the legend, as passed on to \code{\link[graphics]{legend}}
#' @param ... additional arguments passed on to \code{\link[graphics]{plot}}.
#' 
#' @export
#' 
#' 
#' @importFrom graphics plot points abline text legend
#' @importFrom stats median
#' 

plot.dfCV <- function(x, 
                      stat = c('median', 'mean'),
                      subset = 1:length(x$dfs),
                      t.plot = TRUE, t.col = makeTransparent('black', 30),
                      signchange = TRUE,
                      main = "CV Error", 
                      type = 'b', cex = 0.5, lwd = 2, 
                      leg = TRUE,
                      l.pos='top', 
                      ...)
{
  
  stat <- match.arg(stat)
  errs <- x$cv.errors[,subset]
  dfs <- x$dfs[subset]
  
  if('median' == stat)
    y <- apply(errs, 2, stats::median)
  if('mean' == stat)
    y <- apply(errs, 2, mean)
  
  graphics::plot(dfs, y, lwd = lwd, cex = cex, type = type,
                 xlab = "df", ylab = "CV Error", main = main, ...)
  
  # plot individual CV trajectories
  if(isTRUE(t.plot)){
    sapply(1:nrow(errs), function(i){
      graphics::points(dfs, errs[i,], type = type, col = t.col, 
                       lty = 3, cex = cex*.5)
    })
  }
  
  # display '+' or '-' above points to show gradient
  if(isTRUE(signchange))
    graphics::text(dfs[-1], y[-1], 
                   labels = c('-','+')[(diff(y)>0)+1], 
                   pos = 3, font = 2, cex = 2,
                   col = c('chartreuse4','firebrick')[(diff(y)>0)+1])
  
  # show minimum
  graphics::abline(h = min(y), lty = 2, col='royalblue')
  graphics::text(max(dfs), min(y), labels = "min", font = 2, 
                 col = 'royalblue', pos = 3, offset = .25)
  
  if(isTRUE(leg))
    legend(l.pos, bty = 'n', legend = stat, inset = .01,
           lwd = lwd, pt.cex = cex, pch = 1, lty=1, 
           cex = 1.2, text.font = 2)
}
