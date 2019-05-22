#' Plot a timeline segment
#' 
#' Adds horizontal segments with notches at each end to the current plot.
#' 
#' @param x0,x1 coordinates of start and end of the segment.
#' @param y y coordinate of the segment.
#' @param notch.height length of the edge notches.
#' @param padd.x padding ; the segment is shortened by \code{padd.x/2} on each end
#' @param ... arguments passed on to \code{\link[graphics]{segments}}
#' 
#' 
#' @examples
#' \donttest{
#' plot(1:5, type='n', )
#' add.timeline(x0 = c(1,3), x1 = c(4,5), y = c(2,3), lwd=2)
#' }
#' 
#' @importFrom graphics segments
add.timeline <- function(x0, x1, y, notch.height=.5, padd.x=0, ...){
  n <- length(x0)
  if(length(x1)!=n)
    stop("x0 and x1 must have the same length")
  if(length(padd.x)<n)
    padd.x <- rep(padd.x, n)
  
  # apply padding
  x0 <- x0 + padd.x/2
  x1 <- x1 - padd.x/2
  
  # timeline
  graphics::segments(x0 = x0, x1 = x1, 
                     y0 = y, ...)
  
  # notches
  graphics::segments(x0 = x0, 
                     y0 = y - notch.height/2, y1 = y+notch.height/2, ...)
  graphics::segments(x0 = x1, 
                     y0 = y - notch.height/2, y1 = y+notch.height/2, ...)
}


#' Plot the reference dataset time series
#' 
#' Displays the reference datasets' coverage above key devlopmental stages.
#' The data used for plotting is in the \code{\link{ref_tables}} object.
#' 
#' 
#' @export
#' 
#' @examples
#' \donttest{
#' plot_ref_timelines()
#' }
#' 
#' @importFrom utils data
#' @importFrom graphics plot layout text par
plot_ref_timelines <- function(){
  par.save <- graphics::par(mar=c(4,2,3,2), mfrow=c(2,1))
  utils::data("ref_tables", envir = environment())
  graphics::layout(matrix(1:2), heights = c(.6,.4))
  
  ds <- ref_tables$datasets
  
  ## plot post.hatchings
  rf.ph <- ref_tables$post.hatch
  is.ph <- ds$unit=="h post-hatching"
  n.ph <- sum(is.ph)
  y.ph <- 2+(1:n.ph)
  
  graphics::plot(NULL, NULL, 
                 xlab = paste("time,", (ds[is.ph,"unit"])[1]),
                 xlim = range(ds[is.ph, c("from", "to")]), 
                 yaxt = 'n', ylab='',
                 ylim = c(0,3+n.ph),
                 main = "Reference datasets post-hatching")
  
  # datasets
  add.timeline(ds$from[is.ph], ds$to[is.ph], y.ph, lwd=3, 
               notch.height = .4, padd.x = .2)
  
  graphics::text((ds$from[is.ph]+ds$to[is.ph])/2, y.ph, pos=3, 
                 labels = paste(ds$data.name[is.ph], ' (', ds$name[is.ph], ')', sep=''), 
                 font=2)
  
  # reference timepoints
  add.timeline(rf.ph$from, rf.ph$to, 1, lwd=3, col='grey30', 
               notch.height = .25, padd.x = 0)
  
  graphics::text((rf.ph$from+rf.ph$to)/2, 1, offset = 1,
                 pos=3, labels = rf.ph$stage, font=2)
  
  
  ## plot embryo
  rf.e <- ref_tables$embryo
  is.e <- ds$unit=="min from 4C stage"
  n.e <- sum(is.e)
  y.e <- 2+(1:n.e)
  
  graphics::plot(NULL, NULL, 
                 xlab = paste("time,", (ds[is.e,"unit"])[1]),
                 xlim = range(ds[is.e, c("from", "to")]), 
                 yaxt = 'n', ylab='',
                 ylim = c(0,3.5+n.e),
                 main = "Reference datasets for embryo")
  
  # datasets
  add.timeline(ds$from[is.e], ds$to[is.e], y.e, lwd=3, 
               notch.height = .4, padd.x = .2)
  
  graphics::text((ds$from[is.e]+ds$to[is.e])/2, y.e, pos=3, offset = .5,
                 labels = paste(ds$data.name[is.e], ' (', ds$name[is.e], ')', sep=''),
                 font=2)
  
  # reference timepoints
  add.timeline(rf.e$from, rf.e$to, 1, lwd=3, col='grey30', 
               notch.height = .25, padd.x = 0)
  
  ups <- c(1,3,4,5,7) # to plot above, others below
  graphics::text(((rf.e$from+rf.e$to)/2)[ups], 1, offset = .5,
                 pos=3, labels = rf.e$stage[ups], font=2)
  graphics::text(((rf.e$from+rf.e$to)/2)[-ups], 1, offset = .5,
                 pos=1, labels = rf.e$stage[-ups], font=2)
  
  # restore graph. par
  graphics::par(par.save)
}

