#' List the references in a data-package
#' 
#' This function looks for the \code{ref_list} object from the input data-package and prints it.
#' 
#' @param datapkg the name of the data-package to list references from.
#' 
#' @export
#' 
#' @examples 
#' 
#' # list_refs(datapkg = "wormRef")
#' 
list_refs <- function(datapkg){
  requireNamespace(datapkg)
  rl <- "ref_list"
  if(exists(suppressWarnings(data(list = rl, package = datapkg, envir = environment())))){
    dat <- get(data(list = rl, package = datapkg, envir = environment()))
    rm(list = rl)
    print(dat)
    return(invisible(dat))
  } else {
    stop(paste0("Object `",rl,"` not found. Is ", datapkg, " a valid data-package for RAPToR ?"))
  }
}



#' Plot reference timelines
#' 
#' Plots a timeline chart for references of the input data-package.
#' This corresponds to a \code{.plot_refs} function in the data-package.
#' 
#' @param datapkg the name of the data-package to plot references from
#' @param ... extra arguments passed on to the datapkg's plot function
#' 
#' @export
#' 
#' @examples 
#' \donttest{
#' plot_refs("wormRef")
#' }
#' 
#' @importFrom utils data getFromNamespace
plot_refs <- function(datapkg, ...){
  requireNamespace(datapkg)
  
  
  pr <- ".plot_refs"
  if(exists(pr, where = asNamespace(datapkg), mode = 'function')){
    pf <- utils::getFromNamespace(x = pr, ns = datapkg)
    return(invisible(pf(...)))
  } else {
    stop(paste0("Function `",pr,"` not found. Is ", datapkg, " a valid data-package for RAPToR ?"))
  }
}
