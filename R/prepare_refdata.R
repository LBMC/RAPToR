#' Prepare the reference data included in the package
#' 
#' This function loads (and interpolates on) pre-built reference datasets.
#' The available datasets of a valid data-package can be listed with \code{\link{list_refs}}.
#' 
#' @param ref the name of the reference dataset to load (as given in \code{\link{list_refs}}).
#' @param datapkg the name of the data-package to load the rreference from.
#' @param n.inter the resolution of the interpolation, as in \code{seq(start, end, length.out = n.inter)}.
#' 
#' @return a list with the interpolated reference dataset and its associated time series. 
#' 
#' @seealso [list_refs] [ge_im]
#' 
#' @export
#' 
#' @eval ae_example()
#'
#' @importFrom utils getFromNamespace
prepare_refdata <- function(ref, datapkg,  n.inter = 200)
{
  requireNamespace(datapkg)
  
  if(exists(pr, where = asNamespace(datapkg), mode = 'function'){
    prepf <- paste0(".prepref_", ref)
    prepf <- utils::getFromNamespace(x = prepf, ns = datapkg)
    
    return(invisible(prepf(n.inter = n.inter)))
  } else {
    stop(paste0("Reference ", ref, " not found. Check list_refs(datapkg = '", datapkg, "') for a valid reference."))
  }
}


