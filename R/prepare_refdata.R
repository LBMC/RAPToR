#' Prepare the reference data included in the package
#' 
#' This function loads (and interpolates on) the desired reference dataset.
#' The available datasets can be found in the [ref_table].
#' 
#' @param ref the name of the reference dataset to load ; can be abbreviated.
#' @param n.inter the resolution of the interpolation, as passed on to \code{\link{plsr_interpol}}.
#' 
#' @return the interpolated reference dataset, as returned by \code{\link{plsr_interpol}}
#' 
#' @seealso [plsr_interpol]
#' 
#' @export
#' 
#' @eval ae_example()
#'
#' @importFrom utils data
#' @importFrom limma normalizeBetweenArrays
prepare_refdata <- function(ref, n.inter = 200)
{
  utils::data("ref_table", envir = environment())  
  ref <- match.arg(arg = ref, choices = ref_table$name)
  
  dpkg <- ref_table[which(ref_table$name == ref), "data_pkg"]
  # check if needed data package is loaded
  if(!requireNamespace(dpkg, quietly = T)){
    stop(paste0("You must install the ", dpkg, " data package to load the ", ref, " reference"))
  }
  
  prep_func <- do.call(`::`, list(dpkg, paste0('.prepref_', ref)))
  r_i <- prep_func(n.inter = n.inter)

  return(r_i)
}


