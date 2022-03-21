#' Get reference expression profiles matching sample age estimates.
#' 
#' Fetches the (indices of, or) expression profiles of a reference that correspond to the input age estimates.
#' Matching reference data can be used to quantify the developmental signal between 
#' experimental groups where the variable of interest is confounded by development.
#' 
#' @param ref a \code{ref} object (as returned by \link{\code{make_ref}})
#' @param ae_obj an \code{ae} object (as returned by \link{\code{ae}})
#' @param ae_values age estimate values (used if ae_obj is NULL). 
#' @param return.idx if TRUE (default) returns reference indices. If FALSE, returns reference expression matrix
#'
#' @return either a vector of indices (if \code{return.idx} is True) or an expression matrix for selected reference time points.
#' 
#' @export
#' 
#' @eval interpol_example()
#' 
#' @importFrom Rdpack reprompt
#' 
get_refTP <- function(ref, ae_obj=NULL, ae_values=NULL, 
                      return.idx = TRUE){
  if("ref" != class(ref)){
    stop("ref must be an object of class 'ref', as returned by 'make_ref()'.")
  }
  if(is.null(ae_obj)){
    if(is.null(ae_values)){
      stop("One of ae_obj or ae_values must be specified.")
    } else{
      if("ae" != class(ae_obj)){
        stop("ae_obj must be an object of class 'ae', as returned by 'ae()'")
      }
      ae_values <- ae_obj$age.estimates[,1]
    }
  }
  
  idx <- sapply(ae_values, function(t) which.min(abs(ref$time.series-t)))
  if(ret.idx){
    return(idx)
  }
  return(ref$interpGE[,idx])
}




