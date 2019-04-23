#' Prepare the reference data included in the package
#' 
#' This function loads the desired reference dataset from the package and
#' performs the interpolation with optimal parameters.
#' Three available reference datasets are available
#' (embryonic development, larval development and young adult/adult).
#' 
#' @param ref which reference dataset to load. Can be abbreviated.
#' @param n.inter the \code{n.inter} parameter, as passed on to \code{\link{interpol_refdata}}.
#' 
#' @return the interpolated reference dataset, as would \code{\link{interpol_refdata}}
#' 
#' @seealso [interpol_refdata()]
#' 
#' @export
#' 
#' @examples
#' \donttest{ 
#' interpol.oudenaarden <- prepare_refdata(ref="larval")
#' }
#'
prepare_refdata <- function(ref = c("larval_development", "embryonic_development",
                                    "oudenaarden", "hashimshony", 
                                    "young_adult", "reinke"),
                            n.inter = 200)
{
  ref <- match.arg(ref)
  
  if(ref=="larval_development"|ref=="oudenaarden"){
    message("Loading the Oudenaarden reference dataset for larval development")
    data("oud_ref")
    # join the 20 and 25C series (and quantile normalize)
    X <- limma::normalizeBetweenArrays(cbind(oud_ref$X, oud_ref$X.25), method = "quantile")
    # ICA components with relevant time dynamics
    keeps <- (1:20)[-c(11,12,14,18,19)]
    # span values for loess regression of components
    sps <- c(.4,.3,.35,.3,.25,.25,.25,.3,.3,.25,.25, .2,.25, .3,.2)
      
    interp.dat <- interpol_refdata(X[,names(oud_ref$est.time.series)], n.inter,
                                   time.series = oud_ref$est.time.series,
                                   ica.nc = 20, center=T,
                                   keep.c = keeps, span = sps)
  }
  if(ref=="embryonic_development"|ref=="hashimshony"){
    message("Loading the Hashimshony reference dataset for embryonic development")
    data("hash_ref")    
    # ICA components with relevant time dynamics
    keeps <- (1:16)[-c(6,10,14:16)]
    # span values for loess regression of components
    sps <- c(.25,.3,.2,.28,.2,.2,.2,.2,.18,.2,.25)
    
    interp.dat <- interpol_refdata(hash_ref$X, n.inter,
                                   time.series = hash_ref$time.series,
                                   ica.nc = 16, center = T,
                                   keep.c = keeps, span = sps)
  }
  if(ref=="young_adult"|ref=="reinke"){
    message()
    data("reinke_ref")
    # Reinke data is already reconstructed from interpolated data, all 8 first
    # components are good.
    interp.dat <- interpol_refdata(reinke_ref$X, n.inter, 
                                   time.series = reinke_ref$time.series,
                                   ica.nc = 8, center=T,
                                   keep.c = 1:8, span = .1)
  }
  
  return(interp.dat)
}


