#' Prepare the reference data included in the package
#' 
#' This function loads the desired reference dataset from the package and
#' performs the interpolation with optimal parameters.
#' Four available reference datasets are available
#' (embryonic development, larval development and young adult/adult - sterken or reinke).
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
#' @importFrom utils data
#' @importFrom limma normalizeBetweenArrays
prepare_refdata <- function(ref = c("young_adult", "sterken",
                                    "larval_development", "oudenaarden",
                                    "embryonic_development", "hashimshony", 
                                    "reinke"),
                            n.inter = 200)
{
  ref <- match.arg(ref)
  
  if(ref=="larval_development"|ref=="oudenaarden"){
    message("Loading the Oudenaarden reference dataset for larval development")
    utils::data("oud_ref", envir = environment())
    # join the 20 and 25C series (and quantile normalize)
    X <- limma::normalizeBetweenArrays(cbind(oud_ref$X, oud_ref$X.25), 
                                       method = "quantile")
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
    utils::data("hash_ref", envir = environment())    
    # ICA components with relevant time dynamics
    keeps <- (1:16)[-c(6,10,14:16)]
    # span values for loess regression of components
    sps <- c(.25,.3,.2,.28,.2,.2,.2,.2,.18,.2,.25)
    
    interp.dat <- interpol_refdata(hash_ref$X, n.inter,
                                   time.series = hash_ref$time.series,
                                   ica.nc = 16, center = T,
                                   keep.c = keeps, span = sps)
  }
  if(ref=="young_adult"|ref=="sterken"){
    message("Loading the Sterken reference dataset for young adult/adult worms")
    utils::data("oud_ref", envir = environment())
    utils::data("sterken_ref", envir = environment())
    # Interpolation is done together with the Oudenaarden dataset for 
    # better selection of gene expression dynamic components
    ouds <- cbind(oud_ref$X, oud_ref$X.25)[,names(oud_ref$est.time.series)]
    ov <- format_to_ref(sterken_ref$X, ouds, verbose = F)
    X.os <- limma::normalizeBetweenArrays(cbind(ov$refdata, ov$samp), 
                                          method = "quantile")
    
    # ICA components with relevant time dynamics
    keeps <- (1:20)[-c(3,8,11,13,16,18,19,20)]
    # span values for loess regression of components
    sps <- c(0.6, 0.25, 0.4, 0.3, 0.3, 0.15, 0.2, 0.2, 0.25, 0.25, 0.15, 0.2)
    
    interp.dat <- interpol_refdata(X.os, n.inter, 
                                   time.series = c(oud_ref$est.time.series, 
                                                   sterken_ref$time.series),
                                   ica.nc = 20, center = T,
                                   keep.c = keeps, span = sps)
    
    
  }
  if(ref=="reinke"){
    message("Loading the Reinke reference dataset for young adult/adult worms")
    utils::data("reinke_ref", envir = environment())
    # Reinke data is already reconstructed from interpolated data, all 8 first
    # components are good.
    interp.dat <- interpol_refdata(reinke_ref$X, n.inter, 
                                   time.series = reinke_ref$time.series,
                                   ica.nc = 8, center=T,
                                   keep.c = 1:8, span = .1)
    
    
  }
  
  return(interp.dat)
}


