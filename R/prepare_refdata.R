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
#' interpol.larval <- prepare_refdata(ref="larval")
#' }
#'
#' @importFrom utils data
#' @importFrom limma normalizeBetweenArrays
prepare_refdata <- function(ref = c("young_adult", "Cel_YA_adult1", "sterken",
                                    "larval_development", "Cel_larval", "oudenaarden",
                                    "embryonic_development", "Cel_embryo", "hashimshony", 
                                    "Cel_YA_adult2", "reinke"),
                            n.inter = 200)
{
  ref <- match.arg(ref)
  
  if(ref=="larval_development"|ref=="Cel_larval"|ref=="oudenaarden"){
    message("Loading the C. elegans reference dataset for larval development")
    utils::data("Cel_larval", envir = environment())
    # ICA components with relevant time dynamics
    keeps <- (1:20)[-c(2,6,10,14)]
    # span values for loess regression of components
    sps <- c(0.35, 0.30, 0.50, 0.35, 0.25, 0.30, 
             0.25, 0.30, 0.35, 0.30, 0.25, 0.35)
    
    interp.dat <- interpol_refdata(Cel_larval$X, n.inter,
                                   time.series = Cel_larval$time.series,
                                   ica.nc = 20, center=T,
                                   keep.c = keeps, span = sps)
  }
  if(ref=="embryonic_development"|ref=="Cel_embryo"|ref=="hashimshony"){
    message("Loading the C. elegans reference dataset for embryonic development")
    utils::data("Cel_embryo", envir = environment())    
    # ICA components with relevant time dynamics
    keeps <- (1:16)[-c(6,10,14:16)]
    # span values for loess regression of components
    sps <- c(.25,.3,.2,.28,.2,.2,.2,.2,.18,.2,.25)
    
    interp.dat <- interpol_refdata(Cel_embryo$X, n.inter,
                                   time.series = Cel_embryo$time.series,
                                   ica.nc = 16, center = T,
                                   keep.c = keeps, span = sps)
  }
  if(ref=="young_adult"|ref=="Cel_YA_adult1"|ref=="sterken"){
    message("Loading the C. elegans reference dataset for young adult/adult worms")
    utils::data("Cel_larval", envir = environment())
    utils::data("Cel_YA_adult1", envir = environment())
    # Interpolation is done together with the larval dataset for 
    # better selection of gene expression dynamic components
    larvs <- cbind(Cel_larval$X, Cel_larval$X.25)[,names(Cel_larval$est.time.series)]
    ov <- format_to_ref(Cel_YA_adult1$X, larvs, verbose = F)
    X.os <- limma::normalizeBetweenArrays(cbind(ov$refdata, ov$samp), 
                                          method = "quantile")
    
    # ICA components with relevant time dynamics
    keeps <- (1:20)[-c(3, 7, 13, 14, 15, 19, 20)]
    # span values for loess regression of components
    sps <- c(0.6, 0.25, 0.4, 0.3, 0.15, 0.4, 0.3, 0.2, 0.2, 0.25, 0.25, 0.2, 0.15)
    
    interp.dat <- interpol_refdata(X.os, n.inter, 
                                   time.series = c(Cel_larval$est.time.series, 
                                                   Cel_YA_adult1$time.series),
                                   ica.nc = 20, center = T,
                                   keep.c = keeps, span = sps)
    
    
  }
  if(ref=="Cel_YA_adult2"|ref=="reinke"){
    message("Loading the C. elegans reference dataset for young adult/adult worms (2)\nNote : this reference is of much lower quality than Cel_YA_adult1")
    utils::data("Cel_YA_adult2", envir = environment())
    # Reinke data is already reconstructed from interpolated data, all 8 first
    # components are good.
    interp.dat <- interpol_refdata(Cel_YA_adult2$X, n.inter, 
                                   time.series = Cel_YA_adult2$time.series,
                                   ica.nc = 8, center=T,
                                   keep.c = 1:8, span = .1)
    
    
  }
  
  return(interp.dat)
}


