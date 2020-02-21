#' Gene expression interpolation model
#' 
#' Build a model to interpolate on a gene expression dataset. 
#' This can be done either with gam or glm models fit on the components of a PCA or ICA.
#' It's also possible to have a limma model fit directly on the gene expression data.
#' 
#' @param X the gene expression matrix (genes as rows, samples as columns)
#' @param p a dataframe with the pheno data used in the formula (samples as rows) e.g. time, covariates.
#' @param formula the model formula, which must start with 'X ~'. See \code{\link[mgcv]{gam}}, \code{\link[stats]{glm}} or \code{\link[limma]{lmFit}} documentation for specifications.
#' @param method the model to fit, one of c("gam", "glm", "limma").
#' @param dim_red the dimension reduction method to use for interpolation, one of c("pca", "ica"), ignored when method is "limma".
#' @param nc the number of components to extract from \code{X} for interpolation, defaults to \code{ncol(X)}, ignored method is "limma".
#' @param ... extra arguments passed on to model functions.
#'
#' @return a '\code{geim}' model object. This object has its predict method
#' 
#' @export
#' 
ge_im <- function(X, p, formula, 
                  method = c("gam", "glm", "limma"), 
                  dim_red = c("pca", "ica"), 
                  nc = ncol(X), ...)
{
  method <- match.arg(method)
  if("limma" == method){
    dim_red <- NA
    nc <- NA
  } else {
    dim_red <- match.arg(dim_red)
  }
  
  if("gam" == method){
    if("pca" == dim_red){
      m <- .model_gam_pca(X = X, p = p, formula = formula, nc = nc, ...)
    } else if("ica" == dim_red){
      m <- .model_gam_ica(X = X, p = p, formula = formula, nc = nc, ...)
    }
  } else if("glm" == method){
    if("pca" == dim_red){
      m <- .model_glm_pca(X = X, p = p, formula = formula, nc = nc, ...)
    } else if("ica" == dim_red){
      m <- .model_glm_ica(X = X, p = p, formula = formula, nc = nc, ...)
    }
  } else if("limma" == method){
    m <- .model_limma(X = X, p = p, formula = formula)
  } 
  
  class(m) <- "geim"
  attr(m, "formula") <- formula
  attr(m, "method") <- method
  attr(m, "dim_red") <- dim_red
  attr(m, "nc") <- nc
  
  return(m)
}

#' Predictions of a geim model
#' 
#' Returns the model predictions for given data.
#' 
#' @param object a geim object, as returned by \code{\link{ge_im}}
#' @param newdata a dataframe with values of the model variables for prediction, defaults to model input.
#' @param as.c boolean ; if TRUE, returns predictions as components of the dimensionally reduced data.
#' 
#' @export
#' 
predict.geim <- function(object, newdata, as.c = FALSE){
  method <- attr(object, "method")
  dim_red <- attr(object, "dim_red")
  
  if("gam" == method){
    
    if("pca" == dim_red){
      return(.predict_gam_pca(m = object, newdata = newdata, as.pc = as.c))
    } else if("ica" == dim_red){
      return(.predict_gam_ica(m = object, newdata = newdata, as.ic = as.c))
    }
  } else if("glm" == method){
    
    if("pca" == dim_red){
      return(.predict_glm_pca(m = object, newdata = newdata, as.pc = as.c))
    } else if("ica" == dim_red){
      return(.predict_glm_ica(m = object, newdata = newdata, as.ic = as.c))
    }
  } else if("limma" == method){
    
    return(.predict_limma(m = object, newdata = newdata))
  }
}