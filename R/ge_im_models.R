
#' @importFrom stats prcomp as.formula update
#' @import mgcv
.model_gam_pca <- function(X, p, formula, nc = ncol(X), scale = T, center = T, drX = F){
  if(drX){
    pca <- X
    nc <- ncol(pca$rotation)
  }
  else pca <- stats::prcomp(X, rank = nc, scale = scale, center = center)
  
  formula <- stats::as.formula(formula)
  
  formulas <- lapply(seq_len(nc), function(i) stats::update(formula, paste0("PC", i, " ~ .")))
  colnames(pca$rotation) <- paste0("PC", seq_len(nc))
  p <- cbind(p, pca$rotation)
  
  m <- list()
  m$model <- lapply(formulas, mgcv::gam, data = p)
  m$pca <- pca
  
  return(m)
}

#' @import mgcv
.predict_gam_pca <- function(m, newdata, as.pc = FALSE){
  preds <- do.call(cbind, lapply(m$model, predict, newdata = newdata))
  
  if(!as.pc)
    return(tcrossprod(m$pca$x, preds))
  else return(preds)
}



#' @importFrom stats as.formula update 
#' @import mgcv
.model_gam_ica <- function(X, p, formula, nc = ncol(X), center = T, drX = F){
  if(drX){
    ica <- X
    nc <- ncol(ica$M)
  }
  else ica <- ica::icafast(X, nc = nc, center = center)
  
  formula <- stats::as.formula(formula)
  
  formulas <- lapply(seq_len(nc), function(i) stats::update(formula, paste0("IC", i, " ~ .")))
  colnames(ica$M) <- paste0("IC", seq_len(nc))
  p <- cbind(p, ica$M)
  
  m <- list()
  m$model <- lapply(formulas, mgcv::gam, data = p)
  m$ica <- ica
  
  return(m)
}

#' @import mgcv
.predict_gam_ica <- function(m, newdata, as.ic = FALSE){
  preds <- do.call(cbind, lapply(m$model, predict, newdata = newdata))
  
  if(!as.ic)
    return(tcrossprod(m$ica$S, preds))
  else return(preds)
}