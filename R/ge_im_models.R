
#' @importFrom stats prcomp as.formula update
#' @import mgcv
.model_gam_pca <- function(X, p, formula, nc = ncol(X), drX = F){
  if(drX){
    pca <- X
    nc <- ncol(pca$x)
  }
  else {
    tXc <- scale(t(X), scale = FALSE, center = TRUE) 
    pca <- stats::prcomp(tXc, rank = nc, scale = FALSE, center = FALSE)
    pca$gcenters <- attr(tXc, "scaled:center")
  }
  
  formula <- stats::as.formula(formula)
  
  formulas <- lapply(seq_len(nc), function(i) stats::update(formula, paste0("PC", i, " ~ .")))
  colnames(pca$x) <- paste0("PC", seq_len(nc))
  p <- cbind(p, pca$x)
  
  m <- list()
  m$model <- lapply(formulas, mgcv::gam, data = p)
  m$pca <- pca
  
  return(m)
}

#' @import mgcv
.predict_gam_pca <- function(m, newdata, as.pc = FALSE){
  preds <- do.call(cbind, lapply(m$model, predict, newdata = newdata))
  
  if(!as.pc)
    return(apply(tcrossprod(m$pca$rotation, preds), 2, function(co) co + m$pca$gcenters))
  else return(preds)
}



#' @importFrom stats as.formula update 
#' @importFrom ica icafast
#' @import mgcv
.model_gam_ica <- function(X, p, formula, nc = ncol(X), drX = F){
  if(drX){
    ica <- X
    nc <- ncol(ica$M)
  }
  else {
    Xc <- t(scale(t(X), scale = FALSE, center = TRUE))
    ica <- ica::icafast(Xc, nc = nc, center = TRUE)
    ica$gcenters <- attr(Xc, "scaled:center")
  }
  
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
    return(apply(tcrossprod(m$ica$S, preds), 2, function(co) co + m$ica$gcenters))
  else return(preds)
}



#' @importFrom stats prcomp as.formula update glm
.model_glm_pca <- function(X, p, formula, nc = ncol(X), drX = F){
  if(drX){
    pca <- X
    nc <- ncol(pca$x)
  }
  else {
    tXc <- scale(t(X), scale = FALSE, center = TRUE) 
    pca <- stats::prcomp(tXc, rank = nc, scale = FALSE, center = FALSE)
    pca$gcenters <- attr(tXc, "scaled:center")
  }  
  formula <- stats::as.formula(formula)
  
  formulas <- lapply(seq_len(nc), function(i) stats::update(formula, paste0("PC", i, " ~ .")))
  colnames(pca$x) <- paste0("PC", seq_len(nc))
  p <- cbind(p, pca$x)
  
  m <- list()
  m$model <- lapply(formulas, stats::glm, data = p)
  m$pca <- pca
  
  return(m)
}

.predict_glm_pca <- function(m, newdata, as.pc = FALSE){
  preds <- do.call(cbind, lapply(m$model, predict, newdata = newdata))
  
  if(!as.pc)
    return(apply(tcrossprod(m$pca$rotation, preds), 2, function(co) co + m$pca$gcenters))
  else return(preds)
}


#' @importFrom stats as.formula update glm
#' @importFrom ica icafast
.model_glm_ica <- function(X, p, formula, nc = ncol(X), drX = F){
  if(drX){
    ica <- X
    nc <- ncol(ica$M)
  }
  else {
    Xc <- t(scale(t(X), scale = FALSE, center = TRUE))
    ica <- ica::icafast(Xc, nc = nc, center = TRUE)
    ica$gcenters <- attr(Xc, "scaled:center")
  }
  
  formula <- stats::as.formula(formula)
  
  formulas <- lapply(seq_len(nc), function(i) stats::update(formula, paste0("IC", i, " ~ .")))
  colnames(ica$M) <- paste0("IC", seq_len(nc))
  p <- cbind(p, ica$M)
  
  m <- list()
  m$model <- lapply(formulas, stats::glm, data = p)
  m$ica <- ica
  
  return(m)
}

.predict_glm_ica <- function(m, newdata, as.ic = FALSE){
  preds <- do.call(cbind, lapply(m$model, predict, newdata = newdata))
  
  if(!as.ic)
    return(apply(tcrossprod(m$ica$S, preds), 2, function(co) co + m$ica$gcenters))
  else return(preds)
}



#' @importFrom stats as.formula update model.matrix model.frame terms
#' @importFrom limma lmFit
.model_limma <- function(X, p, formula){
  # lmFit takes 'X' as an argument outside of the formula
  formula <- stats::update(stats::as.formula(formula), NULL ~ .)
  
  dsgn <- stats::model.matrix(formula, data = p)
  m <- limma::lmFit(X, design = dsgn)
  
  m$mterms <- stats::terms(stats::model.frame(formula = formula, data = p))
  
  return(m)
}

#' @importFrom stats model.matrix
.predict_limma <- function(m, newdata){
  dsgn <- stats::model.matrix(m$mterms, newdata)
  return(tcrossprod(m$coefficients, dsgn))
}