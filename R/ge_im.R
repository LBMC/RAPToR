#' Gene expression interpolation model
#' 
#' Build a model to interpolate on a gene expression dataset. 
#' This can be done either with gam or glm models fit on the components of a PCA or ICA.
#' It's also possible to have a linear model fit directly (per gene) on the gene expression data (uses limma).
#' 
#' We use components as "eigen genes" to find model parameters fitting the whole 
#' gene set \insertCite{storey2005significance}{RAPToR}.
#' 
#' @param X the gene expression matrix (genes as rows, samples as columns)
#' @param p a dataframe with the phenotypic data used in the formula (samples as rows) e.g. time, covariates.
#' @param formula the model formula, which must start with 'X ~'. See \code{\link[mgcv]{gam}}, \code{\link[stats]{glm}} or \code{\link[limma]{lmFit}} documentation for specifications.
#' @param method the model to fit, one of c("gam", "glm", "limma").
#' @param dim_red the dimension reduction method to use for interpolation, one of c("pca", "ica"), ignored if method is "limma".
#' @param nc the number of components to extract from \code{X} for interpolation, defaults to \code{ncol(X)}, ignored if method is "limma".
#' @param ... extra arguments passed on to model functions.
#'
#' @return a '\code{geim}' model object. This object has its predict method
#' 
#' @export
#' 
#' @eval interpol_example()
#' 
#' @importFrom Rdpack reprompt
#' 
ge_im <- function(X, p, formula, 
                  method = c("gam", "glm", "limma"), 
                  dim_red = c("pca", "ica"), 
                  nc = ncol(X), 
                  ...)
{
  method <- match.arg(method)
  if("limma" == method){
    dim_red <- NA
    nc <- NA
  } else {
    dim_red <- match.arg(dim_red)
  }
  
  vars <- mgcv::interpret.gam(as.formula(formula))$pred.names
  
  # extract time variable (it must be the only numeric variable in the formula)
  t.var <- which(sapply(vars, function(v) is.numeric(p[,v])))[1]
  
  # extract covariates
  cvars <- vars[-t.var]
  
  # keep t.var name
  t.var <- vars[t.var]
  
  # check that (potential) covariates are factors
  if(length(cvars) > 0){
    if(any(!sapply(cvars, function(v) is.factor(p[, v])))){
      stop(paste0("Any specified covariate must be a factor (",paste(cvars, collapse = ", "), ")"))
    }
  }

  if("gam" == method){
    if("pca" == dim_red){
      m <- .model_gam_pca(X = X, p = p, formula = formula, nc = nc, ...)
    } else if("ica" == dim_red){
      if(1==nc){
        stop("nc=1: 1-component ICA is impossible, please use PCA")
      }
      m <- .model_gam_ica(X = X, p = p, formula = formula, nc = nc, ...)
    }
  } else if("glm" == method){
    if("pca" == dim_red){
      m <- .model_glm_pca(X = X, p = p, formula = formula, nc = nc, ...)
    } else if("ica" == dim_red){
      if(1==nc){
        stop("nc=1: 1-component ICA is impossible, please use PCA")
      }
      m <- .model_glm_ica(X = X, p = p, formula = formula, nc = nc, ...)
    }
  } else if("limma" == method){
    m <- .model_limma(X = X, p = p, formula = formula)
  } 
  
  class(m) <- "geim"
  attr(m, "pdata") <- p
  attr(m, "formula") <- formula
  attr(m, "vars") <- list(t.var = t.var, cvars = cvars)
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
#' @param ... ignored (needed to match the S3 standard)
#' 
#' @export
#' 
predict.geim <- function(object, newdata, as.c = FALSE, ...){
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




#' Print a geim object
#' 
#' Prints a \code{geim} object
#' 
#' @param x a \code{geim} object, as returned by \code{\link{ge_im}}.
#' @param ... arguments passed on to \code{\link{print}}
#' 
#' @export
#' 
#' 
print.geim <- function(x, ...){
  ats <- list(f = attr(x, "formula"), m = attr(x, "method"), 
              dr = attr(x, "dim_red"), nc = attr(x, "nc"))
  cat("RAPToR Gene Expression Interpolation Model \n---")
  if(ats$m != "limma"){
    cat("\n\t", casefold(ats$m, upper = T), "fit on", ats$nc,
        casefold(ats$dr, upper = T), "components:")
  } else {
    cat("\n\t", casefold(ats$m, upper = T), " model fit on genes with:")
  }
  cat("\n\t", ats$f)
  cat("\n---\n")
}



#' Plot a GEIM object
#' 
#' Plots geim interpolation on components.
#' 
#' @param x a `geim` object, as returned by \code{\link{ge_im}}.
#' @param ref a \code{ref} object, as returned by \code{\link{make_ref}}.
#' @param ncs which components to plot, defaults to all.
#' @param n.inter interpolation resolution, as in \code{seq(start, end, length.out = n.inter)}. One of \code{ref}, \code{n.inter}, or \code{by.inter} must be specified.
#' @param by.inter interpolation resolution, as in \code{seq(start, end, by = by.inter)}. One of \code{ref}, \code{n.inter}, or \code{by.inter} must be specified.
#' @param cov.levels a named list with potential model covariate levels (e.g batch, strain) to predict as (defaults to first level).
#' @param t.unit an optional string specifying the time unit and t-zero, e.g "h past egg-laying".
#' @param col,col.i color for data and interpolation line respectively.
#' @param show.legend,l.pos whether to show the legend, and its position passed on to \code{\link[graphics]{legend}}
#' @param show.stats whether to print model fit (deviance explained and relative error) and component (variance explained) statistics on the plot.
#' @param ... additional arguments passed on to \code{\link[graphics]{plot}}.
#' 
#' @return Invisibly returns a table with component variance explained and model fit on component indices (R2, deviance explained, relative error).
#' 
#' @export
#' 
#' @eval interpol_example()
#' 
#' @importFrom graphics legend mtext points plot
#' 
plot.geim <- function(x, ref=NULL, ncs=NULL, 
                      n.inter=NULL, by.inter=NULL, cov.levels=NULL, t.unit=NULL, 
                      col=NULL, col.i='red', 
                      show.legend=T, l.pos='topleft', show.stats=T,...){
  p <- attr(x, "pdata")
  t.var <- attr(x, "vars")$t.var # time variable
  tv <- p[, t.var]
  cvars <- attr(x, "vars")$cvars # covariates
  dr <- attr(x, "dim_red")
  nc <- attr(x, "nc")
  
  if(is.null(ncs)){
    ncs <- 1L:nc
  } else if(any(! ncs%in%(1L:nc))){
    stop(paste0("ncs must be among components fit by model (1-",nc,"."))
  }
  
  if(!is.null(ref)){ # if ref is supplied
    if(class(ref)!="ref"){
      stop("'ref' must be a ref object.")
    }
    ts <- ref$time
    cov.levels <- ifelse(is.null(cov.levels), attr(ref, "cov.level"))
    t.unit <- ifelse(is.null(t.unit), attr(ref, "t.unit"))
    
  } else {
    if(is.null(n.inter) & is.null(by.inter)){
      stop("One of ref, n.inter, or by.inter must be specified.")
    }
    cov.levels <- .cov_check(x, cov.levels)
    if(is.null(t.unit)){
      t.unit <- "no unit specified"
    }
    # n/by.inter param handling
    if(!is.null(by.inter)){
      if(!is.numeric(by.inter)){
        stop("by.inter must be a single numeric value.")
      }
      ts <- seq(min(tv), max(tv), by = by.inter)
    } else {
      ts <- seq(min(tv), max(tv), l = n.inter)
    }
    
  }
  # make the new predictor dataframe
  fit_dat <- rep(T, nrow(p)) # keep track of points used for prediction
  l <- length(ts)
  ndat <- data.frame(time = ts)
  colnames(ndat)[1] <- t.var
  for(v in cvars){
    ndat[, v] <- factor(rep(cov.levels[[v]], l=l))
    fit_dat[ (p[,v] != cov.levels[v])] <- F # unused for prediction
  }
  
  if(dr == "pca"){
    yl <- "PC"
    y0df <- x[[dr]]$x
    vexp <- summary(x[[dr]])$importance[2,1L:length(x$model)]
  } else {
    yl <- "IC"
    y0df <- x[[dr]]$M
    vexp <- x[[dr]]$vafs
  }
  user_col <- F
  if(is.null(col)){
    col <- c(makeTransparent(par("fg")), par("fg"))[1+fit_dat]
  } else if(length(col)==1){
    col <- c(makeTransparent(col), col)[1+fit_dat]
  } else {
    warning("User-defined colors: removing elements from legend.")
    user_col <- T
  }
  
  xl <- paste0("Reference time, ", t.unit)
  y1df <- predict(x, ndat, as.c=T)
  
  m_stats <- do.call(rbind, lapply(x$model, function(m){
    c(r2 = 1 - sum(m$residuals^2)/sum(m$y**2),
      deviance.expl = (m$null.deviance - m$deviance)/m$null.deviance,
      relative.err = sum(abs(m$residuals/m$y))/length(m$y))
  }))
  m_stats <- cbind(component.var.exp=vexp, m_stats)
  rownames(m_stats) <- paste0(yl, 1L:length(x$model))
  
  invisible(sapply(ncs, function(i){
    yli <- ifelse(show.stats, paste0(yl, i, ' (', round(100*m_stats[i, "component.var.exp"], 2) , '%)'), paste0(yl, i))
    graphics::plot(tv, y0df[,i], 
         ylab = yli, xlab = xl, 
         lwd=2, col = col, ...)
    graphics::points(ts, y1df[,i], type='l', lwd=2, col = col.i)
    if(i==ncs[1] & show.legend){
      if(user_col){
        graphics::legend(l.pos, bty='n', legend=c("data", "model fit"), 
               lty=c(NA, 1), lwd=c(2,3), 
               col = c(col[1], col.i),
               pch=c(1, NA))
      } else if(all(fit_dat)){
        graphics::legend(l.pos, bty='n', legend=c("data", "model fit"), 
               lty=c(NA, 1), lwd=c(2,3), 
               col = c(col[1], col.i),
               pch=c(1, NA))
      } else {
        graphics::legend(l.pos, bty='n', legend=c("data (all)", "data (selected cov.)", "model fit"), 
               lty=c(NA, NA,1), lwd=c(2,2,3), col = c(unique(col), col.i),
               pch=c(1, 1, NA))
      }
    }
    if(show.stats){
      graphics::mtext(paste0("DE = ", round(m_stats[i,"deviance.expl"], 4)), side = 3, adj = 0,
            at=par("usr")[1], line = 0, cex = .75)
      graphics::mtext(paste0("RE = ", round(m_stats[i,"relative.err"], 4)), side = 3, adj = 1,
            at=par("usr")[2], line = 0, cex = .75)
    }
    
  }))
  
  return(invisible(m_stats))
}