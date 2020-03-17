#' Cross validation comparison of geim models 
#' 
#' This function performs cross-validation (CV) with the aim of finding the optimal model from the given formulas.
#' Parameters are explored through the given list of formulas (*e.g*, 'df=3' or 'df=4' must be specified in the formulas).
#' 
#' The CV training sets are defined to be representative of all variables included in the models. 
#' This is done with a function attributed to GitHub user mrdwab \url{https://gist.github.com/mrdwab/6424112}.
#' 
#' Note that only one method/dimension reduction can be used at a time through this function.
#' 
#' @param X the gene expression matrix (genes as rows, samples as columns)
#' @param p a dataframe with the pheno data used in the formula (samples as rows) e.g. time, covariates.
#' @param formula_list a list of model formulas to compare, which must start with 'X ~' (as passed on to \code{\link{ge_im}}).
#' @param cv.n number of cross-validation repeats.
#' @param cv.s ratio of samples to use for training set. If `cv.s > 1`, then `cv.s` samples are used for the training set.
#' @param method the model type to fit, one of c("gam", "glm", "limma").
#' @param dim_red the dimension reduction method to use for interpolation, one of c("pca", "ica"), ignored when method is "limma".
#' @param nc the number of components to extract from \code{X} for interpolation, defaults to \code{ncol(X)}, ignored method is "limma".
#' @param to_compute the model performance indices to compute during CV (see \code{\link{mperf}})
#' @param nb.cores the number of cores to use for parallel execution.
#' @param verbose boolean ; if TRUE, displays messages of the various steps of the method.
#' @param ... extra arguments passed on to model functions.
#' 
#' 
#' @export
#' 
#' @eval interpol_example()
#' 
#' @importFrom stats as.formula prcomp
#' @importFrom ica icafast
#' @importFrom parallel parLapply makeForkCluster clusterExport stopCluster  
#' 
ge_imCV <- function(X, p, formula_list, cv.n = 50, cv.s = 0.8, 
                    method = c("gam", "glm", "limma"), dim_red = c("pca", "ica"), nc = ncol(X), 
                    to_compute = c("aRE", "MSE", "aRMSE"), 
                    nb.cores = 2, verbose = T, ...)
{
  if(verbose)
    message(paste0("CV on ", length(formula_list), " models. cv.n = ", cv.n, " | cv.s = ", cv.s))
  
  n <- ncol(X)
  nt <- ifelse(cv.s >=1 , cv.s, round(n * cv.s))
  rownames(p) <- NULL
  
  formula_list <- lapply(formula_list, stats::as.formula)
  
  method <- match.arg(method)
  dim_red <- ifelse("limma" == method, NA, match.arg(dim_red))
  if("limma" != method){
    if("pca" == dim_red){
      Xr <- stats::prcomp(X, rank = nc)
    } else if("ica" == dim_red){
      Xr <- ica::icafast(X, nc = nc)
    }
    X <- scale(X)
  } else {
    Xr <- X
  }
  
  if(verbose)
    message("\n...Building training sets")
  
  vlist <- lapply(formula_list, all.vars)
  vlist <- lapply(vlist, `[`, -1) # remove 'X'
  if(sum(duplicated.default(vlist)) != length(vlist) - 1)
    warning("Variable set differs between formulas.\n  The CV training set is built as representative of the first formula's variable set")
  
  plist <- lapply(vlist, function(vl) cbind(p[, vl, drop = F], dummycol = 0L))
  
  
  grpvar <- sapply(vlist[[1]], function(vi) is.factor(p[,vi]))
  
  if(all(!grpvar)){ # if no group variable
    cv.tsets <- lapply(seq_len(cv.n), function(i) sample.int(n, size = nt))
  } else {
    cv.tsets <- lapply(seq_len(cv.n), function(i) as.numeric(rownames(.stratified(plist[[1]], vlist[[1]][grpvar], cv.s))))
  }
  
  if(verbose)
    message("...Setting up cluster")
  cl <- parallel::makeCluster(nb.cores, type = ifelse(.Platform$OS.type == "windows", "PSOCK", "FORK"))
  parallel::clusterExport(cl, varlist = c("plist", "X", "Xr", "method", "dim_red",
                                          "nc", "formula_list", "cv.tsets"), envir = environment())
  
  if(verbose)
    message("...Running CV")
  res <- 
    parLapply(cl = cl,
              # lapply(
              seq_len(cv.n), function(i){
                tset <- cv.tsets[[i]]
                if("limma" != method){
                  if("pca" == dim_red){
                    Xr$rotation <- Xr$rotation[tset, ]
                  }
                  if("ica" == dim_red){
                    Xr$M <- Xr$M[tset, ]
                  }
                  
                } else {
                  Xr <- Xr[, tset]
                  nc <- NA
                }
                
                res <- lapply(seq_along(formula_list), function(j){
                  m <- ge_im(X = Xr, p = plist[[j]][tset,], 
                             formula = formula_list[[j]], 
                             method = method, dim_red = dim_red, drX = TRUE, nc = nc)
                  pred <- predict(m, newdata = plist[[j]])
                  cve <- mperf(X[, -tset], pred[, -tset], to_compute = to_compute, is.t = T)
                  mpf <- mperf(X[, tset], pred[, tset], to_compute = to_compute, is.t = T)
                  return(list(cve = cve, mpf = mpf))
                })
                cve <- do.call(rbind, lapply(res, function(elt){ unlist(elt[["cve"]]) }))
                mpf <- do.call(rbind, lapply(res, function(elt){ unlist(elt[["mpf"]]) }))
                gc()
                return(list(cve = cve, mpf = mpf))
              })
  
  if(verbose)
    message("...Cleanup and formatting")
  
  
  parallel::stopCluster(cl)
  gc(verbose = F)
  
  cve <- simplify2array(lapply(res, function(elt){ as.matrix(elt[["cve"]]) }))
  mpf <- simplify2array(lapply(res, function(elt){ as.matrix(elt[["mpf"]]) }))
  
  cve <- lapply(seq_along(formula_list), function(i){ t(cve[i,,]) })
  mpf <- lapply(seq_along(formula_list), function(i){ t(mpf[i,,]) })
  
  res <- list(formula_list = formula_list, 
              cv.n = cv.n,
              cv.s = cv.s,
              method = method,
              dim_red = dim_red,
              nc = nc,
              ... = ..., # extra model arguments
              mpf = mpf,
              cve = cve)
  class(res) <- "geimCV"
  return(res)
}


#' Plot geimCV results
#' 
#' @param x a geimCV object, as returned by \code{\link{ge_imCV}}
#' @param join.plots boolean ; if TRUE, sets up par(mfrow) appropriately
#' @param to_plot indices to plot, defaults to all computed indices.
#' @param names the labels of the models, defaults to model formulas.
#' @param names.arrange if defined, prints names with shifted height (like 'steps') by given number of labels
#' @param tcol text color for names (only when \code{names.arrange} is defined)
#' @param swarm boolean ; if TRUE, displays individual values on top of boxplots as swarms
#' @param swarmargs list of parameters given to \code{\link[beeswarm]{beeswarm}}
#' @param main title of the plots, pasted with indices on each plot.
#' @param ... extra arguments passed on to \code{\link[graphics]{boxplot}}
#' 
#' @export
#' 
#' @return invisibly returns a list of boxplot objects that can be re-plotted through \code{\link[graphics]{bxp}}.
#' 
#' @eval interpol_example()
#' 
#' @importFrom beeswarm beeswarm
#' @importFrom graphics boxplot mtext par
plot.geimCV <- function(x, join.plots = TRUE, 
                        to_plot = colnames(x$cve[[1]]), 
                        names = as.character(x$formula_list),
                        names.arrange = NULL, tcol = 1,
                        swarm = T, swarmargs = list(pch = 16),
                        main = NULL, ...){
  if(join.plots){
    graphics::par(mfrow = c(2, length(to_plot)))
  }
  blist <- lapply(c("cve", "mpf"), function(ms){
    y <- x[[ms]]
    blist <- lapply(to_plot, function(idx){
      df <- do.call(cbind, lapply(seq_along(x$formula_list), function(i){ y[[i]][, idx]}))
      b <- graphics::boxplot(df, names = names, xaxt = ifelse(is.null(names.arrange), "s", "n"),
                             ylab = idx, main = ifelse(is.null(main), paste(ms, idx, sep = " - "),
                                                       paste0("[",ms, " - ",idx, "] ", main)),
                   ...)
      if(!is.null(names.arrange)){
        graphics::mtext(names, side = 1, at = seq_along(x$formula_list),
                        line = .5 + seq_len(names.arrange) - 1, cex = .9, col = tcol)
      }
      if(swarm){
        do.call(beeswarm::beeswarm, c(df~col(df), add = T, swarmargs))
      }
      class(b) <- "bxp"
      return(b)
    })
    names(blist) <- to_plot
    return(blist)
  })
  names(blist) <- c("cve", "mpf")
  return(invisible(blist))
}


#' Print a geimCV object
#' 
#' Prints a \code{geimCV} object
#' 
#' @param x a \code{geimCV} object, as returned by \code{\link{ge_imCV}}.
#' @param ... arguments passed on to \code{\link{print}}
#' 
#' @export
#' 
print.geimCV <- function(x, ...){
  
  nf <- length(x$formula_list)
  
  
  cat("geimCV object\n---")
  cat("\nNb. compared models: ", nf)
  cat("\nNb. CV steps:        ", x$cv.n)
  cat("\nTraining set ratio:  ", x$cv.s)
  
  
  cat("\nmethod:   ", x$method)
  if(x$method != "limma"){
    
    cat("\ndim_red:  ", x$dim_red)
    cat("\nnc:       ", x$nc)
  }
  cat("\n---\n")
  
  cat("Median CV errors:\n")
  cve <- do.call(rbind, lapply(seq_len(nf), function(i){
    apply(x$cve[[i]], 2, median)
  }))
  rownames(cve) <- paste0(as.character(x$formula_list), "  ")
  print(cve, ...)
  cat("---\n")
  cat("Median model performance:\n")
  cve <- do.call(rbind, lapply(seq_len(nf), function(i){
    apply(x$mpf[[i]], 2, median)
  }))
  rownames(cve) <- paste0(as.character(x$formula_list), "  ")
  print(cve, ...)
}
