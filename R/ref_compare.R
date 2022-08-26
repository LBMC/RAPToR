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
#' @eval rc_example()
#' 
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
    } 
  } else {
    if("ae" != class(ae_obj)){
      stop("ae_obj must be an object of class 'ae', as returned by 'ae()'")
    }
    ae_values <- ae_obj$age.estimates[,1]
  }
  
  idx <- sapply(ae_values, function(t) which.min(abs(ref$time - t)))
  if(return.idx){
    return(idx)
  }
  return(ref$interpGE[,idx])
}

#' Compare expression changes between sample groups to a matching developmental reference.
#' 
#' \code{ref_compare} computes log2-fold changes (logFCs) between given sample groups 
#' and logFCs between matching expression profiles of a developmental reference.
#' By comparing the logFCs of the samples and those of the reference, we can
#' quantify expression changes caused by developmental differences between sample groups
#' rather than by the variable of interest.
#' 
#' We use TPM instead of raw counts to compute logFCs because TPM are more comparable 
#' across datasets, and because reference data is interpolated from and to log(TPM+1). 
#' Because of this, our logFCs may differ from those computed by differential expression (DE) 
#' analysis tools that use raw counts (e.g. \href{https://doi.org/doi:10.18129/B9.bioc.DESeq2}{DESeq2}, 
#' \href{https://doi.org/doi:10.18129/B9.bioc.edgeR}{edgeR}) 
#' 
#' 
#' 
#' @param X sample gene expression matrix (should be comparable to reference, usually log(TPM+1)).
#' @param ref a \code{ref} object (as returned by \link{\code{make_ref}})
#' @param fac factor defining sample groups to compare (e.g. wild-type & mutant), the first factor level is used as control.
#' @param ae_obj an \code{ae} object (as returned by \link{\code{ae}})
#' @param ae_values age estimate values (used if ae_obj is NULL). 
#' @param return.idx if TRUE (default) returns reference indices. If FALSE, returns corresponding reference expression matrix.
#'
#' @return an rcmp object.
#' 
#' @export
#' 
#' @eval rc_example()
#' 
#' @importFrom Rdpack reprompt
#' @importFrom stats lm
#' 
ref_compare <- function(X, ref, fac, 
                        ae_obj=NULL, ae_values=NULL){
  
  if(!is.factor(fac)| length(levels(fac)) < 2){
    stop("fac must be a factor with at least 2 levels specifying groups to compare")
  }
  if(length(fac) != ncol(X)){
    stop("ncol(X) different from length(fac)")
  }
  
  if(is.null(ae_obj)){
    if(is.null(ae_values)){
      stop("One of ae_obj or ae_values must be specified.")
    } else {
      if(!is.numeric(ae_values) | 
         any(ae_values < min(ref$time)) | 
         any(ae_values > max(ref$time))){
        stop("ae_values must be numeric and within reference span.")
      }
    }
  } else {
    if("ae" != class(ae_obj)){
      stop("ae_obj must be an object of class 'ae', as returned by 'ae()'")
    }
    ae_values <- ae_obj$age.estimates[,1]
  }
  if(length(ae_values)!=ncol(X)){
    stop("ncol(X) doesn't match the number of age estimates (ae_obj/ae_values).")
  }
  
  
  ovl <- RAPToR::format_to_ref(X, RAPToR::get_refTP(ref, ae_values = ae_values, return.idx = F), verbose = F)
  
  lm_samp <- stats::lm(log2(exp(t(ovl$samp)))~fac)
  lm_ref <-stats::lm(log2(exp(t(ovl$ref)))~fac)
  
  coefs <- list(samp = t(coef(lm_samp)), ref = t(coef(lm_ref))) 
  
  fac_stats <- list(ae_avg = tapply(ae_values, fac, mean),
                    ae_sd = tapply(ae_values, fac, sd),
                    ae_range = as.data.frame(
                      do.call(cbind, tapply(ae_values, fac, range)),
                      row.names = c("min", "max")
                      ))
  
  res <- list(coefs = coefs, expr = ovl, fac = fac)
  class(res) <- "rcmp"
  attr(res, "fac.stats") <- fac_stats
  attr(res, "ref.info") <- list(t.unit = attr(ref, "t.unit"),
                                metadata = attr(ref, "metadata"))

  return(res)
}



#' Extract sample or reference log-fold changes from rcmp object.
#' 
#' Fetches the log2-fold change (logFC) values for a given group comparison in the samples or in the reference.
#' 
#' @param rc an rcmp object, as returned by \link{\code{ref_compare}}
#' @param l,l0 sample groups to compare. \code{l0} and \code{l} defaults to the first and second levels of \code{fac} respectively.
#' @param verbose if TRUE (default), prints warnings and selected factor levels.
#'
#' @return a dataframe with sample and reference logFCs between groups.
#' 
#' @export
#' 
#' @eval rc_example()
#' 
#' @importFrom Rdpack reprompt
#' 
get_logFC <- function(rc, l = levels(rc$fac)[2], l0 = levels(rc$fac)[1], 
                      verbose=T){
  if("rcmp"!=class(rc)){
    stop("rc must be an 'rcmp' object, as returned by ref_compare.")
  }
  ll <- levels(rc$fac)
  if(!l%in%ll | !l0%in%ll){
    stop("l and l0 must be levels of the group factor (fac).")
  }
  if(verbose){
    message(paste0("Comparing ",l0, " vs. ", l))
  }
  
  if(ll[1] == l0){ # comparison with control, lm coefs are good as is.
    # avg_l0 = beta_0
    # lfc = beta_l 
    lfc <- data.frame(
      samp = rc$coefs$samp[, which(l == ll)],
      ref = rc$coefs$ref[, which(l == ll)]
    )
  } else { # comparison between 2 non-control groups
    # avg_l = beta_0 + beta_l 
    # avg_l0 = beta_0 + beta_l0
    # lfc = (beta_0 + beta_l) - (beta_0 + beta_l0) = beta_l - beta_l0
    lfc <- data.frame(
      samp = rc$coefs$samp[, which(l == ll)] - rc$coefs$samp[, which(l0 == ll)],
      ref = rc$coefs$ref[, which(l == ll)] - rc$coefs$ref[, which(l0 == ll)]
    )
  }
  return(lfc)
}


#' Print an rcmp object
#' 
#' Prints an \code{rcmp} object
#' 
#' @param x an \code{rcmp} object, as returned by \code{\link{ref_compare}}.
#' @param ... arguments passed on to \code{\link{print}}
#' 
#' @return invisibly returns a dataframe wiht the displayed info. 
#' 
#' @export
#' 
#' 
print.rcmp <- function(x, ...){
  lfcs <- lapply(levels(x$fac)[-1],  RAPToR::get_logFC, rc=x, l0=levels(x$fac)[1], verbose=F)
  rs <- unlist(lapply(lfcs, function(lfci) cor(lfci$samp, lfci$ref)))
  ae_avg <- attr(x, "fac.stats")$ae_avg
  df <- as.data.frame(cbind(#fac = levels(x$fac),,
              ref.logFC.r = c(NA, rs),
              ref.logFC.r2 = c(NA, rs^2),
              ae.avg.dif = c(NA, ae_avg[-1] - ae_avg[1])
              ))
  rownames(df) <- sapply(seq_along(levels(x$fac)), function(i){
    paste(levels(x$fac)[i], ifelse(1==i, " (ctrl, n=", " (n="), 
          table(x$fac)[i] ,")", sep="")
  })

  cat("DE comparison with reference\n---\n")
  print(round(df, 3), ...)
  cat("---\n")
  return(invisible(df))
}
