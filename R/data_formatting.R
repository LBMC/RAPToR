#' Format the sample and reference data to match
#'
#' This function formats the sample and reference data to have a matching gene set.
#' This removes genes that are not present in both datasets.
#' 
#' @param samp the sample gene expression matrix with genes as rows and individuals as columns
#' @param refdata the reference matrix, same format as \code{samp}
#' @param na.rm if TRUE, rows with NA are removed
#' @param verbose if TRUE, prints a summary of initial and resulting gene counts
#'
#' @return a list with matching \code{samp} and \code{refdata} as well as \code{inter.genes}, character vector with the matching gene IDs
#' 
#' @export
#' 
format_to_ref <- function(samp, refdata, 
                          na.rm=T, verbose=T)
{
  if(all(!rownames(samp)%in%rownames(refdata))){
    stop("No matching gene IDs (rownames) between sample and reference.")
  }
  if(is.null(rownames(samp))){
    stop("Rownames of samp don't hold gene IDs.")
  }
  
  samp <- as.matrix(samp)
  mode(samp) <- 'numeric'
  
  if(na.rm){
    samp <- na.omit(samp)
    refdata <- na.omit(refdata)
  }

  
  l.r <- nrow(refdata)
  l.s <- nrow(samp)
  
  inter.genes <- intersect(rownames(refdata), rownames(samp))
  
  samp <- samp[inter.genes,, drop=F]
  refdata <- refdata[inter.genes,, drop=F]
  
  if(verbose){
    to.print <- rbind(refdata=l.r, samp=l.s, 
                      intersect.genes=length(inter.genes))
    colnames(to.print) <- "nb.genes"
    print(to.print)
  }
  
  return(list(samp=samp, refdata=refdata, inter.genes=inter.genes))
}






#' Format/convert the sample gene IDs
#'
#' This function converts the gene IDs from given input and aggregates results 
#' if necessary, for resulting IDs to be unique.
#' 
#' @param X the sample gene expression matrix with genes as rows and individuals as columns
#' @param IDs A dataframe holding current IDs for \code{X} and target IDs 
#' @param from,to defaults to 1 and 2 ; columns of \code{IDs} holding current and target IDs of \code{X} respectively
#' @param aggr.fun the function used for \code{\link[stats]{aggregate}}
#' @param verbose if TRUE, prints number of genes kept and aggregated
#'
#' @return a gene expression matrix, with new IDs as rownames
#' 
#' @export
#' 
#' @examples 
#' dat <- Cel_larval$X[1:10,1:3] # get dummy dataset
#' ids <- cbind(id.from=rownames(dat), id.to=rep(1:4, length.out=10)) # make id df
#' format_ids(dat, ids)
#' 
format_ids <- function(X, IDs, from=1, to=2, 
                       aggr.fun=mean, verbose=TRUE){
  requireNamespace("stats", quietly = T)
  n.i <- nrow(X)
  keep.ids <- intersect(rownames(X), IDs[,from])
  if(length(keep.ids)==0)
    stop("No matching IDs between rownames(X) and IDs[,from]")
  
  X.a <- stats::aggregate(X[keep.ids,],
                          by=list(as.factor(IDs[match(keep.ids, IDs[,from]),to])),
                          aggr.fun)
  X <- X.a[,-1]
  rownames(X) <- as.character(X.a[,1])
  n.f <- nrow(X)
  if(verbose){
    message(paste("Kept", length(keep.ids), "out of", n.i, "- aggregated into", n.f))
  }
  return(X)
}

