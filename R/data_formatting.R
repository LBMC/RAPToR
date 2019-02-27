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
  
  samp <- samp[inter.genes,]
  refdata <- refdata[inter.genes,]
  
  if(verbose){
    to.print <- rbind(refdata=l.r, samp=l.s, 
                      intersect.genes=length(inter.genes))
    colnames(to.print) <- "nb.genes"
    print(to.print)
  }
  
  return(list(samp=samp, refdata=refdata, inter.genes=inter.genes))
}
