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



#' Get a GPL200 dataset with WBGene ids
#' 
#' This function either downloads a GPL200 GEO dataset from given accession number and transfers 
#' it to WBGene ids, or takes input Gene expression matrix \emph{of GPL200 format} to do the same thing.
#' **Warning :** Since not all Affymetrix array probes (GPL200 format) have a WBGene equivalent,
#' about 30% of the initial 22625 genes are lost.
#' 
#' @param GEO_id a GEO accession number (e.g 'GSE52747'), must be of GPL200 platform format.
#' @param expr.matrix a raw gene expression matrix of of GPL200 platform format. Is ignored if a GEO accession number is given.
#' @param raw wether to return the raw data or perform log(expr.matrix +1) before output.
#' 
#' @return an expression matrix with WBGene IDs as row names
#' 
#' @export
#' 
#' @examples
#' \donttest{ 
#' gorrepati_wb <- GPL200_to_WB(GEO_id="GSE52747")
#' head(gorrepati_wb)
#' }
#'
#' @importFrom GEOquery getGEO Table
#' @importFrom Biobase exprs
GPL200_to_WB <- function(GEO_id=NULL, expr.matrix=NULL, raw=F)
{
  requireNamespace("GEOquery", quietly = T)
  
  if(is.null(GEO_id) & is.null(expr.matrix)){
    stop("Either a GEO Accession number or an expression matrix must be given")
  }
  if(!is.null(GEO_id) & !is.null(expr.matrix)){
    warning("Ignoring given expr.matrix and downloading from GEO")
    rm(expr.matrix)
  }
  
  # get GEO Platform GPL200 from GEO with GEOquery
  gpl <- GEOquery::getGEO("GPL200", getGPL=FALSE)
  
  # select GPL200 ID and Gene symbols from table
  gpl_sub <-  GEOquery::Table(gpl)[,c(1,13)]
  
  # free space
  rm(gpl)
  
  # extract wb.ids from 'Gene Symbol' field 
  gpl_sub$wb.id <- sapply(gpl_sub$`Gene Symbol`, 
                          function(g){
                            s <- strsplit(g, split = ' /// ')[[1]]
                            ifelse(length(s)>1, s[2], NA)})
  
  gpl_sub <- na.omit(gpl_sub)
  
  
  if(!is.null(GEO_id)){
    gds <- getGEO(GEO_id, GSEMatrix = T, getGPL = F)
    gds <- gds[[1]]
    expr.matrix <- Biobase::exprs(gds)
    
  }

  # subset to the probes we have wb.ids for
  expr.matrix <- expr.matrix[gpl_sub$ID,] 
  
  # set the rownames to wb.ids
  rownames(expr.matrix) <- gpl_sub$wb.id
  
  if(!raw){
    # log transf.
    expr.matrix <- log(expr.matrix+1)
  }
  
  return(expr.matrix)
  
}
