#' Gene ID table
#'
#' A dataframe with WBGene IDs, common public name and Sequence name.
#' This table may be especially useful to convert a dataset's row ids to 
#' WBGene IDs, as needed to find the overlapping geneset with the reference
#' and perform the age estimation.
#'
#' @docType data
#'
#' @usage data(gene_ids)
#'
#' @format a dataframe with corresponding WBGene IDs (\code{WormBase.Gene.ID}), public names (\code{Public.Name}) and sequence names (\code{Sequence.Name}).
#'
#' @keywords datasets
#'
#'
#' @source \href{https://www.wormbase.org/}{WormBase}
#'
#' @examples
#' data(gene_ids)
#' head(gene_ids)
"gene_ids"