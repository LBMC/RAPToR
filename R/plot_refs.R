#' Plot reference timelines
#' 
#' Plots a timeline chart for references of the chosen organism
#' 
#' @param organism the organism name; can be abbreviated.
#' 
#' @export
#' 
#' @examples 
#' \donttest{
#' plot_refs("Caenorhabditis")
#' }
#' 
#' @importFrom utils data
plot_refs <- function(organism){
  utils::data("ref_table", envir = environment())  
  organism <- match.arg(arg = organism, choices = unique(ref_table$organism))
  
  dpkg <- ref_table[which(ref_table$organism == organism)[1], "data_pkg"]
  # check if needed data package is loaded
  if(!requireNamespace(dpkg, quietly = T)){
    stop(paste0("You must install the ", dpkg, " data package to plot the ", organism, " timelines"))
  }
  
  o_code <- ref_table[which(ref_table$organism == organism)[1], "name"]
  o_code <- strsplit(o_code, '_')[[1]][1]
  
  plot_func <- do.call(`::`, list(dpkg, paste0('plot_refs_', o_code)))
                        
  return(plot_func())
}
