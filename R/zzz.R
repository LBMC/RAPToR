.onLoad <- function(libname = find.package("wormAge"), pkgname = "wormAge"){
  # CRAN Note avoidance
  if(getRversion() >= "2.15.1"){
    # 
    utils::globalVariables(c("ref_table"))
  }
  
}