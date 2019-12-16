.onLoad <- function(libname = find.package("RAPToR"), pkgname = "RAPToR"){
  # CRAN Note avoidance
  if(getRversion() >= "2.15.1"){
    # 
    utils::globalVariables(c("ref_table"))
  }
  
}