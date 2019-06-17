.onLoad <- function(libname = find.package("wormAge"), pkgname = "wormAge"){
  # CRAN Note avoidance
  if(getRversion() >= "2.15.1"){
    
    utils::globalVariables(c("Cel_embryo", "Cel_larval",  "Cel_YA_adult1", 
                             #"Cel_YA_adult2", 
                             "ref_tables"))
  }
  
}