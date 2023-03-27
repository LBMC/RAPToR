# setwd('vignettes/')



# Main vignette
rmarkdown::render(input = 'RAPToR.Rmd', output_format = "all", clean = TRUE)

# Showcase vignette
rmarkdown::render(input = 'RAPToR-showcase.Rmd', output_format = "all", clean = TRUE)

# Data pkg vignette 
rmarkdown::render(input = 'RAPToR-datapkgs.Rmd', output_format = "all", clean = TRUE)

# DE correction vignette
rmarkdown::render(input = 'RAPToR-DEcorrection.Rmd', output_format = "all", clean = TRUE)

# Refbuilding vignette
rmarkdown::render(input = 'RAPToR-refbuilding.Rmd', output_format = "all", clean = TRUE)



vlist <- c("RAPToR",
           "RAPToR-showcase",
           "RAPToR-datapkgs",
           "RAPToR-DEcorrection",
           "RAPToR-refbuilding")

# rename the pdf output for separate vignette index entry
file.rename(from =  paste0(vlist, ".pdf"), to =  paste0(vlist, "-pdf.pdf"))

# create vignette .asis files for html and pdf entries
for(v in vlist){
  vpdf <- paste0(v, "-pdf")
  txts <- paste0(
    "%\\VignetteIndexEntry{",c(v, vpdf),"}\n",
    "%\\VignetteEngine{R.rsp::asis}\n")

  cat(txts[1], file = paste0(v, ".html.asis"), append = F)
  cat(txts[2], file = paste0(vpdf, ".pdf.asis"), append = F)
}
