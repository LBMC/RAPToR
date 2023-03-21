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

