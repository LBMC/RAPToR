# setwd('vignettes/')

# DE correction vignette
rmarkdown::render(input = 'RAPToR-DEcorrection.Rmd', output_format = "all", clean = TRUE)

# Showcase vignette
rmarkdown::render(input = 'RAPToR-showcase.Rmd', output_format = "all", clean = TRUE)

# Main vignette
rmarkdown::render(input = 'RAPToR.Rmd', output_format = "all", clean = TRUE)
