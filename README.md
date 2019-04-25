# wormAge R package


This package aims to recover the developmental stage of C. elegans worms based on their gene expression profiles.
Reference time series data is provided from the litterature for this purpose, but one can create their own reference from a time series.

## Installation

To install the package, download the repo, open a terminal in the repo folder and type
```
R CMD INSTALL pkg/
```

Or, using the `devtools` R package, you can also do this in R (with the correct path) :
```
> library(devtools)
> devtools::install("/path_to_downloaded_repo/pkg/")
```

You may need to install the `GEOquery` package separately (see [Bioconductor website](https://bioconductor.org/packages/release/bioc/html/GEOquery.html), or [here](https://bioconductor.org/packages/3.4/bioc/html/GEOquery.html), if you happen to have R<3.5).


## Update info

### v0.4

 - Added a `summary()` and `print()` default functions for `ae` objects
 - Added vignette `prepare_data` on what data to use and example
 - Added vignette `estimate_age` on how to perform a simple estimate on example data

### v0.3
 
 - Added *IC imbalance* information on age estimates to help determine if the estimate 'jumps' between peaks during bootstrap
 - Added possibility to show bootstrap estimates on the `ae` plot function and subset the samples, as in `plot_cor.ae()`
 - Added a reference dataset for young adult to adult worms (Reinke *et al*)
 - Cleaned up all devtools warnings
 - Removed the `corg` object and its plotting function (unused)

### v0.2

 - Made the `cor.gene_expr()` function more efficient
 - Updated default plotting functions for `ae` objects
 - Included `format_to_ref()` calling in the `estimate.worm_age()` function
 - Made prior optional for age estimation
 - Included reference time tables for developmental stages in reference datasets

### v0.1

 - Creation of the package
 - Made the `interpol_refdata()` function
 - Made the `estimate.worm_age()` function and associated plotting functions
 - Included reference datasets for *C. elegans* Embryonic and larval development
 