# wormAge R package


This package aims to recover the developmental stage of C. elegans worms based on their gene expression profiles.
Reference time series data is used for this purpose, but one can create their own reference from a time series.

## Installation

To install the package, download the repo, open a terminal in the repo folder and type
```
R CMD INSTALL .
```

You may need to install the `GEOquery` package separately (see [Bioconductor website](https://bioconductor.org/packages/release/bioc/html/GEOquery.html), or [here](https://bioconductor.org/packages/3.4/bioc/html/GEOquery.html), if you happen to have R<3.5).


## Update info

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
 