# wormAge R package


This package aims to recover the developmental stage of C. elegans worms based on their gene expression profiles.
Reference time series data is provided from the litterature for this purpose, but one can create their own reference from a time series.

## Installation

To install the package, you can use the `devtools` R package. This should be done in your R console :
```
> library(devtools)
> devtools::install_github("LBMC/wormAge", build_vignettes = T)
```

You can also install the package, by downloading the repo, opening a terminal in the repo folder and typing :
```
R CMD INSTALL .
```

You may need to install the `GEOquery` and `limma` packages separately as they are not on CRAN (see [Bioconductor website](https://bioconductor.org/packages/)).

Also, note that the vignette may take a couple of minutes to build.

## Update info

### v0.5

 - Added a reference dataset for young adult to adult worms of better quality than the Reinke dataset (Sterken)
 - Restructured the `ref_tables` object
 - Added a `plot_ref_timelines()` function to plot the datasets' coverage and map key developmental stages
 - Joined both previous vignettes into one `wormAge` vignette and expanded the usage information (still in progress)
 - Updated the `ae` object to include a `$call` output with its call for reproducibility
 - Changed the `show.init_estimates` and `col.i` parameters to `show.prior` and `col.p` in the `ae` plot functions

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
 