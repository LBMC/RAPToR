# wormAge R package

`wormAge` was developped to mediate the issue of developmental differences between samples by estimating the age of samples from their gene expression profiles. 
This is a major problem in the field of *C. elegans*, where many factors can impact development speed.


## Installation

To install the package, you can use the `devtools` R package. This should be done in your R console :
```r
library(devtools)
devtools::install_github("LBMC/wormAge", build_vignettes = T)
```


You may need to install the `GEOquery` and `limma` packages separately as they are not on CRAN (see [Bioconductor website](https://bioconductor.org/packages/)).
This can be done with
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("limma")
```
and
```r
BiocManager::install("GEOquery")
```


## Getting started

Everything you need to know to make this work is detailed in the package's vignette. You can access it from your R console with

```r
vignette("wormAge")
```

### How does it work ?

The method works in a three-step process. 

 1. From a reference gene expression time series (several of which are included in this package), a near-continous, **high-temporal-resolution reference** is built.
 1. A **correlation profile** of your samples against this reference is computed from the gene expression information, the peak of which corresponds to the estimated age.
 1. A **bootstrap procedure** of the previous step on random subsets of genes is performed to give a confidence interval of the estimates.



### What data can be used ?
The `wormAge` package allows you to estimate the developmental age of individual samples from their *gene expression profiles*.
This means that any data providing information on gene expression on a large scale is appropriate : RNA-seq counts (or RPKM), MicroArray, Chips...

**The data must not be gene-centered**, as this destroys the relationship between gene levels within a sample.


## Update info



### v0.7
#### v0.7.2
 - Switched to the Median Absolute Deviation (MAD) instead of 95% interval for the estimate confidence interval. (As a side effect, removed IC imbalance and updated vignette)

#### v0.7.1
 - Fixed a namespace issue
 - Fixed windows compatibility for parallel computing
 
#### v0.7.0 
 - Added a `format_ids()` function to handle transformation from one ID set to another with aggregation of data
 - Updated the post-hatching series with clean, reproducible versions (along with their optimal interpolation parameters in `prepare_refdata()`)
 - `Cel_YA_adult1` reference is now built from joint datasets by Hendriks et al. (2014) and Sterken (2014)
 - Removed the `Cel_YA_adult2` reference (quality too low and reproducibility of data handling was mediocre at best)
 - Fixed a bug in the `ae` plot function with groups and subsets


### v0.6

 - Changed the reference dataset names to be more transparent
 - Added a warning in `estimate_wormage()` for *edge-of-reference* estimates
 - Continued vignette for general usage (still in progress)

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
 