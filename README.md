# `RAPToR` R package

`RAPToR` (**R**eal **A**ge **P**rediction from **T**ranscriptome staging **o**n **R**eference) is a tool to accurately predict the developmental age of individual samples from their gene expression profiles. 

We achieve this by building high-temporal-resolution time series across the development of organisms from multiple available datasets, which we use as a reference to stage samples. 
Inferred age can then be used in multiple ways to 
precisely estimate perturbations effects on developmental timing, 
increase power in differential expression analyses, 
estimate differential expression due to uncontrolled development and, 
most importantly, to recover perturbation specific effects on gene expression even in the extreme scenario when the perturbation is completely confounded by development. 

Please cite our preprint if you use RAPToR in your research:

 - Bulteau R., Francesconi M. Real Age Prediction from the Transcriptome with RAPToR (2021) *bioRxiv* doi: (https://doi.org/10.1101/2021.09.07.459270)

## Installation

To install the package, you can use the `devtools` R package. This should be done in your R console :

```r
library(devtools)
devtools::install_github("LBMC/RAPToR", build_vignettes = T)
```

If you don't have `devtools` installed, you can do the following :
```r
install.packages("devtools")
```





## Getting started

Everything you need to know to make this work is detailed in the package's main vignette. You can access it from your R console with

```r
library(RAPToR)
vignette("RAPToR")
```

### How does it work ?

The method works in a 2-step process. 

 1. From a reference gene expression time series (several of which are included in associated data-packages), a near-continous, **high-temporal-resolution reference** is built.
 1. A **correlation profile** of your samples against this reference is computed from the gene expression information, the peak of which corresponds to the estimated age.  **Bootstrapping on genes** then gives a confidence interval of the estimates.

<center>
<img src="inst/cdoc/tool_overview.png" alt="tool_overview" width="60%"/>
</center>

### What data can be used ?
The `RAPToR` package allows you to estimate the developmental age of individual samples from their *gene expression profiles*.
This means that any method outputting information on gene expression on a large scale is appropriate : RNA-seq (preferably TPM), MicroArray...

**Data must not be gene-centered**, as this destroys the relationship between gene levels within a sample.


## Current available data-packages

Data-packages hold pre-built references for quick & easy usage.

 - [`wormRef` Nematode references](https://www.github.com/LBMC/wormRef) (*C. elegans* development)
 - [`drosoRef` Drosophila references](https://www.github.com/LBMC/drosoRef) (*D. melanogaster* embryo development)
 - [`zebraRef` Zebrafish references](https://www.github.com/LBMC/zebraRef) (*D. rerio* embryo and larval development)
 - [`mouseRef` Mouse references](https://www.github.com/LBMC/mouseRef) (*M. musculus* embryo development)


<br>
<br>
<hr>

## Update info

### v1.1
#### v1.1.4
 - Fixed minor bugs in `plot.ae` 
 - Added links to new data-packages for drosophila, zebrafish and mouse (above).
 
#### v1.1.3 (warning : model construction behavior altered)
 - Fixed a bug with the PCA centering not working as intended (in practice, this removes an erroneous 1st component corresponding to mean gene expression level rather than dynamics). Applied the same changes to ICA. Reconstructed expression matrix is now comparable to input matrix (previously comparable to scaled input matrix). This change impacted the references of the `wormRef` package as well, which has been updated accordingly.
 - Updated all vignettes according to the changes of the above point. 
 - Changed dataset names to be meaningful in the main and refbuilding vignettes (eg. `dshendriks2014` instead of `ds2`). 
 - Detailed some sections of the showcase vignette.
 - Fixed a bug when specifying a color in `plot.ae()` along with `glob.above = TRUE`.
 

#### v1.1.2
 - Created a "showcase" vignette with 3 different uses of `RAPToR` `vignette("RAPToR-showcase")`
 - Updated vignettes with a `data_folder` variable to build the example objects more easily
 - Added code to generate diagnostics plots in the refuilding vignette
 - Cleared up some points of the overall vignette based on user feedback
 - Added print methods for `geim` and `geimCV` objects
 - Added the `scale(X)` within `ge_imCV()` function when using `dim_red`
 - Fixed bug when specifying `prior` but not `prior.params` in `ae()`
 - Cleaned up deprecated elements :
   - Fully deprecated (deleted) `interpol_refdata()` and `estimate.worm_age()`
   - Removed links to deprecated functions in doc
   - Added warnings for deprecation in plsr-linked functions
 
#### v1.1.1
 - Created a detailed vignette on reference building, with examples `vignette("RAPToR-refbuilding")`
 - Created a vignette on data-package building `vignette("RAPToR-datapkgs")`
 
#### v1.1.0
 - Further split the references from the main package : 
   - removed `ref_table` object (available data-packages are now listed in this README)
   - added `list_refs()` function to list the references in given data-package
   - `plot_refs()` now takes `datapkg` as an argument
   - `prepare_refdata()` now takes `datapkg` as an argument
 
 - Major changes to the gene expression interpolation method :
   - deprecated `plsr_interpol()` and `df_CV()`
   - added `ge_im()` and `ge_imCV()` functions, with an improved, more flexible model interface for reference building (see vignettes)
   - Started a vignette specifically on reference building (`vignette("RAPToR-refbuilding")`)
 
 - Updated general vignette with a broad overview on reference building. 
 
### v1.0 Major Update
 - Changed the name of the package from `wormAge` to `RAPToR`. (As the tool is not limited to *C.elegans*)
 
 - Split the references from the main tool. 
 The references are now in seperate data-packages (*e.g.* `wormRef` will need to be installed separately to access *C. elegans* pre-built references). 
 All available references and their associated data-package can be seen via the `ref_table` object.

 - Changed the gene expression interpolation method. Now using Partial Least Square Regression (PLSR), with the `plsr_interpol()` function. `interpol_refdata()` is now deprecated as a consequence.
 - Added the `df_CV()` function to help estimate the parameters for PLSR interpolation when building references
 
 - Deprecated `estimate.worm_age()`, renamed to `ae()`
 - Updated the vignette with new examples, incorporating all changes to the package. 
 

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
 
