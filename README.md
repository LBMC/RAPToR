# wormAge R package


This package aims to recover the developmental stage of C. elegans worms based on their gene expression profiles.
Reference time series data is used for this purpose, but one can create their own reference from a time series.

## Installation

To install the package, download the repo, open a terminal in the repo folder and type
```
R CMD INSTALL .
```

You may need to install the `GEOquery` package separately (see (https://bioconductor.org/packages/release/bioc/html/GEOquery.html)[Bioconductor website], or (https://bioconductor.org/packages/3.4/bioc/html/GEOquery.html)[here], if you happen to have R<3.5).


