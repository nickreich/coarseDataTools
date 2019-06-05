coarseDataTools
===============

<!-- badges: start -->
[![CRAN](http://cranlogs.r-pkg.org/badges/coarseDataTools)](https://cran.r-project.org/package=coarseDataTools)
[![Travis build status](https://travis-ci.org/nickreich/coarseDataTools.svg?branch=master)](https://travis-ci.org/nickreich/coarseDataTools)
<!-- badges: end -->

This is the repository for the coarseDataTools R package. We use this as a development space. The most recent, stable version of the package can be downloaded either from CRAN (https://cran.r-project.org/package=coarseDataTools) or from the most recent release version on github (https://github.com/nickreich/coarseDataTools/releases).

This package contains functions to analyze coarsely observed data.
    Specifically, it contains functions to (1) fit parametric accelerated
    failure time models to interval-censored survival time data, and (2)
    estimate the case-fatality ratio in scenarios with underreporting.
    This package's development was motivated by applications to infectious
    disease: in particular, problems with estimating the incubation period and
    the case fatality ratio of a given disease. Sample data files are included
    in the package.


As of March 2016, coarseDataTools imports functions from the MCMCpack package, which in turn imports functions from the graph and Rgraphviz packages. Both of these packages have been removed from CRAN, but can be installed from Bioconductor using the following code:

```
source("https://bioconductor.org/biocLite.R")
biocLite("graph")
biocLite("Rgraphviz")
```
