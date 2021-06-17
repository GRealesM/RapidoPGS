# RápidoPGS <a href='https://grealesm.github.io/RapidoPGS/'><img src='man/figures/logo.png' align="right" height="190.5" /></a>

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/RapidoPGS)](https://cran.r-project.org/package=RapidoPGS)
<!-- badges: end -->

**A *rápido* and lightweight method to compute Polygenic Risk Scores.**

**Last update:** 2021-06-09

**Current version:** 2.1.0.9003

This package allows to quickly (*rápido* is Spanish for "fast") compute polygenic scores (PGS) from case-control or quantitative trait GWAS summary statistic datasets, without the need of an external validation dataset.

# Background

You can find a description of the ideas behind RápidoPGS, as well as technical details in our preprint:

> [Reales G, Vigorito E, Kelemen M, & Wallace C (2021) RápidoPGS: A rapid polygenic score calculator for summary GWAS data without validation dataset. *BioRxiv*.](https://www.biorxiv.org/content/10.1101/2020.07.24.220392v2)

# News

* In version 2.1.0 we added a functionality to `rapidopgs_multi()`, which now allows users to use their own LD matrices instead of computing them on the go from a reference panel. For European datasets, we recommend downloading UK Biobank LD matrices kindly provided by Privé et al., which can be accessed [here](https://figshare.com/articles/dataset/European_LD_reference/13034123).


# Installation

RápidoPGS (1.0.2) is now available on CRAN. You can install it by typing the code below.
```
install.packages("RapidoPGS")
```
Since we are constantly improving the package, and CRAN version is currently outdated, we recommend to install the development version instead.

## Development version

There's also a development version, that can be installed from GitHub.
```
library(remotes)
install_github('GRealesM/RapidoPGS')
```

### A note on dependencies

RápidoPGS has some dependencies that aren't available directly from CRAN, so must be installed a bit differently.

**bigsnpr**

Current `bigsnpr` development version (v1.6.7 as of 2021-03-08) is available from GitHub:
```
remotes::install_github('privefl/bigsnpr')
```


**coloc**

We used a development branch of `coloc` package, which can be installed by typing:
```
remotes::install_github('chr1swallace/coloc', ref="susie")
```


**GenomicRanges**

`GenomicRanges` package is a Bioconductor package. Please type:
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GenomicRanges")
```


# Documentation

Full documentation and vignettes are available on the website (click on the cat if you're at the GitHub repo).

At the moment, vignettes cover `rapidopgs_single()` only, but we'll try to add a tutorial for `rapidopgs_multi()` soon.

