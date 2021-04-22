# GenomicOriginsScripts <img src="man/figures/logo.png" align="right" alt="" width="120" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/k-hench/GenomicOriginsScripts/workflows/R-CMD-check/badge.svg)](https://github.com/k-hench/GenomicOriginsScripts/actions)
[![DOI](https://zenodo.org/badge/208829675.svg)](https://zenodo.org/badge/latestdoi/208829675)
<!-- badges: end -->

The **R** package **GenomicOriginsScripts** provides the scripts needed to reproduce the figures and the supplementary figures of the paper *Ancestral variation, hybridization and modularity fuel a marine radiation* by Hench, Helmkampf, McMillan and Puebla.

## Dependencies 

**GenomicOriginsScripts** depends on several non-CRAN R-packages.
To be able to install the package successfully, the following packages will also need to be installed:

```r
install.packages("remotes")
remotes::install_bioc("rtracklayer")
remotes::install_github("YuLab-SMU/ggtree")
remotes::install_github("k-hench/hypogen")
remotes::install_github("k-hench/hypoimg")
```

## Install

To install **GenomicOriginsScripts** please run:

```r
remotes::install_github("k-hench/GenomicOriginsScripts")
```

This package is intended as a supplement to the [code acompanying the publication](https://k-hench.github.io/hamlet_radiation/) and is used extensively in **R** scripts that produce the figures of that publication.

---
