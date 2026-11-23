
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sptmtracks

<!-- badges: start -->

[![GitHub
issues](https://img.shields.io/github/issues/tjgu/sPTMTracks)](https://github.com/tjgu/sPTMTracks/issues)
[![GitHub
pulls](https://img.shields.io/github/issues-pr/tjgu/sPTMTracks)](https://github.com/tjgu/sPTMTracks/pulls)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![Bioc release
status](http://www.bioconductor.org/shields/build/release/bioc/sptmtracks.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/sptmtracks)
[![Bioc devel
status](http://www.bioconductor.org/shields/build/devel/bioc/sptmtracks.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/sptmtracks)
[![Bioc downloads
rank](https://bioconductor.org/shields/downloads/release/sptmtracks.svg)](http://bioconductor.org/packages/stats/bioc/sptmtracks/)
[![Bioc
support](https://bioconductor.org/shields/posts/sptmtracks.svg)](https://support.bioconductor.org/tag/sptmtracks)
[![Bioc
history](https://bioconductor.org/shields/years-in-bioc/sptmtracks.svg)](https://bioconductor.org/packages/release/bioc/html/sptmtracks.html#since)
[![Bioc last
commit](https://bioconductor.org/shields/lastcommit/devel/bioc/sptmtracks.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/sptmtracks/)
[![Bioc
dependencies](https://bioconductor.org/shields/dependencies/release/sptmtracks.svg)](https://bioconductor.org/packages/release/bioc/html/sptmtracks.html#since)
<!-- badges: end -->

Understanding how post-translational modifications (PTMs), including
glycosylation, influence protein function and disease requires
integrating three complementary dimensions: which residues are modified,
how modification levels vary across patient samples, and where modified
residues occur within protein sequence and three-dimensional structure.
Existing PTM visualization tools are often limited to pre-annotated
sites or sequence-level displays and generally do not integrate
quantitative mass spectrometry (MS) data with structural context and
patient-level clinical variables. sPTMtracks is an Rshiny dashboard
which addresses this gap by providing an interface for interactive
visualization of site-specific PTMs that integrates quantitative MS
measurements with the linear protein sequence and AlphaFold-predicted
three-dimensional structure.

## Installation instructions

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `sptmtracks` from
[Bioconductor](http://bioconductor.org/) using the following code:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("sptmtracks")
```

And the development version from
[GitHub](https://github.com/tjgu/sPTMTracks) with:

``` r
BiocManager::install("tjgu/sPTMTracks")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library("sptmtracks")
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" alt="" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub!

## Citation

Below is the citation output from using `citation('sptmtracks')` in R.
Please run this yourself to check for any updates on how to cite
**sptmtracks**.

``` r
# print(citation('sptmtracks'), bibtex = TRUE)
```

Please note that the `sptmtracks` was only made possible thanks to many
other R and bioinformatics software authors, which are cited either in
the vignettes and/or the paper(s) describing this package.

## Code of Conduct

Please note that the `sptmtracks` project is released with a
[Contributor Code of
Conduct](http://bioconductor.org/about/code-of-conduct/). By
contributing to this project, you agree to abide by its terms.

## Development tools

- Continuous code testing is possible thanks to [GitHub
  actions](https://www.tidyverse.org/blog/2020/04/usethis-1-6-0/)
  through *[usethis](https://CRAN.R-project.org/package=usethis)*,
  *[remotes](https://CRAN.R-project.org/package=remotes)*, and
  *[rcmdcheck](https://CRAN.R-project.org/package=rcmdcheck)* customized
  to use [Bioconductor’s docker
  containers](https://www.bioconductor.org/help/docker/) and
  *[BiocCheck](https://bioconductor.org/packages/3.22/BiocCheck)*.
- Code coverage assessment is possible thanks to
  [codecov](https://codecov.io/gh) and
  *[covr](https://CRAN.R-project.org/package=covr)*.
- The [documentation website](http://tjgu.github.io/sptmtracks) is
  automatically updated thanks to
  *[pkgdown](https://CRAN.R-project.org/package=pkgdown)*.
- The code is styled automatically thanks to
  *[styler](https://CRAN.R-project.org/package=styler)*.
- The documentation is formatted thanks to
  *[devtools](https://CRAN.R-project.org/package=devtools)* and
  *[roxygen2](https://CRAN.R-project.org/package=roxygen2)*.

For more details, check the `dev` directory.

This package was developed using
*[biocthis](https://bioconductor.org/packages/3.22/biocthis)*.
