
# mtest

<!-- badges: start -->
<!-- badges: end -->

The goal of mtest is to study contingency tables. It is related to Fisher's 
exact test and to Barnard's exact test. The main difference with Barnard's test
is that the nuisance parameter is not estimated, but integrated over all 
possible values. 

## Installation

You can install the released version of mtest from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("mtest")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(mtest)
m.test(list(c(10, 7), c(11, 6)))
tailed.m.test(list(c(10, 7), c(11, 6)))
```

