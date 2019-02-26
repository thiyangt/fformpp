
[![Project Status: Active ? The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active) [![Build Status](https://travis-ci.org/thiyangt/fformpp.svg?branch=master)](https://travis-ci.org/thiyangt/fformpp.svg?branch=masterr)

[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)

<!-- README.md is generated from README.Rmd. Please edit that file -->
fformpp
=======

Installation
------------

The linked packages [flutils](https://github.com/feng-li/flutils) and [movingknots](https://github.com/feng-li/movingknots) should be installed prior to the installation of fformpp.

``` r
# install.packages("devtools")
devtools::install_github("thiyangt/fformpp")
library(fformpp)
```

Usage
-----

**Load packages**

``` r
library(methods)
library(MASS)
library(Matrix)
library(mvtnorm)
library(flutils)
library(fformpp)
```

**Load example dataset**

``` r
data(features.df)
data(forecast.error)
features_mat <- as.matrix(features.df)
forecast.error <- as.matrix(forecast.error)
```

**Fit surface regression model**

``` r
## This will take time.
n <- dim(forecast.error)[1]
p <- dim(forecast.error)[2]

fformpp.model <- fit_fformpp(feamat=features_mat, accmat=forecast.error, 
                             sknots=2, aknots=2,
                            fix.s=0, fix.a=0, fix.shrinkage=1:p, fix.covariance=0,
                            fix.coefficients=0, n.iter=100,
                            knot.moving.algorithm="Random-Walk",
                            ptype=c("identity", "identity", "identity"),
                            prior.knots=n)
```

**Predict forecast error on new data**
