
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

### Following example illustrates how package functionalities work

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
## This will take time. This model is saved in the package.
n <- dim(forecast.error)[1]
p <- dim(forecast.error)[2]

fformpp.model <- fit_fformpp(feamat=features_mat, accmat=forecast.error, 
                             sknots=2, aknots=2,
                            fix.s=0, fix.a=0, fix.shrinkage=1:p,            fix.covariance=0,
                            fix.coefficients=0, n.iter=100,
                            knot.moving.algorithm="Random-Walk",
                            ptype=c("identity", "identity", "identity"),
                            prior.knots=n)

```

**Predict forecast error on new data**

``` r
data("fformpp.model")
data("forecast.error.m1")
data("features.df.m1")
predict.m1 <- predict_fformpp(fformpp.model, features.df.m1, c("ets", "arima", "rw", "rwd", "wn", "theta", "nn"), log=FALSE)
head(predict.m1)
#>           ets    arima       rw      rwd        wn    theta       nn
#> [1,] 5.067975 5.084571 5.280992 4.269131 16.842958 4.423122 4.683714
#> [2,] 2.147037 1.969738 1.804267 2.071608  7.439973 1.419295 3.000198
#> [3,] 3.935932 3.386601 4.048462 3.912119 12.918865 3.477314 2.834796
#> [4,] 2.204292 3.186613 2.232075 2.487253  3.292003 2.190116 1.770282
#> [5,] 5.239636 4.177136 4.646175 4.225863 15.408236 4.102705 5.800446
#> [6,] 4.393721 4.688228 5.206441 3.583013 15.980581 4.161738 3.784316
```
