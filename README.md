
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
## This will take time. This model is saved in the package. THe fitted model is  saved into the package for easy references.
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

**Generate forecast from the model with minimum forecast error**

``` r
library(Mcomp)
#> Loading required package: forecast
yearlym1 <- subset(M1, "yearly")
min.fcasterror <- individual_forecast(predicted=predict.m1, 
                                      accmat=cal_MASE, 
                                      real.error=forecast.error.m1, 
                                      tslist=yearlym1)
min.fcasterror
#> $models
#>   [1] 4 6 7 7 6 4 2 2 1 4 7 4 2 2 1 4 4 6 2 4 2 2 7 5 5 4 4 2 4 5 5 2 4 2 5
#>  [36] 5 5 4 6 4 4 6 6 4 6 4 6 4 1 4 7 4 5 4 2 2 5 4 2 5 2 3 4 4 2 2 2 4 1 4
#>  [71] 4 4 4 4 2 7 3 6 4 5 4 1 3 6 1 4 2 2 4 2 3 7 4 5 4 4 4 4 3 6 6 1 7 2 2
#> [106] 2 2 2 4 4 4 7 2 2 1 2 4 1 2 1 1 1 2 4 6 2 4 5 2 4 4 2 1 2 7 5 7 3 2 6
#> [141] 4 1 2 4 2 2 7 1 4 2 3 7 2 2 2 2 2 2 2 2 2 1 1 3 5 6 2 2 3 2 4 7 1 1 1
#> [176] 7 1 6 4 5 5
#> 
#> $minmase
#>   [1] 10.52761210  6.22546310 10.78527697  5.61344769 11.95670159
#>   [6]  5.69921169  7.60049192  6.84308990  9.25760811  1.23793267
#>  [11]  1.90898198 10.47077133  8.49703094  4.36913758  9.08427659
#>  [16]  1.20990669  1.61439275  4.07341544  6.05683571  2.23599924
#>  [21]  3.31794454  1.47552380  7.17717193  1.80107203  2.61004398
#>  [26]  1.15436877  1.25614120  1.21864476  0.53540100  1.81680195
#>  [31]  2.63421268  1.93192037  1.81862167  6.40541465  3.00446864
#>  [36]  7.38636819  1.91196768  0.61164722  3.74870621  2.69644785
#>  [41]  5.10341668 10.63481352  6.63412267  1.91040248  3.16719783
#>  [46]  8.57037535 53.73586632  1.64307480  0.75496556  1.79037148
#>  [51]  2.25905723  2.80014760  1.59886537  1.82402405  0.25387047
#>  [56]  2.62445201  1.22689951  0.84038360  1.84706443  2.67198877
#>  [61]  4.17243587  1.08453387  1.92452209  3.38475037  1.56600461
#>  [66]  5.50159241  2.28432579  1.53760361  1.74543256  2.12268809
#>  [71] 12.65958968  3.42245954  0.87133384  1.69685836  4.07161189
#>  [76] 44.10931880  4.21540729  6.63412267 12.74831543  1.01788772
#>  [81]  7.80688378  3.37963810  3.14644538  2.62920810  2.60940689
#>  [86]  2.73203531  2.70743534  1.93668578  2.89743590  0.82972321
#>  [91]  4.90137993  1.77773793  2.54933692  5.20886397  3.84657661
#>  [96]  1.86886761  2.18306011  2.24875931  0.66914498  3.94785422
#> [101]  2.47925480  3.77163176  4.68042455  1.57899998  4.55048752
#> [106]  0.97375938  2.84432185  3.90624211  1.34031353  5.37710665
#> [111]  7.10430933  2.54130598  2.99541982  0.59876543  0.60126910
#> [116]  2.08003675  1.87950000  0.65390173  3.09992000  0.88503436
#> [121]  0.72735092  1.54732972  1.40139123  1.64788713  4.31723596
#> [126]  1.09723495  2.93115557  0.68603491  1.86323404  4.76832406
#> [131]  2.10952495  2.58773523  7.97419436  1.82866492  7.52213857
#> [136]  4.16099770  6.19047049  0.90349483  5.26498644  3.67682938
#> [141]  0.84371342  0.09952735  1.95227324  4.60402154  1.80620861
#> [146]  2.51906507  5.10632031  7.35639284  2.29253704  3.47254800
#> [151]  3.14646007  5.56015396  4.16840894  0.42340177  1.74576093
#> [156]  1.01524242  0.97887135  0.93418975  0.95015958  1.60601348
#> [161]  1.79944413  3.63726892  0.94672217  1.02927451  1.32050710
#> [166]  8.07552668  0.41025641  8.25973952  6.51058974  5.55066667
#> [171]  5.05639113  1.64134710  1.55287633  1.60152285  1.17437947
#> [176]  1.57217926  1.20168108  2.70161792  4.25364412  1.68330558
#> [181]  4.38261326
#> 
#> $summary
#>        our_method      ets    arima       rw      rwd        wn    theta
#> mean     3.869669 3.771245 3.473598 4.893131 3.489743 10.006127 4.189472
#> median   2.549337 2.323728 2.191300 3.771522 2.292537  9.391286 3.154610
#>              nn
#> mean   4.602657
#> median 3.111168
```
