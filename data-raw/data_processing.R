##pkg
library(devtools)
library(tidyverse)
library(magrittr)

## read data
load("data-raw/classlabelM1Y.rda")
load("data-raw/featuresM1Y.rda")

## include data to the package
##classlabel
clm1y <- classlabelM1Y$accuracy
forecast.error <- clm1y[seq(1, nrow(clm1y), by = 2), ]
forecast.error <- data.frame(forecast.error)
use_data(forecast.error)
##features
featuresM1Y <- featuresM1Y %>%
                        select(trend, ur_pp, spikiness, beta, diff1y_acf1,
                               linearity, curvature, N)
use_data(featuresM1Y)
