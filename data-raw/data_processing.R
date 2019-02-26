## pkg
library(devtools)
library(tidyverse)
library(magrittr)

## read data
load("data-raw/classlabelM3Y.rda")
load("data-raw/featuresM3Y.rda")

## include data to the package
## classlabel
clm3y <- classlabelM3Y$accuracy
forecast.error <- clm3y[seq(1, nrow(clm3y), by = 2), ]
forecast.error <- data.frame(forecast.error)
use_data(forecast.error)
## features
features.df <- featuresM3Y
use_data(features.df)

## Include the fitted model
load("data-raw/fformpp.model.rda")
use_data(fformpp.model)


## Include data corresponds to M1
load("data-raw/classlabelM1Y.rda")
clm1y <- classlabelM1Y$accuracy
forecast.error.m1 <- clm1y[seq(1, nrow(clm1y), by = 2), ]
forecast.error.m1 <- data.frame(forecast.error.m1)
use_data(forecast.error.m1)
load("data-raw/featuresM1Y.rda")
features.df.m1 <-featuresM1Y
use_data(features.df.m1)
