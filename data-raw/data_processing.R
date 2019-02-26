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
