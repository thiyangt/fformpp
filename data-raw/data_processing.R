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

## Save forecast from M1
library(seer)
library(Mcomp)
yearlym1.test <- subset(M1, "yearly")[1:2]
accuracy_info <- fcast_accuracy(tslist=yearlym1.test,
                                models= c("ets","arima","rw","rwd","wn","theta","nn"),
                                database ="M1", cal_MASE, h=6, length_out = 1, fcast_save = TRUE)
fcast_m1 <- accuracy_info$forecasts
save(fcast_m1, file="data-raw/fcast_m1.rda")
library(devtools)
use_data("fcast_m1")



aa <- yearlym1.test[1:2]

pp <- predict.m1[1:2,]
pp[1,] <- rep(1,7)

min.fcasterror <- individual_forecast(predicted= pp,
                                      accmat=cal_MASE,
                                      real.error=accuracy_info$accuracy,
                                      tslist=aa,
                                      forecast_list=accuracy_info$forecasts, h=6)

predicted=pp
accmat=cal_MASE
real.error=accuracy_info$accuracy
tslist=aa
forecast_list=accuracy_info$forecasts
