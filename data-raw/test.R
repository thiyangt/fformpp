# Model specifications: with 5000 iterations, without making transformations to features
#----------------------------------------------------------------------------------------
## ---- setup
rm(list = ls())
gc()

## ---- pkg
library(tidyverse)
library(methods)
library(MASS)
library(Matrix)
library(mvtnorm)
library(flutils)
library(movingknots)
library(fformpp)

## ---- data
load("Beijing/features.df.rda")
load("Beijing/forecast.error.rda")
features_mat <- as.matrix(features.df)
forecast.error <- as.matrix(forecast.error)

## ----fit surface regression models
n <- dim(forecast.error)[1]
p <- dim(forecast.error)[2] 

fformpp.model <- fit_fformpp(feamat=features_mat,
                                accmat=forecast.error, 
                                sknots = 2, aknots = 2,  # arguments for surface and additive splines
                                fix.s = 0, fix.a = 0,  # fix parameters, 0 means all are updated
                                fix.shrinkage=1:p, # shrinkages for covaiates are not updated
                                fix.covariance = 0, 
                                fix.coefficients = 0,
                                n.iter = 100,
                                knot.moving.algorithm = "Random-Walk",
                                ptype = c("identity", "identity", "identity"), 
                                prior.knots=n)
save(fformpp.model, file="Beijing/fformpp.model.rda")