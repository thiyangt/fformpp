#' Function to fit Efficient Bayesian Surface Regression Models.
#' This function is written based on the example codes in movingknots
#' package by Dr. Feng Li.
#'
#' @param feamat matrix of features, rows corresponds to each time series,
#' columns coresponds to features
#' @param accmat matrix of forecast errors from each method,
#' rows represent time series, columns represent forecast algorithms
#' @param sknots the dimension of knots for surface, default 2
#' @param aknots no. of knots used in each covariates for the additive part, default 2
#' @param fix.s number of knots to be fixed in the surface components, 0 means all are
#' updated
#' @param fix.a number of knots to be fixed in the additive components, default is 0 which means all are updated
#' @param fix.shrinkage, number of shrinkage covariates not to be updated, defalut is, 1:p,
#' @param fix.covariance, number of knots to be fixed in the covariance, default is 0, all are updated
#' @param fix.coefficients, number of knots to be fixed in the coefficients, default is 0, all are updated
#' @param n.iter number of ierations
#' @param knot.moving.algorithm, select either "KStepNewton" or "Random-Walk", to fasten the
#' code use Random-Walk
#' @param ptype For fixing gprior, This could be c("X'X", "identity", "identity") or c("identity", "identity", "identity")
#' @param prior.knots this could be n, log(n) or 1 # to set priors for knots

