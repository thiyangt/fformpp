#' @importFrom tibble tibble
NULL
#'features.df
#'
#' @description features calculated on M3 competition data
#' @format A data frame with 645 rows and 25 variables
#' @examples
#' data(features.df)
#' head(features.df)
"features.df"

#' forecast.error
#'
#' @desciption MASE values calculated based on six forecast algorithms
#' @format A data frame with 645 rows with 7 variables
#' \describe{
#' \item{ets}{ets algorithm in the forecast package}
#' \item{arima}{auto.arima algorithm in the forecast package}
#' \item{rw}{random walk}
#' \item{rwd}{random walk with drift}
#' \item{wn}{white noise process}
#' \item{theta}{theta approach}
#' \item{nn}{neural network approach}
#' }
"forecast.error"

#' fformpp.model
#' @description fitted model from fit_fformpp
#' @format list containing fitted model and spline arguments
"fformpp.model"

#' features.df.m1
#' @description features calculated on M1 competition data
#' @format dataframe features calculated on M1 competition data
"features.df.m1"

#' forecast.error.m1
#' @description forecast error calculated on M1 data
#' @format dataframe features calculated on M1 competition data
"forecast.error.m1"

#' fcast_m1
#' @description forecasts of M1 data
#' @format list
"fcast_m1"
