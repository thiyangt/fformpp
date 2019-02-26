#' @importFrom tibble tibble
NULL
#'features.df
#'
#' @description Six time series features calculated on M1 competition data
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
"forecast.error"
