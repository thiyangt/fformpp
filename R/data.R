#' @importFrom tibble tibble
NULL
#'featuresM1Y
#'
#' @description Six time series features calculated on M1 competition data
#' @format A data frame with 181 rows and 8 variables
#' \describe{
#' \item{trend}{strength of trend}
#' \item{ur_pp}{test statistic of Phillip-Perron test }
#' \item{spikiness}{Spikiness}
#' \item{beta}{parameter estimate of beta in ETS(A,A,N)}
#' \item{diff1y_acf1}{first ACF value of the differenced series}
#' \item{linearity}{strength of linearity in a time series}
#' \item{curvature}{strength of curvature}
#' \item{N}{length of the series}
#' }
#' @examples
#' data(featuresM1Y)
#' head(featuresM1Y)
"featuresM1Y"

#' forecast.error
#'
#' @desciption MASE values calculated based on six forecast algorithms
#' @format A data frame with 181 rows with 7 variables
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
