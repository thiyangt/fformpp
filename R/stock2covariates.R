##' Construct sensible covariates for stock market data.
##'
##' The historical finance data can be obtained from "Yahoo! Finance".
##'
##' @title Construct time series covariates.
##'
##' @param file "character".  The location of the CSV format file from Yahoo!  Finance
##' which can be the http address.
##'
##' @param g "vector".  Parameters for geometrically averaging.
##'
##' @param kappa "scaler".
##'
##' @param ma "vector" with integer entries.  Moving average indicator for the returns .
##'
##' @return "data.frame".  The returned data frame only contains the non-NA data due to
##' the moving average procedure.
##'
##' @references Geweke, and Keane (2007) p.274
##'
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##'
##' @note Created: Thu Jan 05 13:23:35 CET 2012; Current: Mon Sep 28 22:01:20 CST 2015.
##' @export
stock2covariates <- function(file, g = c(0.95, 0.80),
                             kappa = 1, ma = c(1, 5, 20))
{
  ## Read the CSV file from disk and subtract the raw covariates
  ## TODO: the comma of the csv file makes troubles.
  DataRaw <- read.csv(file, sep = ",", header=TRUE)
  nObs <- nrow(DataRaw) # Number of observations

  ## The Yahoo csv format was ordered with time descending order. Reorder it
  ## so that the first observation is from the oldest time.
  Data <- DataRaw[nObs:1, , drop = FALSE]

  ## Use adjusted close price instead of close price, for explanation see
  ## https://help.yahoo.com/kb/finance/historical-prices-sln2311.html

  Close = Data[,"Adj.Close"]
  High = Data[,"High"]
  Low = Data[, "Low"]

  ## Reserve the storage
  CloseAbs <- matrix(NA, nObs, length(g),
                     dimnames = list(NULL,
                                     paste("CloseAbs", g*100, sep = "")))
  CloseSqrt <- matrix(NA, nObs, length(g),
                      dimnames = list(NULL,
                                      paste("CloseSqrt", g*100, sep = "")))
  MaxMin <- matrix(NA, nObs, length(g),
                   dimnames = list(NULL,
                                   paste("MaxMin", g*100, sep = "")))
  RMA <- matrix(NA, nObs, length(ma),
                dimnames = list(NULL,
                                paste("RMA", ma, sep = "")))

  ## The Returns
  ## 100log(p(t)/p(t-1)), where p(t) is closing price
  Close1 <- c(NA, Close[1:(nObs-1)])
  Returns <- 100*(log(Close)-log(Close1))

  ## Time indices
  tInit <- 2 # The initial point by taking away the first lag.
  tIdx <- (max(ma, 3)+2):nObs

  ## Loop over all the observations to construct the covariates. NOTE: This is the
  ## utility function which does not require very high efficiency--loops are OK at the
  ## moment.
  for(t in tIdx)
  {
    for(i in 1:length(ma))
    {
      RMA[t, i] <- sum(Returns[(t-1):(t-ma[i])])
    }
    for(j in 1:length(g))
    {
      g.s <- g[j]^((t-2-tInit):0)
      CloseAbs[t, j] <- (1-g[j])*sum(g.s*abs(Returns[tInit:(t-2)])^kappa)
      CloseSqrt[t, j] <- sqrt((1-g[j])*
                              sum(g.s*(Returns[tInit:(t-2)]^2)^kappa))
      MaxMin[t, j] <- (1-g[j])*
      sum(g.s*(100*(log(High[tInit:(t-2)])-log(Low[tInit:(t-2)])))^kappa)
    }
  }

  X <- cbind(RMA, CloseAbs, CloseSqrt, MaxMin)[tIdx, , drop = FALSE]


  Y <- matrix(Returns[tIdx])
  ##ID <- as.Date(Data[tIdx, "Date"], "%Y-%m-%d")
  ID <- as.character(Data[tIdx, "Date"])

  rownames(X) <- ID
  rownames(Y) <- ID

  out <- list(ID = ID, X = X, Y = Y)
  return(out)
}
