#' Function to predict forecast error measures from the
#' the fitted fformpp
#'
#' This function is written based on the example codes in movingknots
#' package by Dr. Feng Li. Please see https://github.com/feng-li/movingknots for more details.
#'
#' @param model fitted fformpp models, output of the function fit_fformpp
#' @param feature.df feature matrix of test data
#' @param model.names vector of names of the forecast algorithms, similar to the order
#' of accmat argument in fit_ebmsr
#' @importFrom magrittr %>%
#' @importFrom stats median
#' @export
predict_fformpp <- function(model, feature.df, model.names){

  # Preparation of the testing file
  x.testing <- feature.df %>% as.matrix()
  colnames(x.testing) <- NULL

  x.testing <- flutils::StdData(x.testing, method = "norm-0-1")[["data"]]
  ## Which cross validation fold used
  iCross <- 1

  ###----------------------------------------------------------------------------
  ### The scripts
  ###----------------------------------------------------------------------------

  ## OUT.IF.testing <- matrix(NA, nTest, length(RdataFiles))
  nTest <- nrow(x.testing)
  nDim <- dim(model$out.fitted[["Params"]][["coefficients"]])[2]
  nIter <- dim(model$out.fitted[["Params"]][["coefficients"]])[3]
  nCross <- dim(model$out.fitted[["Params"]][["coefficients"]])[4]
  Y.pred <- array(NA, c(nTest,nDim, nIter))

  for(iCross in 1:nCross)
  {
    for(i in 1:nIter)
    {
      knots.ilst <- movingknots::knots_mat2list(model$out.fitted[["Params"]][["knots"]][, , i, iCross], model$spline.args)
      ## knots.s.mat[(1+q.s*(i-1)):(i*q.s), ] <- knots.ilst[["thinplate.s"]]
      ## knots.a.mat[(1+q.a1*(i-1)):(i*q.a1), ] <- matrix(knots.ilst[["thinplate.a"]], q.a1, 2)
      X.i <- d.matrix(x.testing, knots.ilst, model$spline.args)
      B.i <- matrix(model$out.fitted[["Params"]][["coefficients"]][, , i, iCross], ,nDim)
      Y.pred[, , i] <- X.i %*% B.i # Transformed scale,  but should be OK here.
    }
  }
  pred.mean <- apply(Y.pred, c(1, 2), mean)
  colnames(pred.mean) <- model.names
  if(log==TRUE){pred.mean <- exp(pred.mean)}
  colnames(real.MASE) <- model.names

  return(pred.mean)

}
#' @examples
#' predictions <- predict_fformpp(fformpp, feamat, c("arima","ets","rw","rwd", "theta", "nn"))
