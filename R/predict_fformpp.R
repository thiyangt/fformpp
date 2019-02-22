#' Function to predict forecast error measures from the
#' the fitted fformpp
#'
#' This function is written based on the example codes in movingknots
#' package by Dr. Feng Li. Please see https://github.com/feng-li/movingknots for more details.
#'
#' @param model fitted fformpp models, output of the function fit_fformpp
#' @param feature.df feature matrix of test data
#' @param model.names vector of names of the forecast algorithms, similar to the order
#' pf accmat argument in fit_ebmsr
#' @param real.MASE optional, matrix of MASE values for all algorithms for test data
#' @param log if log transformation is used to convert Y values to real line
#' @importFrom magrittr %>%
#' @importFrom stats median
#' @export
predict_fformpp <- function(model, feature.df, model.names, real.MASE, log=TRUE){

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
      knots.ilst <- knots_mat2list(model$out.fitted[["Params"]][["knots"]][, , i, iCross], model$spline.args)
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

  if(real.MASE=="NA"){
    return(pred.mean)
    } else { # if you have the matrix that calculated pedictions from all models

  ## Comparison with real MASE
  fm <- apply(pred.mean, 1, which.min)
  MASEOpt <- rep(NA, length(fm))
  for (i in 1: length(fm))
  {
    MASEOpt[i] <- real.MASE[i, fm[i]]
  }

  summary_prediction_results <- rbind(c(mean(MASEOpt), apply(real.MASE, 2, mean)),
                                      c(median(MASEOpt), apply(real.MASE, 2, median)))
  dimnames(summary_prediction_results) <- list(c("mean", "median"), c("our_method", colnames(real.MASE)))

  return(list(predictions=pred.mean,minmase=MASEOpt, summary=summary_prediction_results))
    }
}
#' @examples
#' predictions <- predict_fformpp(fformpp, feamat, c("arima","ets","rw","rwd", "theta", "nn"), real.MASE=accmat, log=FALSE)