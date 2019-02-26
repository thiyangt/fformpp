#' Function to predict forecast error measures from the
#' the fitted fformpp
#'
#' This function is written based on the example codes in movingknots
#' package by Dr. Feng Li. Please see https://github.com/feng-li/movingknots for more details.
#' @param predicted predicted forecast errors from predict_fformpp
#' @param real.error optional, matrix of MASE values for all algorithms for test data
#' @param accmat function to compute forecast accuracy
#' @param tslist list of time series, as in the format of Mcomp object
#' @param forecast_list true forecast from different models
#' @importFrom magrittr %>%
#' @importFrom stats median
#' @export
individual_forecast <- function(predicted, accmat=NULL, real.error=NULL, tslist=TRUE, forecast_list=NULL){

    ## Comparison with real MASE
    #fm <- apply(pred.mean, 1, which.min)
    fm <- apply(predicted, 1, function(x) which(x == min(x, na.rm = TRUE))) ## this can track multiple minimum values
    MASEOpt <- rep(NA, length(fm))
    for (i in 1: length(fm))
    {
      if (length(fm[[i]])==1){
        MASEOpt[i] <- real.error[i, fm[[i]]]
      } else {
        min_model <- colnames(predicted)[fm[[i]]]
        comb_fcast_components <- matrix(NA, ncol=length(min_model))
        for(j in 1:length(min_model)){
          comb_fcast_componet[,j] <- forecast_list[min_model[j]][,i]
        }
        comb_fcast <- rowMeans(comb_fcast_components)
        real <- real.error[i, fm[[i]]]
        training <- tslist[[i]]$x
        test <- tslist[[i]]$xx
        MASEOpt[i] <- accmat(training, test, forecast=comb_fcast)
      }
    }

    summary_prediction_results <- rbind(c(mean(MASEOpt), apply(real.error, 2, mean)),
                                        c(median(MASEOpt), apply(real.error, 2, median)))
    dimnames(summary_prediction_results) <- list(c("mean", "median"), c("our_method", colnames(real.error)))

    return(list(models=fm, minmase=MASEOpt, summary=summary_prediction_results))
}
