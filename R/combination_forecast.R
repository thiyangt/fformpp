#' Function to generate combination forecast
#'
#' Select "n" number of algorithms based on minimum forecast error
#' @param ncomp number of algorithms use to combine
#' @param predicted prediction value matrix
#' @param real.error optional, matrix of MASE values for all algorithms for test data
#' @param accmat function to compute forecast accuracy
#' @param tslist list of time series, as in the format of Mcomp object
#' @param forecast_list true forecast from different models
#' @param h length of forecast horizon
#' @param weights weighted by MASE
#' @return list containing the forecasts and summaries
#' @export
combination_forecast <- function(predicted, ncomp=2, accmat=NULL, real.error=NULL, tslist=TRUE, forecast_list=NULL, h=NULL, weights=TRUE){
  tpredicted <- t(predicted)
  pred.list <- lapply(seq_len(ncol(tpredicted)), function(i) tpredicted[,i])
  fm <- lapply(pred.list, function(temp, ncomp) {
    which(temp %in% sort(unique(temp))[1:ncomp])
  }, ncomp=ncomp)

  fm.value <- lapply(pred.list, function(temp, ncomp) {
    weights <- temp[temp %in% sort(unique(temp))[1:ncomp]]
    weights/sum(weights)
  }, ncomp=ncomp)

  MASEOpt <- rep(NA, length(fm))
  for (i in 1: length(fm)){

      min_model <- colnames(predicted)[fm[[i]]]
      comb_fcast_components <- matrix(NA, ncol=length(fm[[i]]), nrow=h)
      for(j in 1:length(fm[[i]])){
        if(weights==TRUE){
        comb_fcast_components[,j] <- forecast_list[[min_model[j]]][,i]*fm.value[[i]][j]
        } else {
        comb_fcast_components[,j] <- forecast_list[[min_model[j]]][,i]
        }
      }


      comb_fcast <- apply(comb_fcast_components, 1, mean)
      real <- real.error[i, fm[[i]]]
      training <- tslist[[i]]$x
      test <- tslist[[i]]$xx
      MASEOpt[i] <- accmat(training, test, forecast=comb_fcast)
    }


  summary_prediction_results <- rbind(c(mean(MASEOpt), apply(real.error, 2, mean)),
                                      c(median(MASEOpt), apply(real.error, 2, median)))
  dimnames(summary_prediction_results) <- list(c("mean", "median"), c("our_method_comb", colnames(real.error)))

  return(list(models=fm, minmase=MASEOpt, summary=summary_prediction_results))

}
