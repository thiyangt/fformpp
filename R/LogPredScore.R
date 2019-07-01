##' This is the LPDS
##'
##' Details.
##' @name
##' @title
##' @param Y
##' @param x
##' @param logpost.fun.name
##' @param crossvaid.struc
##' @param splineArgs
##' @param priorArgs
##' @param OUT.Params
##' @param Params_Transform
##' @param burn.in
##' @return
##' @references
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note First version: Tue Nov 30 23:18:18 CET 2010;
##'       Current:       Tue Nov 30 23:18:25 CET 2010.
##' TODO: Don't use the full draws when nIter is large.
##'       Still doubt the nseLPDS
##' @export
LogPredScore <- function(Y, x, logpost.fun.name, crossvaid.struc, splineArgs, priorArgs,
                 OUT.Params, Params_Transform, burn.in, LPDS.sampleProp)
{

  n <- dim(Y)[1] # no. of obs.
  nIter <- dim(OUT.Params[[1]])[4] # no. of iterations
  num.burn.in <- floor(nIter*burn.in)
  nUsed <- nIter - num.burn.in

  ## The sample indices for LPDS after burn-in
  ## Make a decreasing sequence and sort it by increasing order.
  ## Just to make sure last draw is allays used
  LPDS.sampleIdx <- sort(seq(nIter, (nIter-nUsed+1), by = -round(1/LPDS.sampleProp))) # The exact
                                        # size may not same as the result from sample
                                        # proportion.
  nSample <- length(LPDS.sampleIdx)

  nCross <- length(crossvalid.struc[["training"]])

  if(nCross  == 1)
    {
      LPDS  = NA
      nseLPDS = NA
      logPredMatrix <- NA
    }
  else if(nCross  != 1)
    {

      ## Detect how many folds used
      if(length(crossvalid.struc[["training"]][[nCross]]) +
         length(crossvalid.struc[["testing"]][[nCross]])  == n)
        { nFold <- nCross}
      else
        { nFold <- nCross -1}

      ## The log predictive matrix
      logPredMatrix <- matrix(NA, nSample, nFold)

      ## Calculate the predictive densities for all folds
      for(iCross in 1:nFold)
        {
         # iTraining <- crossvalid.struc[["training"]][[iCross]]
          iTesting <- crossvalid.struc[["testing"]][[iCross]]

          ## Y.iTraining <- Y[iTraining, , drop = FALSE]
          ## x.iTraining <- x[iTraining, , drop = FALSE]

          Y.iTesting <- Y[iTesting, , drop = FALSE]
          x.iTesting <- x[iTesting, , drop = FALSE]

          which.j <- 0
          for(j in LPDS.sampleIdx) ## Just the likelihood function with posterior samples
            {
              Params.j <- lapply(OUT.Params, function(x) apply(x[,,, j, iCross, drop =
                                                                 FALSE], c(1, 2), "["))
              caller.log.like <- call(logpost.fun.name,Y = Y.iTesting, x = x.iTesting,
                                      Params = Params.j, callParam = list(id =
                                                           "likelihood"), priorArgs =
                                      priorArgs, splineArgs = splineArgs, Params_Transform
                                      = Params_Transform)
              log.like <- eval(caller.log.like)

              which.j <- which.j + 1
              logPredMatrix[which.j, iCross] <- log.like

              ## Simple progress bar
              ## progressbar(((iCross-1)*nSample + which.j), nFold*nSample)
            }

        }

      ## Calculate the LPDS
      scaleFactors <- apply(logPredMatrix, 2, median)
      scaleMatrix <- matrix(scaleFactors, nSample, nFold)
      expPredMatrix <- exp(logPredMatrix-scaleMatrix)
      expPredMatrix[is.infinite(expPredMatrix)] <- NA # TODO: Think about it carefully
      expPredMean <- colMeans(expPredMatrix, na.rm = TRUE)
      LPDS <- mean(scaleFactors + log(expPredMean))

      if (sum(is.infinite(expPredMatrix)) > nSample*.05) # If at least 5% are infinity,
                                        # maybe due to unsuccessful MCMC.
        {
          LPDS <- NA
          warning("Too many infinite desities produced, check MCMC convergence!")
        }
      ## Calculate the numerical se. for LPDS. See Box Jenkins (2007)  p.31
      ACFxx <- matrix(0, 1, nFold)
      for(k in 1:nFold)
        {
          acfSample.tmp0 <- expPredMatrix[, k]
          acfSample.tmp <- acfSample.tmp0[acfSample.tmp0<1e100] # Just let it numerically stable.

          predACF <- acf(acfSample.tmp, plot = FALSE,  na.action = na.pass)$acf
          nlagk <- length(predACF)
          for(kk in 0:(nlagk-1))
            {
              ACFxx[k] <- ACFxx[k] + (1 - kk/nlagk)*predACF[kk +1]
            }
        }

      expPredMatrix.tmp <- expPredMatrix
      expPredMatrix.tmp[expPredMatrix.tmp>1e100] <- NA # numerically stable
      predVar <- apply(expPredMatrix.tmp, 2, var, na.rm = TRUE)

      ## var.MeanexpPredMatrix <- predVar/nUsed*(1+2*ACFxx)
      var.MeanexpPredMatrix <- predVar/nSample*(1+2*ACFxx)

      nvarLPDS <- 1/nFold^2*sum(1/(expPredMean)^2*var.MeanexpPredMatrix)
      nseLPDS <- sqrt(nvarLPDS)
    }

  out <- list(LPDS = LPDS, nseLPDS = nseLPDS, logPredMatrix = logPredMatrix, LPDS.sampleIdx = LPDS.sampleIdx)
  return(out)
}
