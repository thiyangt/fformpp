##' Stochastic MCMC using Stochastic gradient Langevin dynamics
##'
##' <description>
##' @param param.cur NA
##' @param logpost.fun.name MA
##' @param Params NA
##' @param Y NA
##' @param x0 NA
##' @param callParam NA
##' @param splineArgs NA
##' @param priorArgs NA
##' @param Params_Transfor NA
##' @return NA
##' @references Ma, Chen, & Fox 2015 A Complete Recipe for Stochastic Gradient MCMC.
##' @author Feng Li, Central University of Finance and Economics.
##' @export
SGLD = function(param.cur,
                logpost.fun.name,
                Params,
                Y,
                x0,
                callParam,
                splineArgs,
                priorArgs,
                algArgs,
                Params_Transform)
{
    minibatchProp = algArgs[["minibatchProp"]]
    nEpoch= algArgs[["nEpoch"]]
    calMHAccRate = algArgs[["calMHAccRate"]] # if TRUE, SGLD will be the proposal of
                                             # MH. See Welling & Teh (2011), p 3.
    stepsizeSeq = algArgs[["stepsizeSeq"]]

    nPar = length(param.cur)
    nObs = dim(Y)[1]

    nIterations = ceiling(1/ minibatchProp) # no. runs with one epoch.
    nRuns = (nIterations * nEpoch)
    param.out = matrix(NA, nPar, nRuns)
    accept.probMat = matrix(NA, 1, nRuns)

    for(iRun in 1:nRuns)
    {
        ## Redo a splitting when finishing one epoch
        if(iRun %in% seq(1, nRuns, by = nIterations))
        {
            suppressWarnings(subsetIdxLst <- split(sample(1:nObs, size = nObs), 1:nIterations))
        }
        whichIteration = iRun %% nIterations
        whichIteration[whichIteration == 0] = nIterations

        subsetIdx = subsetIdxLst[[whichIteration]]
        epsilon = stepsizeSeq[iRun]
        minibatchSizeNew = length(subsetIdx) # The minibatch size is not exactly same as
                                        # required due to unequal splitting.
        Params[[callParam$id]][callParam$subset] <- param.cur
        caller = call(gradhess.fun.name,
                      Params = Params,
                      hessMethod = "skip",
                      Y = Y[subsetIdx, , drop = FALSE],
                      x0 = x0[subsetIdx, , drop = FALSE],
                      callParam = callParam,
                      splineArgs = splineArgs,
                      priorArgs = priorArgs,
                      Params_Transform = Params_Transform)
        gradhess <- eval(caller)

        gradObs.logLik <- gradhess$gradObs.logLik
        gradObs.logPri <- gradhess$gradObs.logPri

        gradObsSub = (gradObs.logLik / minibatchSizeNew * nObs + gradObs.logPri)

        ## Stochastic Gradient Riemannian Langevin Dynamics (SGRLD). Notation followed in Ma 2015
        ## fisherInfoMatG = diag(as.vector(1 / param.cur), nrow = nPar)
        ## diffusionMatD = diag(as.vector(param.cur), nrow = nPar)

        fisherInfoMatG = diag(1, nrow = nPar)
        diffusionMatD = diag(1, nrow = nPar)

        correctionGamma = matrix(0, nPar, 1)
        potentialUgrad = -gradObsSub

        noise = matrix(rmvnorm(n = 1, mean = matrix(0, nrow = 1, ncol = nPar),
                               sigma = 2 * epsilon * diffusionMatD)) # row vector -> col vector

        param.prop = (param.cur - epsilon * (diffusionMatD %*% potentialUgrad +
                                             correctionGamma) + noise)


        if(calMHAccRate == TRUE)
        {
            stop("Not yet done!")
        }
        else
        {
            ## No Metropolis-Hastings calibration
            param.cur = param.prop
            param.out[, iRun] = param.prop
            accept.probMat[, iRun] = 1
        }
    }

    out <- list(param.out = param.out, accept.prob = mean(accept.probMat))
    return(out)

}
