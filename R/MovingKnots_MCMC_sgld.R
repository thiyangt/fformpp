#' @export
MovingKnots_MCMC_sgld <- function(gradhess.fun.name,
                             logpost.fun.name,
                             nIter,
                             Params,
                             Params4Gibbs,
                             Params.sub.struc,
                             Y,
                             x0,
                             splineArgs,
                             priorArgs,
                             Params_Transform,
                             propMethods,
                             algArgs,
                             crossvalid.struc,
                             OUT.Params,
                             OUT.accept.probs,
                             burn.in,
                             LPDS.sampleProp,
                             track.MCMC){
    ##----------------------------------------------------------------------------------------
    ## Update Knots locations (subsets),  shrinkage and covariance jointly
    ##----------------------------------------------------------------------------------------
    Running.date <- Sys.time() # The stating time
    Start.Time <- proc.time() # The CPU time

    cat("Updating Knots, Shrinkages, and Covariance >>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n\n")
    MCMCPropFun <- function(iCross)
    {
        ## Training sample
        Y.iCross <- Y[crossvalid.struc$training[[iCross]], , drop = FALSE]
        x.iCross <- x0[crossvalid.struc$training[[iCross]], , drop = FALSE]
        for(iIter in 1:nIter) # loop nIter times
        {
            ## a <- proc.time()
            for (iPar in Params4Gibbs) # loop over parameters
            {
                for(iSub in 1:length(Params.sub.struc[[iPar]])) # update subsets
                {
                    Sub4iPar <- Params.sub.struc[[iPar]][[iSub]]
                    if(!is.null(Sub4iPar)) # update if subsets are all fixed
                    {
                        param.cur <- matrix(Params[[iPar]][Sub4iPar])
                        algArgs.cur = algArgs[[iPar]]

                        ## Special case to pass stepsize sequence.
                        if(tolower(propMethods[[iPar]]) %in% c("sgld"))
                        {
                            nInner = length(algArgs.cur[["stepsizeSeq"]]) / nIter
                            algArgs.cur[["stepsizeSeq"]] = (algArgs.cur[["stepsizeSeq"]]
                                [((iIter - 1)* nInner + 1):(iIter * nInner)])
                        }

                        out.iSub <- MHPropMain_sgld(param.cur = param.cur,
                                               gradhess.fun.name = gradhess.fun.name,
                                               logpost.fun.name = logpost.fun.name,
                                               Params = Params,
                                               Y.iCross = Y.iCross,
                                               x.iCross = x.iCross,
                                               callParam = list(id = iPar, subset = Sub4iPar),
                                               splineArgs = splineArgs,
                                               priorArgs = priorArgs,
                                               algArgs = algArgs[[iPar]],
                                               Params_Transform = Params_Transform,
                                               propMethod = propMethods[[iPar]])

                        ## Update the parameters in the parameters list.
                        param.cur.outMat = out.iSub$param.out

                        ## Take the last one if inner loops are used in e.g. SGLD
                        Params[[iPar]][Sub4iPar] <- param.cur.outMat[, ncol(param.cur.outMat)]
                        ## Save the acceptance probability
                        OUT.accept.probs[[iPar]][iSub, iIter, iCross] <- out.iSub$accept.prob

                        ## Save the updated parameters for current iteration.
                        OUT.Params[[iPar]][Sub4iPar, , , iIter, iCross] <- param.cur.outMat
                    }

                }
            } # for (iPar in Params4Gibbs)
            ## Track the iterations
            if(track.MCMC)
            {
                if(nIter>1000)
                {
                    interval = 0.20
                }
                else
                {
                    interval = 1/nIter
                }
                MCMC.trajectory(iIter, nIter, iCross, OUT.accept.probs, interval = .10)
            }
            ## print(proc.time()-a)
        }
        return(list(OUT.Params = OUT.Params, OUT.accept.probs = OUT.accept.probs))
    } # for(iCross in 1:nCross)

    nCross <- length(crossvalid.struc$training)

    cl <- getDefaultCluster()
    if(length(cl) != 0)
    {
        clusterExport(cl, "nCross")
        sink.parallel(cl)
        MCMCPropOut = parLapply(cl = cl, X = as.list(1:nCross), fun = MCMCPropFun)
        sink.parallel(cl, file = NULL)
    }
    else
    {
        MCMCPropOut = lapply(X = as.list(1:nCross), FUN = MCMCPropFun)
    }

    ## Collecting results
    for(iCross in 1:nCross){
        for (iPar in Params4Gibbs) # loop over parameters
        {
            OUT.Params[[iPar]][,,,, iCross] = MCMCPropOut[[iCross]][["OUT.Params"]][[iPar]][,,,, iCross]
            OUT.accept.probs[[iPar]][,, iCross] = MCMCPropOut[[iCross]][["OUT.accept.probs"]][[iPar]][,, iCross]
        }
    }

    ##----------------------------------------------------------------------------------------
    ## Sample coefficients from Normal
    ##----------------------------------------------------------------------------------------
    cat("Updating Coefficients >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n\n")

    ## Average inner draws but keep the output structure.
    OUT.Params1Inner = lapply(OUT.Params, function(x){
        x4 = apply(x, c(1, 2, 4, 5), mean)
        dimx = dim(x)
        dimx[3] = 1
        array(x4, dimx)
    })

    OUT.Params <- linear_post4coef(Y = Y,
                                   x0 = x0,
                                   OUT.Params = OUT.Params1Inner,
                                   crossvalid.struc = crossvalid.struc,
                                   nCross = nCross,
                                   nIter = nIter,
                                   splineArgs = splineArgs,
                                   priorArgs = priorArgs,
                                   Params_Transform = Params_Transform)

    ##----------------------------------------------------------------------------------------
    ## Compute the predictive density and LPDS
    ##----------------------------------------------------------------------------------------
    ## Update the LPDS
    cat("Updating LPDS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n\n")
    OUT.LPDS <- LogPredScore(Y = Y,
                             x = x0,
                             logpost.fun.name = logpost.fun.name,
                             crossvaid.struc = crossvaid.struc,
                             splineArgs = splineArgs,
                             priorArgs = priorArgs,
                             OUT.Params = OUT.Params,
                             Params_Transform = Params_Transform,
                             burn.in = burn.in,
                             LPDS.sampleProp = LPDS.sampleProp)

    cat("LPDS:", round(OUT.LPDS$LPDS, 2), "( n.se:", round(OUT.LPDS$nseLPDS, 2), ")",
        "\n\n")

    Finish.Time <- proc.time()

##########################################################################################
###                                 Post Analysis
##########################################################################################

    ##----------------------------------------------------------------------------------------
    ## Computing posterior modes
    ##----------------------------------------------------------------------------------------

    ## Drop the burn-in
    num.burn.in <- floor(nIter*burn.in)

    ## TODO: the posterior mode should be averaged with stepsize? Willing & Teh 2011, Eq(11).
    OUT.Params.mode <- lapply(OUT.Params, function(x)
        apply(x[,,, (num.burn.in+1):nIter,, drop=FALSE], c(1, 2, 5), mean))

    OUT.Params.sd <- lapply(OUT.Params, function(x)
        apply(x[,,, (num.burn.in+1):nIter,, drop=FALSE], c(1,  2, 5), sd))

    OUT.Params.ineff <- lapply(OUT.Params, function(x)
        apply(x[,,,(num.burn.in+1):nIter,, drop=FALSE], c(1, 2, 5), ineff))

    OUT.accept.probs.mean <- lapply(OUT.accept.probs, function(x)
        apply(x[, (num.burn.in+1):nIter, ,drop = FALSE],c(1, 3), mean))

    ##----------------------------------------------------------------------------------------
    ## Collecting system information
    ##----------------------------------------------------------------------------------------
    SYS.INFO <- list(Running.date = Running.date,
                     Elapsed.Time <- Finish.Time - Start.Time,
                     OMP_NUM_THREADS = as.numeric(Sys.getenv("OMP_NUM_THREADS")),
                     System.info = as.list(Sys.info()))

##########################################################################################
    ##                                Save output
##########################################################################################
    ## save important output to global environment. "<<-"
    OUT <- list()
    OUT[["Params"]] <- OUT.Params
    OUT[["Params.mode"]] <- OUT.Params.mode
    OUT[["Params.sd"]] <- OUT.Params.sd
    OUT[["Params.ineff"]] <- OUT.Params.ineff

    OUT[["accept.probs"]] <- OUT.accept.probs
    OUT[["accept.probs.mean"]] <- OUT.accept.probs.mean

    OUT[["LPDS"]] <- OUT.LPDS
    OUT[["SYS.INFO"]] <- SYS.INFO

    return(OUT)
}
