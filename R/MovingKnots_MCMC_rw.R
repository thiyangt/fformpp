#' @export
MovingKnots_MCMC_rw <- function(gradhess.fun.name, logpost.fun.name, nNewtonSteps, nIter,
                             Params, Params4Gibbs, Params.sub.struc, hessMethods, Y, x0, splineArgs,
                             priorArgs, MH.prop.df, Params_Transform, propMethods,
                             crossvalid.struc, OUT.Params, OUT.accept.probs, burn.in,
                             LPDS.sampleProp, track.MCMC)
{
  ##----------------------------------------------------------------------------------------
  ## Update Knots locations (subsets),  shrinkage and covariance jointly
  ##----------------------------------------------------------------------------------------
  Running.date <- Sys.time() # The stating time
  Start.Time <- proc.time() # The CPU time

  cat("Updating Knots, Shrinkages, and Covariance >>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n\n")

  for(iCross in 1:nCross) # loop over subsets of data. TODO: Parallel Computing?
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
            out.iSub <- MHPropMain_rw(param.cur = param.cur,
                                   gradhess.fun.name = gradhess.fun.name,
                                   logpost.fun.name = logpost.fun.name,
                                   nNewtonStep = nNewtonSteps[[iPar]],
                                   Params = Params,
                                   hessMethod = hessMethods[[iPar]],
                                   Y.iCross = Y.iCross,
                                   x.iCross = x.iCross,
                                   callParam = list(id = iPar, subset = Sub4iPar),
                                   splineArgs = splineArgs,
                                   priorArgs = priorArgs,
                                   prop.df = MH.prop.df[[iPar]],
                                   Params_Transform = Params_Transform,
                                   propMethod = propMethods[[iPar]])
            ## Update the parameters in the parameters list.
            Params[[iPar]][Sub4iPar] <- out.iSub$param.out
            ## Save the acceptance probability
            OUT.accept.probs[[iPar]][iSub, iIter, iCross] <- out.iSub$accept.prob
          }
        }
        ## Save the updated paramters for current iteration.
        OUT.Params[[iPar]][,  , iIter, iCross] <- Params[[iPar]]
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
  } # for(iCross in 1:nCross)

  ##----------------------------------------------------------------------------------------
  ## Sample coefficients from Normal
  ##----------------------------------------------------------------------------------------
  cat("Updating Coefficients >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n\n")


  OUT.Params <- linear_post4coef(Y, x0, OUT.Params, crossvalid.struc, nCross, nIter,
                                 splineArgs, priorArgs, Params_Transform)

  ##----------------------------------------------------------------------------------------
  ## Compute the predictive density and LPDS
  ##----------------------------------------------------------------------------------------

  ## Update the LPDS
  cat("Updating LPDS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n\n")
  OUT.LPDS <- LogPredScore(Y = Y, x = x0, logpost.fun.name = logpost.fun.name,
                           crossvaid.struc = crossvaid.struc, splineArgs = splineArgs,
                           priorArgs = priorArgs, OUT.Params = OUT.Params, Params_Transform
                           = Params_Transform, burn.in = burn.in, LPDS.sampleProp = LPDS.sampleProp)
  cat("LPDS:", round(OUT.LPDS$LPDS, 2), "( n.se:", round(OUT.LPDS$nseLPDS, 2), ")",
      "\n\n")

  Finish.Time <- proc.time()

  ##########################################################################################
  ##                                 Post Analysis
  ##########################################################################################

  ##----------------------------------------------------------------------------------------
  ## Computing posterior modes
  ##----------------------------------------------------------------------------------------

  ## Drop the burn-in
  num.burn.in <- floor(nIter*burn.in)

  OUT.Params.mode <- lapply(OUT.Params,
                            function(x) apply(x[, , (num.burn.in+1):nIter,drop=FALSE, ],
                                              c(1, 2, 4), mean))
  OUT.Params.sd <- lapply(OUT.Params,
                          function(x) apply(x[, , (num.burn.in+1):nIter,drop=FALSE, ],
                                            c(1, 2, 4), sd))
  OUT.Params.ineff <- lapply(OUT.Params,
                             function(x) apply(x[, , nIter:(num.burn.in+1),drop=FALSE, ],
                                               c(1, 2, 4), ineff))

  OUT.accept.probs.mean <- lapply(OUT.accept.probs,
                                  function(x) apply(x[, (num.burn.in+1):nIter, ,drop = FALSE],c(1, 3), mean))

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
