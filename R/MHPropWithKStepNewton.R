##' Metropolis–Hastings algorithm with K-step Newton method for the spline model.
##'
##' Details.
##' @name
##' @title
##' @param param.cur
##' @param gradhess.fun.name
##' @param logpost.fun.name
##' @param nNewtonStep
##' @param Params
##' @param hessMethod
##' @param Y
##' @param x
##' @param callParam
##' @param splineArgs
##' @param priorArgs
##' @param prop.df
##' @param Params_Transform
##' @return
##' @references
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note First version: Thu Feb 17 14:03:14 CET 2011;
##'       Current:       Wed Feb 23 14:43:35 CET 2011.
##' DEPENDS: mvtnorm
##' TODO: replace the old multivariate t functions to mvtnorm functions
##' @export
MHPropWithKStepNewton <- function(param.cur, gradhess.fun.name, logpost.fun.name,
                                  nNewtonStep, Params, hessMethod, Y, x0, callParam,
                                  splineArgs, priorArgs, prop.df, Params_Transform)
{
    ## All the propose are made at the level of transformed stage. No need further
    ## transformation
    ## print(callParam)
    ## browser()
    KStepNewton1 <- KStepNewtonMove(param.cur = param.cur,
                                    gradhess.fun.name = gradhess.fun.name,
                                    KStep = nNewtonStep,
                                    Params = Params,
                                    hessMethod = hessMethod,
                                    Y = Y,
                                    x0 = x0,
                                    callParam = callParam,
                                    splineArgs = splineArgs,
                                    priorArgs = priorArgs,
                                    Params_Transform = Params_Transform)

    param.cur.prop <- KStepNewton1$param.cur
    HessObs.cur.prop <- KStepNewton1$hessObs.cur
    invHessObs.cur.prop <- KStepNewton1$invHessObs.cur

    if(any(is.na(HessObs.cur.prop)) ||
       is(try(chol(-invHessObs.cur.prop), silent=T), "try-error")) # Something is wrong,  reject it.
    {
        logjump.cur2prop <- NaN
        param.prop <- NaN
    }
    else
    {
        # param.prop <- RndMultiT(param.cur.prop, -invHessObs.cur.prop, prop.df) # propose a draw
                                        # from multivariate t-distribution.
        param.prop <- (param.cur.prop +
                       t(rmvt(sigma = -invHessObs.cur.prop,
                              n = 1, df = prop.df,
                              method = "chol")))

        # logjump.cur2prop <- DensMultiT(param.prop, param.cur.prop, -HessObs.cur.prop, prop.df)
                                        # the jump density from current draw to propose
                                        # draw.
        logjump.cur2prop = dmvt(x = matrix(param.cur.prop - param.prop, 1, ),
                                sigma = -invHessObs.cur.prop,
                                df = prop.df, log = TRUE)


    }

    if (any(is.na(param.prop))) # Previous proposal unsuccessful
    {      HessObs.prop.prop <- NaN    }
    else # all are right
    {
        KStepNewton2 <- KStepNewtonMove(param.cur = param.prop,
                                        gradhess.fun.name = gradhess.fun.name,
                                        KStep =nNewtonStep,
                                        Params = Params,
                                        hessMethod = hessMethod,
                                        Y = Y, x0 = x0,
                                        callParam = callParam,
                                        splineArgs = splineArgs,
                                        priorArgs = priorArgs,
                                        Params_Transform = Params_Transform)

        param.prop.prop <- KStepNewton2$param.cur
        HessObs.prop.prop <- KStepNewton2$hessObs.cur
        invHessObs.prop.prop <- KStepNewton2$invHessObs.cur
    }

    if(any(is.na(HessObs.prop.prop))) # Something is wrong at KStepNewton2,  reject it.
    {
        logpost.cur <- NaN
        logpost.prop <- NaN
        logjump.prop2cur <- NaN
    }
    else # all are right
    {

        ## logjump.prop2cur <- DensMultiT(param.cur, param.prop.prop, -invHessObs.prop.prop,
        ##                                prop.df) # the jump density from propose draw to current
                                        # draw.
        logjump.prop2cur = dmvt(x = matrix(param.prop.prop - param.cur, 1, ),
                                sigma = -invHessObs.prop.prop,
                                df = prop.df, log = TRUE)


        Params.prop <- Params
        Params.prop[[callParam$id]][callParam$subset] <- param.prop
        Params.cur <- Params
        Params.cur[[callParam$id]][callParam$subset] <- param.cur

        caller.prop <- call(logpost.fun.name,
                            Y = Y,
                            x0 = x0,
                            Params = Params.prop,
                            callParam = callParam ,
                            priorArgs = priorArgs,
                            splineArgs = splineArgs,
                            Params_Transform = Params_Transform)
        caller.cur <- call(logpost.fun.name,
                           Y = Y,
                           x0 = x0,
                           Params = Params.cur,
                           callParam = callParam,
                           priorArgs = priorArgs,
                           splineArgs = splineArgs,
                           Params_Transform = Params_Transform)

        logpost.prop <-  eval(caller.prop) # the logpost density for the proposed draw.
        logpost.cur <- eval(caller.cur) # the logpost density for the current draw.
    }

    ## compute the MH ratio.
    log.r <- logpost.prop - logpost.cur + logjump.prop2cur - logjump.cur2prop
    r <- exp(log.r)

    accept.prob <- min(1, as.numeric(r)) # the acceptance probability.

    ## make decision to update or keep current draw.
    if(!is.na(accept.prob) && runif(1) < accept.prob)
    {
        param.out <- param.prop  # keep update
    }
    else
    {
        param.out <- param.cur # keep current
        accept.prob <- 0 # set acceptance probs to zero.
    }
    out <- list(param.out = param.out, accept.prob = accept.prob)

    ## cat("prop:", param.prop, "cur:", param.cur, "\n")
    return(out)
}
