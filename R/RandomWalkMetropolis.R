##' Random walk Metropolis algorithm.
##'
##' <details>
##' @title Random Walk Metropolis
##' @param param.cur
##' @param logpost.fun.name
##' @param Params
##' @param Y
##' @param x0
##' @param callParam
##' @param splineArgs
##' @param priorArgs
##' @param Params_Transform
##' @return
##' @references Li Villani 2012
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: Fri Aug 31 13:44:15 CEST 2012;
##'       Current: Fri Aug 31 13:44:21 CEST 2012.
##' @export
RandomWalkMetropolis <- function(
    param.cur,
    logpost.fun.name,
    Params,
    Y,
    x0,
    callParam,
    splineArgs,
    priorArgs,
    Params_Transform)
  {
    ## All the propose are made at the level of transformed stage.
    ## No need further transformation

    ## Random walk proposal
    nprop <- length(param.cur)
    error.prop.mean <- matrix(0, nrow = 1, ncol = nprop)
    error.prop.sigma <- diag(0.2, nrow = nprop, ncol = nprop) # diag() is slow
                                        # prop sd is hard coded. use a better
                                        # strategy

    error.prop <- rmvnorm(n = 1, mean = error.prop.mean,
                          sigma = error.prop.sigma) # 1-by-nprop

    param.prop <- param.cur + matrix(error.prop)

    ## The jump densities
    ## Symmetric with random walk case
    logjump.cur2prop <- 0 # the jump density from current draw to propose draw.
    logjump.prop2cur <- 0 # the jump density from propose draw to current draw.

    ## Calculate the ratio of the densities
    Params.prop <- Params
    Params.prop[[callParam$id]][callParam$subset] <- param.prop
    Params.cur <- Params
    Params.cur[[callParam$id]][callParam$subset] <- param.cur

    caller.prop <- call(logpost.fun.name,Y = Y, x0 = x0, Params = Params.prop,callParam =
                        callParam ,  priorArgs = priorArgs, splineArgs = splineArgs,
                        Params_Transform = Params_Transform)
    caller.cur <- call(logpost.fun.name, Y = Y,  x0 = x0, Params = Params.cur,callParam =
                       callParam, priorArgs = priorArgs, splineArgs = splineArgs,
                       Params_Transform = Params_Transform)

    logpost.prop <-  eval(caller.prop) # the logpost density for the proposed draw.
    logpost.cur <- eval(caller.cur) # the logpost density for the current draw.

    log.r <- logpost.prop - logpost.cur + logjump.prop2cur - logjump.cur2prop
    r <- exp(log.r)     ## The MH ratio.

    accept.prob <- min(1, r) # the acceptance probability.

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

    ##cat("prop:", param.prop, "cur:", param.cur, "\n")
    return(out)

  }
