##' Random walk Metropolis–Hastings algorithm for Sigma
##'
##' Details.
##' @name
##' @title
##' @param param.cur
##' @param logpost.fun.name
##' @param Params
##' @param Y
##' @param x
##' @param callParam
##' @param splineArgs
##' @param priorArgs
##' @param Params_Transform
##' @return
##' @references
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note First version: Mon Nov 15 18:23:37 CET 2010;
##'       Current:       Wed Feb 23 09:26:07 CET 2011.
##' DEPENDS:
##' @export
MHPropWithIWishart <- function(param.cur, logpost.fun.name, Params, Y, x0, callParam,
                               splineArgs, priorArgs, Params_Transform)
{
  # browser()
  ## Update param.cur to the Params list
  Params[[callParam$id]][callParam$subset] <- param.cur

  ## The current Parameter and transform back to original scale
  param.cur.TB <- par.transform(param.cur, method = Params_Transform[[callParam$id]])
  Sigma.cur <- vech2m(param.cur.TB) # symmetric form

  ## Calculate the posterior V and df given current parameters.
  par.cur.list <- linear_IWishart(Params, Y, x0, splineArgs, priorArgs, Params_Transform)

  ## browser()
  V.cur.prop <- par.cur.list[["V"]]
  prop.df <- par.cur.list[["df"]]

  ## Check if they are all right
  if(is.na(V.cur.prop) || is.na(prop.df)) # bad proposal
    {
      out <- list(param.out = param.cur, accept.prob = 0) # keep current value and quit
      return(out)
    }

  ## Propose a draw
  Sigma.prop <- riwishart(df = prop.df, V = V.cur.prop) # propose one draw, full matrix

  ## print(Sigma.prop)
  ## Calculate the jumping densities
  ## Random walk Metropolis Hasting
  logjump.cur2prop <- diwishart(Sigma.prop, df = prop.df,  V = V.cur.prop, log = TRUE) # the jump density from current draw to propose draw.

  logjump.prop2cur <- diwishart(Sigma.cur, df = prop.df, V = V.cur.prop, log = TRUE) # the
                                        # jump density from propose draw to current draw.

  ## Calculate the posterior densities for both current and proposal draw.
  ## Transform the scale to the final scale used in Params
  Params.prop <- Params
  param.prop <- par.transform2(vech(Sigma.prop), Params_Transform[[callParam$id]])
  Params.prop[[callParam$id]][callParam$subset] <- param.prop

  Params.cur <- Params # The current Params list

  caller.prop <- call(logpost.fun.name,Y = Y, x0 = x0, Params = Params.prop,callParam =
                      callParam ,  priorArgs = priorArgs, splineArgs = splineArgs,
                      Params_Transform = Params_Transform)
  caller.cur <- call(logpost.fun.name, Y = Y,  x0 = x0, Params = Params.cur,callParam =
                     callParam, priorArgs = priorArgs, splineArgs = splineArgs,
                     Params_Transform = Params_Transform)

  logpost.prop <-  eval(caller.prop) # the logpost density for the proposed draw.
  logpost.cur <- eval(caller.cur) # the logpost density for the current draw.

  ## compute the MH ratio
  log.r <- logpost.prop - logpost.cur + logjump.prop2cur - logjump.cur2prop
  r <- exp(log.r)

  ## Compute the acceptance probability
  accept.prob <- min(1, r) # the acceptance probability.

  ## make desicion to update or keep current draw.
  if(!is.na(accept.prob) && runif(1) < accept.prob)
    {
      param.out <- param.prop  # keep update, output with final scale
    }
  else
    {
      param.out <- param.cur # keep current, output with final scale
      accept.prob <- 0 # set acceptance probs to zero.
    }

  out <- list(param.out = param.out, accept.prob = accept.prob)
  return(out)
}


## V <- diag(c(0.2, 0.3, 0.4), 3)
## df <- 20

## n <- 1000
## mat <- array(0, c(3, 3, n))
## for(i in 1:n)
##   {
##     mat[, , i] <- riwishart(df, V)
##   }
