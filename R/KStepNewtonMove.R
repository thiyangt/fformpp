##' This is the simple Newton move for spline models without dimession changes.
##'
##'
##' @name
##' @title
##' @param param.cur "matrix".
##'         The initial values for the Newton update.
##' @param gradhess.fun.name
##' @param KStep
##' @param Params
##' @param hessMethod
##' @param Y
##' @param x
##' @param callParam
##' @param splineArgs
##' @param priorArgs
##'        priorArgs$prior_type:
##'        priorArgs$ka0:
##'        priorArgs$M:
##'        priorArgs$S0:
##'        priorArgs$n0:
##'        priorArgs$mu0:
##'        priorArgs$Sigma0:
##'        priorArgs$ka.mu0:
##'        priorArgs$ka.Sigma0:
##' @return "list". See bellow.
##' \item   {gradObs.cur}
##'         {"matrix". The gradient}
##' \item   {invHessObs.cur}
##'         {"matrix". The inverse hessian matrix.}
##' \item   {param.cur}
##'         {"matrix". The updated paramters after K step Newton integrations.}
##' @references
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note First version: Wed Sep 29 17:18:22 CEST 2010;
##'       Current:       Wed Sep 29 17:18:34 CEST 2010.
##' TODO: Error handling when hessian is ill behaved.
##'       Think about how to add shrinkage in.
##'       And variable selection
##' @export
KStepNewtonMove <- function(param.cur, gradhess.fun.name, KStep, Params,
                            hessMethod, Y, x0, callParam, splineArgs,
                            priorArgs, Params_Transform)
{

    ## Don't have to consider the linkage. It was taken care of out before calling this function

    ## K step to approach the mode plus one more step to update the gradient and hessian at
    ## K:th step.
    for(k in 1:(KStep+1))
    {
        ## Update the parameter with current updated results. The first update of course
        ## is the initial value.
        Params[[callParam$id]][callParam$subset] <- param.cur
        caller <- call(gradhess.fun.name, Params = Params, hessMethod = hessMethod, Y = Y, x0
                       = x0, callParam = callParam, splineArgs = splineArgs,
                       priorArgs = priorArgs, Params_Transform = Params_Transform)  # Creat the gradien object. Will be passed to eval().
        gradhess <- eval(caller) # Calculate the hessian at current draw. If the
                                 # parameters are good enough, do one more update the
                                 # gradient and return.

        gradObs.cur <- gradhess$gradObs
        hessObs.cur <- gradhess$hessObs

        ## Check if gradient or Hessian is bad
        if(any(is.na(gradObs.cur)) || any(is.na(hessObs.cur))) # bad, return NaN and quit, this draw will be
                                        # discarded in MH step
        {
            out <- list(gradObs.cur = NaN, hessObs.cur = NaN, invHessObs.cur = NaN,
                        param.cur = NaN)
        }
        else # good, do k-step newton
        {
            invHessObs.cur <- ginv(as.matrix(hessObs.cur)) # Convert Matrix class to matrix

            if((k <= KStep)) # if need to update newton steps
            {
                ## cat(k, "step: param.cur: ", param.cur, "\n")
                ## cat(k, "step: gradObs.cur", gradObs.cur, "\n")
                ## cat(k, "step: -invHessObs.cur", -invHessObs.cur, "\n")
                ## if (is.na(gradObs.cur)) browser()
                param.cur <- param.cur-invHessObs.cur%*%gradObs.cur # update the parameters
            }
            else if(k == (KStep+1)) # The parameters at (K+1):th step. Fine enough to make a output
            {
                param.out <- param.cur # let last update be the final out put
                out <- list(gradObs.cur = gradObs.cur,
                            hessObs.cur = hessObs.cur,
                            invHessObs.cur = invHessObs.cur,
                            param.cur = param.out)
            }
        }
    }
    return(out)
}

## don't put return within for loop,  it will make an error,  see this example
## for(i in 1:10)
## {
##   cat("at stage", i, "\n")

##   if(i <= 5)
##     {
##       cat("i = ", i, "continuous...\n")
##     }
##   else if(i  == 6)
##     {
##       return("End\n")
##     }
## }
