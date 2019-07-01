##' The conditional and joint log posterior function
##'
##' Details.
##' @param Y
##' @param x
##' @param Params
##' @param callParam
##' @param splineArgs "list".
##' @param ParamsTransArgs "list"
##' @param Params
##'         Params$xi:
##'         Parmas$K:
##' @param priorArgs
##'         priorArgs$prior_type:
##'         priorArgs$n0:
##'         priorArgs$S0:
##'         priorArgs$mu:
##'         priorArgs$M:
##'         priorArgs$mu0:
##'         priorArgs$Sigma0:
##'         priorArgs$ka0:
##'         priorArgs$ka.mu0:
##'         priorArgs$ka.Sigma0:
##'         priorArgs$Sigma.mu0:
##'         priorArgs$Sigma.Sigma0:
##' @param callParam
##'        callParam$id: It should be able to obtain conditianl posterior or joint
##'        posterior or just likelihood.
##'        callParam$subset: If the paramater set is larg, we may only update a subset of
##'        them.
##'
##' @return "scalar".
##' @references
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note First version: Sat Nov 06 22:23:57 CET 2010;
##'       Current:       Wed Dec 29 15:55:51 CET 2010.
##' @export
glm_logpost <- function(Y, x, Params, callParam, splineArgs, priorArgs, Params_Transform)
{

  ## Transform back when Params has been transformed.
  ParamsTB <- mapply(par.transform, par = Params, method = Params_Transform, SIMPLIFY =
                     FALSE)

  ## Get the knots name
  comp <- splineArgs[["comp"]]
  knots.name <- comp[! comp %in% c("intercept", "covariates")]

  ## Get the parameters
  diag.K <- ParamsTB[["shrinkages"]]
  Sigma <- vech2m(ParamsTB[["covariance"]])
  B <- ParamsTB[["coefficients"]]
  knots<- ParamsTB[["knots"]]
  knots.list <- knots_mat2list(ParamsTB[["knots"]], splineArgs)

  ## Pre-compute essential parts
  X <- d.matrix(x = x,knots = knots.list,splineArgs = splineArgs) # The design matrix.

  ## Only return the surface mean
  if("surface-mean" %in% callParam$id)
    {
      out <- X%*%B
      return(out)
    }

  dim.x <- dim(x)
  n <- dim.x[1] # no. of obs
  p <- dim(Y)[2] # multivariate if p > 1
  q <- dim(X)[2] # no. of covs including knots and intercept.

  diag.K.list <- lapply(apply(matrix(diag.K, p), 2, list), unlist)

  Sigma.inv <- ginv(Sigma) # inverse of Sigma
  P4X <- crossprod(X) # X'X where X is the design matrix

  q.knots <- sapply(knots.list, nrow) # no. of knots used for surface, and additive
  q.i <- c(q - sum(q.knots), q.knots) # no. covs used in each components,  cov,  surface,
                                        # additive

  ## Extract the prior settings
  P.mats.all <- P.matrix(X, q.i, priorArgs) # The P matrices and X matrices, list
  P.mats <- P.mats.all[["P"]]
  X.mats <- P.mats.all[["X"]]
  P.type <- priorArgs$P.type # The type of P matrices of the prior

  mu <- priorArgs$coefficients.mu0
  n0 <- priorArgs$covariance.df0
  S0 <- priorArgs$covariance.S0

  ## Storage
  logprior <- list()
  logpost <- list()

###----------------------------------------------------------------------------
### THE LIKELIHOOD
###----------------------------------------------------------------------------

  ## The likelihood for the Gaussian case TODO: make it general,  write a
  ## separate function if needed.
  loglik<- -n/2*determinant(2*pi*Sigma)$modulus[1] -
    1/2* tr(Sigma.inv %*% crossprod(Y-X%*%B))

  ## The full log likelihood(Can be use for calculating LPDS)
  if("likelihood" %in% callParam$id)
    {
      out <- loglik
      return(out)
    }

###----------------------------------------------------------------------------
### The priors w.r.t. differnt conditions
###----------------------------------------------------------------------------

  ## Remember to use the final scale,  since the prior are set on the final scale,
  ## e.g. the shrinkages are estimated with a log link.
  if("knots" %in% callParam$id) ## The prior for the knots
    {
      ## Get the priors parameters
      pri.type <- priorArgs$knots.priType
      pri.mean <- priorArgs$knots.mu0
      pri.covariance <- priorArgs$knots.Sigma0
      pri.shrinkage <- priorArgs$knots.c
      logprior[["knots"]] <- log_prior(B = Params[["knots"]], priorArgs = list(prior_type
                                                                = pri.type, mean =
                                                                pri.mean, covariance =
                                                                pri.covariance, shrinkage
                                                                = pri.shrinkage))
    }
  if("shrinkages" %in% callParam$id) ## The prior for the  shrinkage
    {
      pri.type <- priorArgs$shrinkages.priType
      pri.mean <- priorArgs$shrinkages.mu0
      pri.covariance <- priorArgs$shrinkages.Sigma0
      pri.shrinkage <- priorArgs$shrinkages.c

      logprior[["shrinkages"]] <- log_prior(B = Params[["shrinkages"]], priorArgs =
                                            list(prior_type = pri.type, mean = pri.mean,
                                                 covariance = pri.covariance, shrinkage =
                                                 pri.shrinkage))
    }
  if("covariance" %in% callParam$id)
    {
      pri.type <- priorArgs$shrinkages.priType
      pri.mean <- priorArgs$shrinkages.mu0
      pri.covariance <- priorArgs$shrinkages.Sigma0
      pri.shrinkage <- priorArgs$shrinkages.c

      logprior[["covariance"]] <- log_prior(B = Params[["covariance"]], priorArgs =
                                            list(prior_type = pri.type, mean = pri.mean,
                                                 covariance = pri.covariance, shrinkage =
                                                 pri.shrinkage)) ## FIXME: check!
    }

###----------------------------------------------------------------------------
### The posterior
###----------------------------------------------------------------------------

  ## The marginal posterior (without coefficients)
  logpost <- loglik + sum(unlist(logprior))

  out <- sum(unlist(logpost))

  return(out)
}

##----------------------------------------------------------------------------------------
## TESTS: PASSED
##----------------------------------------------------------------------------------------

## linear_logpost(Y, x, Params, callParam = list(id = c("knots")), splineArgs, priorArgs,
## ParamsTransArgs)
