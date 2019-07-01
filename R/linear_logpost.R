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
linear_logpost <- function(Y, x0, Params, callParam, splineArgs, priorArgs, Params_Transform)
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
  X <- d.matrix(x0,knots.list,splineArgs) # The design matrix.

  ## Return the surface mean and quit
  if("surface-mean" %in% callParam$id)
    {
      out <- X%*%B
      return(out)
    }

  dim.x0 <- dim(x0)
  n <- dim.x0[1] # no. of obs
  p <- dim(Y)[2] # multivariate if p > 1
  q <- dim(X)[2] # no. of covs including knots and intercept.

  diag.K.list <- lapply(apply(matrix(diag.K, p), 2, list), unlist)

  Sigma.inv <- ginv(Sigma) # inverse of Sigma
  P4X <- crossprod(X) # X'X where X is the design matrix

  q.knots <- sapply(knots.list, nrow) # no. of knots used for surface, and additive
  q.i <- c(q - sum(q.knots), q.knots) # no. covs used in each components,  cov,  surface,
                                        # additive

  ## The prior settings
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

  ## The full log likelihood(Can be use for calculating LPDS)
  if("likelihood" %in% callParam$id)
    {
      out <- -n/2*determinant(2*pi*Sigma)$modulus[1] -
        1/2* tr(Sigma.inv %*% crossprod(Y-X%*%B))
      return(out)
    }
  else  # Use the Marginal likelihood
    {
      Sigma4beta.inv <- Sigma4betaFun(diag.K, Sigma, P.mats, inverse = TRUE)
      Sigma4beta.tilde.inv <- Sigma.inv %x% P4X + Sigma4beta.inv


      ## Check if the design matrix and the covariance matrix are singular

      if(is.singular(P4X) || is.singular(Sigma4beta.tilde.inv))
        {
          out <- NaN
          return(out)
        }

      ## Not singular, continuous
      Sigma4beta.tilde <- Matrix::solve(Sigma4beta.tilde.inv)
      beta.tilde <- Sigma4beta.tilde %*% (matrix(crossprod(X, Y) %*% Sigma.inv) +
                                          Sigma4beta.inv %*% mu)
      B.tilde <- matrix(beta.tilde, q, p)

      E.tilde <- Y-X %*% B.tilde # 2 The residual
      S.tilde <- crossprod(E.tilde)/n  # Resd' * Resd
      d <- beta.tilde - mu# 3

      ## Part 1:
      q.k <- rep(q.i, each = p)
      SumqlogDet.K <- sum(q.k*log(diag.K))

      logDet.P <- sapply(P.mats, function(x) Matrix::determinant(x)$modulus[1])
      SumplogDet.P <- sum(p*logDet.P)

      out.margi.1 <- -(SumqlogDet.K - SumplogDet.P)/2

      ## Part 2:
      out.margi.2 <- -(n+n0+p+q+1)/2*determinant(Sigma)$modulus[1]

      ## Part 3
      out.margi.3 <- (-1/2*tr(Sigma.inv%*%(n0*S0 + n*S.tilde)) -
                      1/2*Matrix::crossprod(d, Sigma4beta.inv) %*% d)

      ## Part 4:
      out.margi.4 <- -1/2*Matrix::determinant(Sigma4beta.tilde.inv)$modulus[1]

      loglike.margi <- as.matrix(out.margi.1 + out.margi.2 + out.margi.3 + out.margi.4)
    }

  ## The priors w.r.t. differnt conditions
  ## Remember to use the final scale,  since the prior are set on the final scale,
  ## e.g. the shrinkages are estimated with a log link.
  if ("knots" %in% callParam$id) ## The prior for the knots
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
  if ("shrinkages" %in% callParam$id) ## The prior for the  shrinkage
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
  if ("covariance" %in% callParam$id)
    {
      logprior[["covariance"]] <- 0  # Pre-specified.
    }

  ## The marginal posterior (without coefficients)
  logpost[["margi"]] <- as.numeric(loglike.margi + sum(unlist(logprior)))

  ## Conditional posterior for the coefficients
  if("coefficients" %in% callParam$id)
    {
      beta <- matrix(B, 1)
      Norm.Sigma0 <- Sigma4beta.tilde
      Norm.Sigma <- (Norm.Sigma0 + t(Norm.Sigma0))/2

      logpost[["coefficients"]] <- dmvnorm(x = beta, mean = beta.tilde, sigma =
                                           Norm.Sigma, log = TRUE)
    }


    out <- sum(unlist(logpost))

  return(out)
}

##----------------------------------------------------------------------------------------
## TESTS: PASSED
##----------------------------------------------------------------------------------------

## linear_logpost(Y, x, Params, callParam = list(id = c("knots")), splineArgs, priorArgs,
## ParamsTransArgs)
