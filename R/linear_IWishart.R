##' Calculate the posterior df and location matrix V from the inverse Wishart distribution.
##'
##' <details>
##' @title
##' @return
##' @references
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note First version: Tue Feb 22 19:49:35 CET 2011;
##'       Current:       Tue Feb 22 19:49:42 CET 2011.
##' @export
linear_IWishart <- function(Params, Y, x0, splineArgs, priorArgs, Params_Transform)
{
  ## Transform back when Params has been transformed.
  ParamsTB <- mapply(par.transform, par = Params, method = Params_Transform, SIMPLIFY =
                     FALSE)

  ## Get the knots name
  ## comp <- splineArgs[["comp"]]
  ## knots.comp <- comp[! comp %in% c("intercept", "covariates")] # use in making the
  ## gradient

  ## Get the parameters w.r.t. the model
  diag.K <- ParamsTB[["shrinkages"]]
  Sigma <- vech2m(ParamsTB[["covariance"]])
  ## B <- ParamsTB[["coefficients"]]
  knots.mat <- ParamsTB[["knots"]]
  knots.list <- knots_mat2list(knots.mat, splineArgs)

  ## Pre-compute essential parts
  X <- d.matrix(x0,knots.list,splineArgs) # The design matrix.

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
  cum.q.i <- c(0, cumsum(q.i))

  ## The prior settings
  P.mats.all <- P.matrix(X, q.i, priorArgs) # The P matrices and X matrices, list
  P.mats <- P.mats.all[["P"]]
  ## X.mats <- P.mats.all[["X"]]
  P.type <- priorArgs$P.type # The type of P matrices of the prior

  mu <- priorArgs$coefficients.mu0
  n0 <- priorArgs$covariance.df0
  S0 <- priorArgs$covariance.S0

  Sigma4beta.inv <- Sigma4betaFun(diag.K, Sigma, P.mats, inverse = TRUE)
  Sigma4beta.tilde.inv <- Sigma.inv %x% P4X + Sigma4beta.inv

  if(is.singular(Sigma4beta.tilde.inv))
    {
      out <- list(df = NaN, V = NaN)
      return(out)
    }

  Sigma4beta.tilde <- ginv(as.matrix(Sigma4beta.tilde.inv))

  beta.tilde <- Sigma4beta.tilde %*% (matrix(crossprod(X, Y) %*% Sigma.inv) +
                                      Sigma4beta.tilde.inv %*% mu)
  B.tilde <- matrix(beta.tilde, q)
  M <- matrix(mu, q)
  D <- B.tilde-M

  E.tilde <- Y-X %*% B.tilde # 2 The residual
  S.tilde <- crossprod(E.tilde)/n  # Resd' * Resd

  ## Loop to obtain G.tilde
  G.tilde <- 0
  for(i in 1:length(q.i))
    {
      K.i0 <- diag.K.list[[i]]
      ## K.i2 <- diag(1/sqrt(K.i0), length(K.i0))
      diag.K.i2 <- 1/sqrt(K.i0)

      P.i <- P.mats[[i]]
      idx4D <- (cum.q.i[i]+1):(cum.q.i[i+1])
      D.i <- D[idx4D, ,drop = FALSE]
      if(P.type[i] == "identity") # Sparse
          {
            G.tilde.i <- diag.K.i2 %d*% crossprod(D.i) %*d% diag.K.i2
          }
        else # not sparse
          {
            G.tilde.i <- diag.K.i2 %d*% crossprod(D.i, P.i) %*% D.i %*d% diag.K.i2
          }
      G.tilde <- G.tilde + G.tilde.i
    }

  ## The posterior parameters
  V <- n0*S0 + n * S.tilde + G.tilde
  df <- n + n0

  ## browser()

  out <- list(df = df, V = V)
  return(out)
}

##----------------------------------------------------------------------------------------
## TESTS: PASSED
##----------------------------------------------------------------------------------------
## a <- linear_IWishart(Params, Y, x, splineArgs, priorArgs, Params_Transform)
