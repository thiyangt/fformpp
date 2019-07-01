##' @export
glm_gradhess <- function(Params, hessMethod, Y, x, callParam, splineArgs, priorArgs,
                            Params_Transform)
{
##----------------------------------------------------------------------------------------
  ## initialize and pre-computing
##----------------------------------------------------------------------------------------
  ## Transform back when Params has been transformed.
  ParamsTB <- mapply(par.transform, par = Params, method = Params_Transform, SIMPLIFY =
                     FALSE)

  ## Get the knots name
  comp <- splineArgs[["comp"]]
  knots.comp <- comp[! comp %in% c("intercept", "covariates")] # use in making the gradient

  ## Get the parameters w.r.t. the model
  diag.K <- ParamsTB[["shrinkages"]]
  Sigma <- vech2m(ParamsTB[["covariance"]])
  B <- ParamsTB[["coefficients"]]
  knots.mat <- ParamsTB[["knots"]]
  knots.list <- knots_mat2list(knots.mat, splineArgs)

  ## Pre-compute essential parts
  X <- d.matrix(x,knots.list,splineArgs) # The design matrix.

  dim.x <- dim(x)
  n <- dim.x[1] # no. of obs
  p <- dim(Y)[2] # multivariate if p > 1
  q <- dim(X)[2] # no. of covs including knots and intercept.

  diag.K.list <- lapply(apply(matrix(diag.K, p), 2, list), unlist)
  Sigma.inv <- solve(Sigma) # inverse of Sigma
  P4X <- crossprod(X) # X'X where X is the design matrix

  q.knots <- sapply(knots.list, nrow) # no. of knots used for surface, and additive
  q.i <- c(q - sum(q.knots), q.knots) # no. covs used in each components,  cov,  surface,
                                      # additive

  ## The prior settings
  P.mats.all <- P.matrix(X, q.i, priorArgs) # The P matrices and X matrices, list
  P.mats <- P.mats.all[["P"]]
  X.mats <- P.mats.all[["X"]]
  P.type <- priorArgs$P.type # The type of P matrices of the prior

  mu <- priorArgs$coefficients.mu0 # for B

  ## Boundary check
  if(knots_check_boundary(P4X, method = "singular") == "bad") # bad boundary, return NaN
                                        # and quit
    {
      out <- list(gradObs = NaN, hessObs = NaN)
      return(out)
    }
  ## Good and continuous

##----------------------------------------------------------------------------------------
  ## The gradient and hessian with respect to the knots and shrinkage
##----------------------------------------------------------------------------------------

  ## Conditional gradient for knots (surface and/or additive)
  if("knots" %in% callParam$id)
    {
##----------------------------------------------------------------------------------------
      ## gradient for likelihood part
##----------------------------------------------------------------------------------------

      ## Final gradient for marginal part.
      gradObs.margi <- t(gradObs.margi0) # transform to a col

##----------------------------------------------------------------------------------------
      ## gradient and hessian for the prior
##----------------------------------------------------------------------------------------

      ## Gradient and Hessian for prior
      pri.type <- priorArgs$knots.priType
      pri.mean <- priorArgs$knots.mu0
      pri.covariance <- priorArgs$knots.Sigma0
      pri.shrinkage <- priorArgs$knots.c

      gradHessObsPri <- deriv_prior(Params[["knots"]], priorArgs = list(mean = pri.mean,
                                                         covariance = pri.covariance,
                                                         shrinkage = pri.shrinkage,
                                                         prior_type = pri.type))

      ## Pick gradient and hessian part for the knots (subset)
      gradObs.pri <- gradHessObsPri[["gradObsPri"]][subset.idx, ,drop = FALSE]
      hessObs.pri <- sub.hessian(gradHessObsPri[["hessObsPri"]], subset.idx)
##----------------------------------------------------------------------------------------
      ## Hessian matrix for marginal part
##----------------------------------------------------------------------------------------
      if(hessMethod == "exact") # Use the exact Hessian
        {
          hessObs <- "Write the exact hessian here"
        }
      else # call the approximation of Hessian
        {
          hessObs.margi <- hessian_approx(gradient = gradObs.margi, method = hessMethod)
        }
##----------------------------------------------------------------------------------------
      ## The final gradient and Hessian
##----------------------------------------------------------------------------------------
      gradObs = gradObs.margi + gradObs.pri
      hessObs = hessObs.margi + hessObs.pri
      out <- list(gradObs = gradObs, hessObs = hessObs)
      return(out)
    }
  else if("shrinkages" %in% callParam$id) ## gradient for shrinkage K.
    {

##----------------------------------------------------------------------------------------
      ## gradient for the likelihood part
##----------------------------------------------------------------------------------------


      ## Final gradient for marginal part.
      gradObs.margi <- t(gradObs.margi0) # transform to a col

##----------------------------------------------------------------------------------------
      ## gradient and hessian for the prior
##----------------------------------------------------------------------------------------

      ## Gradient and Hessian for prior
      pri.type <- priorArgs$shrinkages.priType
      pri.mean <- priorArgs$shrinkages.mu0
      pri.covariance <- priorArgs$shrinkages.Sigma0
      pri.shrinkage <- priorArgs$shrinkages.c

      gradHessObsPri <- deriv_prior(Params[["shrinkages"]], priorArgs = list(mean = pri.mean,
                                                 covariance = pri.covariance, shrinkage =
                                                 pri.shrinkage, prior_type = pri.type))

      ## Pick gradient and hessian part for the knots (subset)
      gradObs.pri <- gradHessObsPri[["gradObsPri"]][subset.idx, ,drop = FALSE]
      hessObs.pri <- sub.hessian(gradHessObsPri[["hessObsPri"]], subset.idx)
##----------------------------------------------------------------------------------------
      ## Hessian (prior + marginal likelihood)
##----------------------------------------------------------------------------------------
      if(hessMethod == "exact") # Use the exact Hessian
        {
          hessObs <- "Write the exact hessian here"
        }
      else # call the approximation of Hessian
        {
          hessObs.margi <- hessian_approx(gradient = gradObs.margi, method = hessMethod)
        }

##----------------------------------------------------------------------------------------
      ## The final output
##----------------------------------------------------------------------------------------
      gradObs = gradObs.margi + gradObs.pri
      hessObs = hessObs.margi + hessObs.pri
      ## cat("hessObs.marig", hessObs.margi, "hessObs.pri", hessObs.pri, "\n")

      out <- list(gradObs = gradObs, hessObs = hessObs)
      return(out)
    }
  else if("covariates" %in% callParam$id) ## gradient for covariates.
    {

    }
  else if("covariance" %in% callParam$id) ## gradient for covariance.
    {

    }
  else
    {
      stop("Wrong argument for callParam !")
    }
}

## Gradient w.r.t. vecB
## Be aware that this will give a q--by--1 matrix
## Fri Mar 26 13:38:51 CET 2010
##' @export
gradient_vecB <- function(B,Sigma,x,xi,l0,l,link,gradient.prior.vecB)
  {

    p <-dim(Sigma)[1]

    X <- d.matrix(x,xi,l0,l)
    XB <- X%*%B ## Linear Predictor
    Sigma_1 <- solve(Sigma)

    mu <- Mu(X, B, link)
    residual_t<- t(Y - mu)

    grad.vecB.out <- -1/2*matrix(diag(p), nrow = 1) %*% (diag(p) %x% (Sigma_1 %*% residual_t) +
                           K.X(p,p,Sigma_1%x%residual_t,t=FALSE)) %*%
                             delta.link(X,B,link) %*% (diag(p)%x% X) + gradient.prior.vecB
    return(t(grad.vecB.out))
  }

## The gradient with respect to vech Sigma
## Mon Mar 29 09:17:32 CEST 2010
## p--by--1 matrix

##' @export
grad_vech_Sigma <- function(B,Sigma,x,xi,l0,l,link,gradient.prior.Sigma)
{
  p <-dim(Sigma)[1]
  n <- dim(x)[1]

  X <- d.matrix(x,xi,l0,l)
  XB <- X%*%B ## Linear Predictor
  Sigma_1 <- solve(Sigma)

  mu <- Mu(X, B, link)
  residual<- Y - mu

  grad.Sigma <- -n*p/2*log(2*pi)*Sigma_1 + 1/2*Sigma_1 %*% t(residual) %*% residual %*% Sigma_1 + gradient.prior.Sigma
  grad.vech.Sigma.out <- matrix(grad.Sigma[!upper.tri(grad.Sigma)],ncol=1)

  return(grad.vech.Sigma.out)
}


## gradient w.r.t. xi
##' @export
gradient_xi <- function(Y,x,xi,l0,l,n0,S0,B,ka,gradient.prior.xi)
  {

    X <- d.matrix(x,xi,l0,l)
    n <- dim(x)[1]
    p <- dim(Y)[2]
    q <- dim(X)[2]

    P <- t(X)%*%X
    P_1 <- solve(P)

    XP_1<- X%*%P # 6
    B_tilde <- 1/(1+ka)*P_1%*%(t(X)%*%Y+ka*P%*%M)
    Q_YXB <- t(Y-X%*%B_tilde) # 2
    S_tilde <- Q_YXB%*%t(Q_YXB)/n

    B_tilde_M <- B_tilde-M # 3

    S_tilde_S0 <- n0*S0+n*S_tilde+ka*t(B_tilde_M)%*%P%*%B_tilde_M #

    # grad.tmp0 <- t(delta.xi(x,xi,l0,l))
    grad.tmp0 <- t(delta.xi(x0,splineArgs))

    grad.tmp1 <- matrix(XP_1,ncol=1) # vec
    grad.tmp2 <- matrix(solve(S_tilde_S0),ncol=1) # vec

    Q_MB <- t(ka*M-(ka+1)*B_tilde) # 1
    Q_YXMXB <- t(Y+ka*X%*%M-(ka+1)*X%*%B_tilde) # 4
    Q_XYPM <- t(t(X)%*%Y-P%*%M) # 5


    Q.tmp1.0 <-  Q_MB%x%(Q_YXB%*%XP_1%*%t(X))
    Q.tmp1 <- 1/(ka+1)*(Q.tmp1.0 + K.X(p,p,Q.tmp1.0,t=FALSE))

    Q.tmp2.0 <- Q_YXMXB%x%(Q_YXB%*%XP_1)
    Q.tmp2 <- 1/(ka+1)*(K.X(n,q,(Q.tmp2.0+K.X(p,p,Q.tmp2.0,t=FALSE)),t=TRUE))

    Q.tmp3.0 <- (t(Y)-t(X%*%M))%x%t(B_tilde_M)
    Q.tmp3.1 <- t(M)%x%(t(X%*%B_tilde_M))
    Q.tmp3 <- ka/(ka+1)*(K.X(n,q,Q.tmp3.0,t=TRUE)-Q.tmp3.1)

    Q.tmp4.0 <- Q_MB%x%(Q_XYPM%*%t(XP_1))
    Q.tmp4 <- ka/(ka+1)^2*K.X(p,p,Q.tmp4.0,t=FALSE)

    Q.tmp5.0 <- (Q_XYPM%*%P_1)%x%Q_YXMXB
    Q.tmp5 <- ka/(ka+1)^2*Q.tmp5.0

    Q.tmp6.0 <- t(B_tilde)%x%Q_YXB
    Q.tmp6 <- Q.tmp6.0+K.X(p,p,Q.tmp6.0,t=FALSE)

    grad.tmpQ <- Q.tmp1 +Q.tmp2+Q.tmp3+Q.tmp4+Q.tmp5+Q.tmp6
    gradient.marginal.xi <- - grad.tmp0%*%grad.tmp1 -
                             (n+n0)/2*grad.tmp0%*%t(grad.tmpQ)%*%grad.tmp2

    gradient.xi <- gradient.marginal.xi + gradient.prior.xi

    return(gradient.xi)
  }


## Gradient w.r.t xi (conditional method)
## Be aware that this will give k--by-1 matrix
## Fri Mar 26 15:39:01 CET 2010
##' @export
gradient_xi_condi <- function(B,Sigma,x,xi,l0,l,link,gradient.prior.xi)
{
  p <-dim(Sigma)[1]
  n <- dim(x)[1]

  X <- d.matrix(x,xi,l0,l)
  XB <- X%*%B ## Linear Predictor
  Sigma_1 <- solve(Sigma)

  mu <- Mu(X, B, link)
  residual_t<- t(Y - mu)

  grad.xi.out <- -1/2*matrix(diag(p), nrow = 1) %*% (diag(p) %x% (Sigma_1 %*% residual_t) +
                                          K.X(p,p,Sigma_1%x%residual_t,t=FALSE)) %*%
                                            delta.link(X,B,link) %*% (t(B)%x% diag(n)) %*%
                                              delta.xi(x0,splineArgs) + gradient.prior.xi
    return(t(grad.xi.out))
}
