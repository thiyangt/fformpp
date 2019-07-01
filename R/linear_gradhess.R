##' This is the gradient and hessian matrix for the simple spline model with conjugate
##' priors.
##'
##' This function will only return the gradient and hessian part for the *labeled* knots.
##' It is possiable to provide the full gradient and hessian matrix if the parameter
##' otherArgs$subsets contains the full knots labels.
##'
##' @name linear_gradhess
##' @title Gradient and Hessian matrix for the "marginal posterior" for knots in the
##'        linear model.
##'
##' @param Params "list". Contains the matrices for the parameters.
##'         Params$knots: knots
##'         Params$shrinkages: The shrinkages
##'         Params$Sigma: The variance
##' @param hessMethod "character".
##'         Method to be used in the Hessian approximation.
##' @param Y "matrix".
##'         The response matrix, we consider the multivariate case.
##' @param x0 "matrix".
##'         The covaritates, you should *NOT* provide the intercept.
##'         It will be added automatically if necessary.
##' @param callParam
##' @param splineArgs "list".
##'         Parameters for the spline options to pass to the function, where
##'         splineArgs$method should be the spline method for x, see also d.matrix().
##'         splineArgs$method: "character". The method of splines, which can be
##'         "thinplate".
##'         splineArgs$withInt: "logical". If TRUE, will return design matrix with
##'         intercept; if FALSE, design matrix without intercept.
##' @param priorArgs
##' @param Params_Transform
##' @param callParams "character".
##'         callParams$id: the calling tag for the gradient. Possible values is "xi". "ka"
##'         callParams$subset: not all updated.
##' @param priorArgs
##'         priorArgs$prior_type: gradient type for prior
##'         priorArgs$M: mean of B,  q-by-q
##'         priorArgs$n0: df. of inverse wishart distribution.
##'         priorArgs$S0: Covariance matrix from the prior of B.
##'         priorArgs$xi.mu0: Mean of the knots locations from the prior
##'         priorArgs$xi.Sigma0: Covariance matrix from the prior of knots
##'         priorArgs$K.mu0: mean for ka
##'         priorArgs$K.Sigma0: variance for ka
##'
##' @return "list",  see bellow.
##' \item   {gradObs}
##'         {"matrix". n*p-by-1 The observed gradient for xi.}
##' \item   {hessObs}
##'         {"matrix". The obsered Hessian matrix for xi.}
##'
##' @references Li and Villani's notes
##' @author Feng Li, Dept. of Statistics, Stockholm University, Sweden.
##' @note First version: Fri Aug 27 10:36:38 CEST 2010.
##'       Current:	 Wed Sep 15 23:57:50 CEST 2010.
##'
##' DEPENDS:
##'
##' @export
linear_gradhess <- function(Params, hessMethod, Y, x0, callParam, splineArgs, priorArgs,
                            Params_Transform)
{

    ##----------------------------------------------------------------------------------------
    ## initialize and pre-computing
    ##----------------------------------------------------------------------------------------
    ## Transform back when Params has been transformed.
    ParamsTB <- mapply(par.transform, par = Params,
                       method = Params_Transform, SIMPLIFY = FALSE)

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

    ## The prior settings
    P.mats.all <- P.matrix(X, q.i, priorArgs) # The P matrices and X matrices, list
    P.mats <- P.mats.all[["P"]]
    X.mats <- P.mats.all[["X"]]
    P.type <- priorArgs$P.type # The type of P matrices of the prior

    mu <- priorArgs$coefficients.mu0 # for B
    ## browser()
    ## Boundary check
    ## if(knots_check_boundary(P4X, method = "singular") == "bad")
    ## {
    ##     ## bad boundary, return NaN and quit
    ##     out <- list(gradObs = NaN, hessObs = NaN)
    ##     return(out)
    ## }
    ## Good and continuous

    ##----------------------------------------------------------------------------------------
    ## The gradient and hessian with respect to the knots and shrinkage
    ##----------------------------------------------------------------------------------------

    ## Conditional gradient for knots (surface and/or additive)
    if("knots" %in% callParam$id)
    {
        ##----------------------------------------------------------------------------------------
        ## gradient for the marginal part
        ##----------------------------------------------------------------------------------------
        subset.idx <- as.vector(callParam[["subset"]])

        ## TODO: Consider a special case for subset If both additive and surface parts are
        ## included but one part is fixed. Don't waste time to calculate the gradient for
        ## this part since it is not used at all.

        delta.knots <- delta.xi(x0, knots.list, splineArgs)
        Xmats.delta.knots.lst <- Xmats.x.delta.xi(X.mats, delta.knots, q.i, P.type)
        X.delta.knots.lst <- X.x.delta.xi(X, delta.knots, q.i)
        n.par4knots <- sapply(delta.knots, function(x) dim(x)[2]) # no. of knots in each
                                        # components
        Sigma4beta.inv <- Sigma4betaFun(diag.K, Sigma, P.mats, inverse = TRUE)
        Sigma4beta.tilde.inv <- Sigma.inv %x% P4X + Sigma4beta.inv
        ## if(is.singular(Sigma4beta.tilde.inv))
        ##   {
        ##     out <- list(gradObs = NaN, hessObs = NaN)
        ##     return(out)
        ##   }
        ## Sigma4beta.tilde <- ginv(as.matrix(Sigma4beta.tilde.inv))

        Sigma4beta.tilde <- ginv(as.matrix(Sigma4beta.tilde.inv))

        Y.Sigma.inv <- Y %*% Sigma.inv
        XT.Y.Sigma.inv <- crossprod(X, Y.Sigma.inv)

        beta.tilde <- Sigma4beta.tilde %*% (matrix(XT.Y.Sigma.inv) +
                                            Sigma4beta.tilde.inv %*% mu)
        B.tilde <- matrix(beta.tilde, q)
        d <- beta.tilde - mu

        E.tilde <- Y - X %*% B.tilde
        E.tildeT.X <- crossprod(E.tilde, X)

        ## Part 1: PASSED
        gradObs.tmp1.0 <- delta.sumlogdetP(P.mats, q.i, Xmats.delta.knots.lst, priorArgs,
                                           n.par4knots)
        gradObs.tmp1 <- p/2*gradObs.tmp1.0

        ## Part 2: PASSED

        gradObs.tmp2.0 <- Mat.x.AT.k.I.x.K.x.delta.knots(
            Mat = Sigma4beta.tilde, A = Y.Sigma.inv, delta.knots, p, q,
            q.i) # PASSED


        if(all(mu == 0)) # gradient for beta can be simplified
        {
            aT2 <- crossprod(matrix(XT.Y.Sigma.inv), Sigma4beta.tilde)
            B2 <- Sigma4beta.tilde
            aT2.k.B2 <- aT2 %x% B2
            aTfor2 <- -aT2.k.B2
        }
        else # no way
        {
            aT1 <- t(mu) # 1-by-pq

            ## B1 <- diag1(p*q)

            B1 <- Diagonal(p*q)

            aT2 <- crossprod((matrix(XT.Y.Sigma.inv) + Sigma4beta.inv %*% mu),
                             Sigma4beta.tilde) # 1-by-p*q
            B2 <- Sigma4beta.tilde # pq-by-pq
            aT2.k.B2 <- aT2 %x% B2 ## pq-by-ppqq   TODO: HUGE object.
            aTfor2 <- Sigma4beta.tilde %*% (aT1 %x% B1) - aT2.k.B2
        }
        gradObs.tmp2.1 <- aT.x.DvecSigma4beta.inv(aTfor2, Sigma.inv, diag.K.list,
                                                  Xmats.delta.knots.lst, P.type, p, q, q.i,
                                                  n.par4knots) # 0 if P = identity,  OK for
                                        # 0.


        gradObs.tmp2.2 <- Mat.x.DvecSigma.inv.k.XTX(aT2.k.B2, Sigma.inv, X.delta.knots.lst,
                                                    p, q, q.i) ## FIXME: Slow

        gradObs.tmp2.beta.tilde <- gradObs.tmp2.0 + gradObs.tmp2.1 - gradObs.tmp2.2
                                        # Gradient for beta.tilde

        ## gradObs.tmp2.4 <- (diag1(p) %x% E.tildeT.X) %*% gradObs.tmp2.beta.tilde

        gradObs.tmp2.4 <- (Diagonal(p) %x% E.tildeT.X) %*% gradObs.tmp2.beta.tilde # dgeMatrix


        gradObs.tmp2.5.1 <- (t(B.tilde) %x% t(E.tilde))

        gradObs.tmp2.5 <- Mat.delta.xi(gradObs.tmp2.5.1, delta.knots, n, p, q.i) # FIXME: slow

        gradObs.tmp2.6 <- gradObs.tmp2.4 + gradObs.tmp2.5

        gradObs.tmp2.nS.tilde <- - gradObs.tmp2.6 - K.X(p, p, gradObs.tmp2.6, t = FALSE)

        gradObs.tmp2 <- -1/2*Matrix::crossprod(matrix(Sigma.inv),  gradObs.tmp2.nS.tilde)

        ## Part 3:
        gradObs.tmp3.1 <- 2*Matrix::crossprod(d, Sigma4beta.inv) %*% gradObs.tmp2.beta.tilde

        ddT0 <- matrix(Matrix::tcrossprod(d), 1) # row vector
        Sigma4beta.tilde0 <- matrix(Sigma4beta.tilde, 1)

        aTfor3 <- ddT0 + Sigma4beta.tilde0 # merged from Part 4.
        gradObs.tmp3.2 <- aT.x.DvecSigma4beta.inv(aTfor3, Sigma.inv, diag.K.list,
                                                  Xmats.delta.knots.lst, P.type, p, q, q.i,
                                                  n.par4knots)
        gradObs.tmp3 <- - (gradObs.tmp3.1 + gradObs.tmp3.2)/2

        ## Part 4:
        gradObs.tmp4.1<- Mat.x.DvecSigma.inv.k.XTX(Sigma4beta.tilde0, Sigma.inv,
                                                   X.delta.knots.lst, p, q, q.i)
        gradObs.tmp4 <- - 1/2*gradObs.tmp4.1

        ## The final gradObs
        gradObs.full <- gradObs.tmp1 + gradObs.tmp2 + gradObs.tmp3 + gradObs.tmp4

        ## Transform gradient according to the linkage by the chain rule
        ## only update the subset
        gradObs.orig.sub <- gradObs.full[, subset.idx, drop = FALSE]
        Param.sub <- knots.list[subset.idx]

        gradObs.logLik0 <- grad.x.deriv_link(gradObs.orig.sub, Param.sub,
                                            Params_Transform[["knots"]])

        ## Final gradient for marginal part.
        gradObs.logLik <- as.matrix(Matrix::t(gradObs.logLik0)) # transform to a col


        ##----------------------------------------------------------------------------------------
        ## gradient and hessian for the prior
        ##----------------------------------------------------------------------------------------

        ## Gradient and Hessian for prior
        pri.type <- priorArgs$knots.priType
        pri.mean <- priorArgs$knots.mu0
        pri.covariance <- priorArgs$knots.Sigma0
        pri.shrinkage <- priorArgs$knots.c

        gradHessObsPri <- deriv_prior(
            Params[["knots"]],
            priorArgs = list(mean = pri.mean,
                             covariance = pri.covariance,
                             shrinkage = pri.shrinkage,
                             prior_type = pri.type),
            hessMethod = hessMethod)

        ## Pick gradient and hessian part for the knots (subset)
        ## Pick gradient and hessian part for the knots (subset)
        gradObs.logPri <- gradHessObsPri[["gradObsPri"]][subset.idx, ,drop = FALSE]
        gradObs = gradObs.logLik + gradObs.logPri
        ##----------------------------------------------------------------------------------------
        ## Hessian (prior + marginal likelihood)
        ##----------------------------------------------------------------------------------------
        if(hessMethod == "exact") # Use the exact Hessian
        {
            ## hessObs <- "Write the exact hessian here"
            stop("Write the exact Hessian here.")
        }
        else # call the approximation of Hessian
        {
            ## hessObs.pri <- sub.hessian(gradHessObsPri[["hessObsPri"]], subset.idx)
            hessObs <- hessian_approx(gradient = gradObs, method = hessMethod)
        }

        ## Check if hessian is good.
        ## if(is.singular(hessObs))
        ## {

        ##     hessObs <- hessian_approx(gradient = gradObs, method = "identity")
        ##     # hessObs = -diag(1, nrow = length(gradObs))
        ##     warning("Singular Hessian matrix ocurred. Repace with identity Hessian.")
        ## }
        ##----------------------------------------------------------------------------------------
        ## The final gradient and Hessian
        ##----------------------------------------------------------------------------------------
        out <- list(gradObs = gradObs,
                    hessObs = hessObs,
                    gradObs.logLik = gradObs.logLik,
                    gradObs.logPri = gradObs.logPri)
        return(out)
    }
    else if("shrinkages" %in% callParam$id) ## gradient for shrinkage K.
    {

        ##----------------------------------------------------------------------------------------
        ## gradient for the marginal part
        ##----------------------------------------------------------------------------------------

        ## Essential computing
        subset.idx <- as.vector(callParam[["subset"]])

        ## P4X.inv <- solve(P4X) # inverse of X'X

        P4X.inv <- ginv(P4X)

        Sigma4beta.inv <- Sigma4betaFun(diag.K, Sigma, P.mats, inverse = TRUE)
        Sigma4beta.tilde.inv <- Sigma.inv %x% P4X + Sigma4beta.inv
        ## print(Sigma4beta.tilde.inv)

        ##      browser()
        ## Check if the shrinkages make the variance singular
        ## if(is.singular(Sigma4beta.tilde.inv))
        ##   {
        ##     out <- list(gradObs = NaN, hessObs = NaN)
        ##     return(out)
        ##   }
        ## Sigma4beta.tilde <- ginv(as.matrix(Sigma4beta.tilde.inv))

        Sigma4beta.tilde <- ginv(as.matrix(Sigma4beta.tilde.inv))

        beta.tilde <- Sigma4beta.tilde %*% (matrix(crossprod(X, Y) %*% Sigma.inv) +
                                            Sigma4beta.inv %*% mu)
        B.tilde <- matrix(beta.tilde, q, p)

        E.tilde <- Y-X %*% B.tilde # The residual
        S.tilde <- crossprod(E.tilde)/n  # Resd' * Resd
        d <- beta.tilde - mu

        ## Part 1:
        gradObs.tmp1 <- -1/2*delta.sumqlogdetK(q.i, diag.K) # PASSED

        ## Part 2:
        ## gradObs.tmp2.1 <- diag1(p) %x% crossprod(E.tilde, X)

        gradObs.tmp2.1 <- Diagonal(p) %x% crossprod(E.tilde, X)

        gradObs.tmp2.2 <- -(gradObs.tmp2.1 + K.X(p, p, gradObs.tmp2.1, t = FALSE))

        if(all(mu == 0)) # Gradient for beta can be simplified
        {
            gradObs.tmp2.3.0 <- matrix(crossprod(X, Y) %*% Sigma.inv)
            gradObs.tmp2.3.2 <- crossprod(gradObs.tmp2.3.0, Sigma4beta.tilde)
            gradObs.tmp2.3 <- -gradObs.tmp2.3.2 %x% Sigma4beta.tilde
        }
        else # bad luck.
        {
            gradObs.tmp2.3.0 <- matrix(crossprod(X, Y) %*% Sigma.inv)
            gradObs.tmp2.3.1 <- Sigma4beta.inv %*% mu
            gradObs.tmp2.3.2 <- crossprod((gradObs.tmp2.3.0+gradObs.tmp2.3.1),
                                          Sigma4beta.tilde)
            gradObs.tmp2.3.3 <- - gradObs.tmp2.3.2 %x% Sigma4beta.tilde


            gradObs.tmp2.3.4 <- Sigma4beta.tilde %*% (t(mu) %x% Diagnal(p*q))
            gradObs.tmp2.3 <- gradObs.tmp2.3.3 + gradObs.tmp2.3.4
        }

        delta4K.list <- mapply(delta.vecPartiSigma4beta.inv, diag.K_i = diag.K.list,
                               Sigma.inv = list(Sigma.inv), P_i = P.mats, SIMPLIFY = FALSE)
                                        # get the gradient for each part of the model
                                        # w.r.t. K
        delta4K.beta.tilde <- Mat.x.delta.vecSigma4beta.inv(gradObs.tmp2.3, delta4K.list, p,
                                                            q, q.i)
                                        # The gradient for tilde beta.

        gradObs.tmp2.nS.tilde <- gradObs.tmp2.2 %*% delta4K.beta.tilde
        gradObs.tmp2 <- -1/2*matrix(Sigma.inv, 1)  %*% gradObs.tmp2.nS.tilde

        ## Part 3:
        gradObs.tmp3.1 <- 2*Matrix::crossprod(d, Sigma4beta.inv) %*% delta4K.beta.tilde # notice:
                                        # second part merged to Part 4
        gradObs.tmp3 <- -1/2*gradObs.tmp3.1

        ## Part 4:
        Mat4.tmp <- Matrix::tcrossprod(d) + Sigma4beta.tilde
        Mat4 <- matrix(Mat4.tmp, 1)

        gradObs.tmp4.1 <- Mat.x.delta.vecSigma4beta.inv(Mat4, delta4K.list, p, q, q.i)
        gradObs.tmp4 <- -1/2*gradObs.tmp4.1


        ## The full gradient
        gradObs.full <- gradObs.tmp1 + gradObs.tmp2 + gradObs.tmp3 + gradObs.tmp4

        ##      print(gradObs.full)
        ## Transform gradient according to the linkage by the chain rule
        ## only update the subset
        gradObs.orig.sub <- gradObs.full[, subset.idx, drop = FALSE]
        Param.sub <- diag.K[subset.idx]

        gradObs.logLik0 <- grad.x.deriv_link(gradObs.orig.sub, Param.sub,
                                            Params_Transform[["shrinkages"]])

        ## Final gradient for marginal part.
        gradObs.logLik <- as.matrix(Matrix::t(gradObs.logLik0)) # transform to a col
                                        # traditional dense matrix.

        ##----------------------------------------------------------------------------------------
        ## gradient and hessian for the prior
        ##----------------------------------------------------------------------------------------

        ## Gradient and Hessian for prior
        pri.type <- priorArgs$shrinkages.priType
        pri.mean <- priorArgs$shrinkages.mu0
        pri.covariance <- priorArgs$shrinkages.Sigma0
        pri.shrinkage <- priorArgs$shrinkages.c

        gradHessObsPri <- deriv_prior(
            Params[["shrinkages"]],
            priorArgs = list(mean = pri.mean,
                             covariance = pri.covariance,
                             shrinkage = pri.shrinkage,
                             prior_type = pri.type),
            hessMethod =  hessMethod)

        ## Pick gradient and hessian part for the knots (subset)
        gradObs.logPri <- gradHessObsPri[["gradObsPri"]][subset.idx, ,drop = FALSE]
        gradObs = gradObs.logLik + gradObs.logPri
        ##----------------------------------------------------------------------------------------
        ## Hessian (prior + marginal likelihood)
        ##----------------------------------------------------------------------------------------
        if(hessMethod == "exact") # Use the exact Hessian
        {
            ## hessObs <- "Write the exact hessian here"
            stop("Write the exact Hessian here.")
        }
        else # call the approximation of Hessian
        {
            ## hessObs.pri <- sub.hessian(gradHessObsPri[["hessObsPri"]], subset.idx)
            hessObs <- hessian_approx(gradient = gradObs, method = hessMethod)
        }

        ## Check if hessian is good.
        ## if(is.singular(hessObs))
        ## {
        ##     hessObs <- hessian_approx(gradient = gradObs, method = "identity")
        ##     # hessObs = diag(1, nrow = length(gradObs))
        ##     warning("Singular Hessian matrix ocurred. Repace with identity Hessian.")
        ## }
        ##----------------------------------------------------------------------------------------
        ## The final output
        ##----------------------------------------------------------------------------------------
        ## cat("hessObs.marig", hessObs.margi, "hessObs.pri", hessObs.pri, "\n")

        out <- list(gradObs = gradObs,
                    hessObs = hessObs,
                    gradObs.logLik = gradObs.logLik,
                    gradObs.logPri = gradObs.logPri)
        return(out)
    }
    else
    {
        stop("Wrong argument for callParam !")
    }
}
##----------------------------------------------------------------------------------------
## TESTS:
##----------------------------------------------------------------------------------------
## linear_gradhess(Params, hessMethod = "outer", Y, x, callParam = list(id = "shrinkages", subset = 2:3), splineArgs,
##                 priorArgs, Params_Transform)
