##' Direct sample the coefficients from normal distribution.
##'
##' <details>
##' @title
##' @param Y
##' @param x
##' @param OUT.Params
##' @param crossvalid.struc
##' @param nCross
##' @param nIter
##' @param splineArgs
##' @param priorArgs
##' @param Params_Transform
##' @return
##' @references
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note First version: Fri Mar 04 20:21:33 CET 2011;
##'       Current:       Fri Mar 04 20:21:41 CET 2011.
##' TODO: Use apply?
##' @export
linear_post4coef <- function(Y, x0, OUT.Params, crossvalid.struc, nCross, nIter,
                             splineArgs, priorArgs, Params_Transform)
{
    for(jCross in 1:nCross)
    {
        Y.jCross <- Y[crossvalid.struc$training[[jCross]], , drop = FALSE]
        x.jCross <- x[crossvalid.struc$training[[jCross]], , drop = FALSE]

        for(jIter in 1:nIter) # loop nIter times
        {
            ## Pick up values from Gibbs
            knots.mat <- OUT.Params[["knots"]][, , jIter, jCross]
            diag.K <- OUT.Params[["shrinkages"]][, , jIter, jCross]
            covariance <- OUT.Params[["covariance"]][, , jIter, jCross]

            knots.mat.TB <- par.transform(par = knots.mat, method =
                                              Params_Transform[["knots"]])
            diag.K.TB <- par.transform(par = diag.K, method =
                                           Params_Transform[["shrinkages"]])
            Sigma.TB <- par.transform(par = covariance, method =
                                          Params_Transform[["covariance"]])

            Sigma <- vech2m(Sigma.TB)

            Sigma.inv <- ginv(Sigma)
            knots.list <- knots_mat2list(knots.mat, splineArgs)

            X <- d.matrix(x0,knots.list,splineArgs)
            q = dim(X)[2] # total covariates in design matrix
            q.knots <- sapply(knots.list, nrow)
            q.i <- c(q - sum(q.knots), q.knots) # size for X_0(with/without intercept), X_s,  X_a
            # q <- sum(q.i)
            p <- dim(Sigma)[1]

            P4X <- crossprod(X)
            P.mats.all <- P.matrix(X, q.i, priorArgs) # The P matrices and X matrices, list
            P.mats <- P.mats.all[["P"]]

            mu <- priorArgs$coefficients.mu0

            ## K.tmp <- diag(diag.K.tmp, length(diag.K.tmp))


            Sigma4beta.inv <- Sigma4betaFun(diag.K.TB, Sigma, P.mats, inverse = TRUE)
            Sigma4beta.tilde.inv <- Sigma.inv %x% P4X + Sigma4beta.inv

            Sigma4beta.tilde <- ginv(as.matrix(Sigma4beta.tilde.inv))


            Y.Sigma.inv <- Y %*% Sigma.inv
            XT.Y.Sigma.inv <- crossprod(X, Y.Sigma.inv)

            beta.tilde <- Sigma4beta.tilde %*% (matrix(XT.Y.Sigma.inv) +
                                                    Sigma4beta.tilde.inv %*% mu)

            ## Sigma4beta.tilde should be symmetric. But might be not numerically after 1e-5
            Norm.Sigma <- (Sigma4beta.tilde + t(Sigma4beta.tilde))/2

            out.beta <- rmvnorm(mean = beta.tilde, n = 1, sigma = Norm.Sigma)

            ## Save to OUT.Params
            OUT.Params[["coefficients"]][, , jIter, jCross] <- out.beta

            ## Simple progress bar
            ## progressbar(((jCross-1)*nIter + jIter), nCross*nIter)

        }
    }
    return(OUT.Params)
}
