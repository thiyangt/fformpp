##' This is the gradient for a single partition of beta's variance.
##'
##' <details>
##' @title
##' @param diag.K_i
##' @param Sigma.inv
##' @param P_i
##' @return "matrix"
##' @references
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note First version: Wed Jan 12 13:45:56 CET 2011; Current: Wed Jan 12 13:46:00 CET 2011.
##' @export
delta.vecPartiSigma4beta.inv <- function(diag.K_i, Sigma.inv, P_i)
  {
    p <- dim(Sigma.inv)[1]
    q_i <- dim(P_i)[1]
    out <- matrix(0, p^2*q_i^2, p)

    ## K^(-1/2) Sigma.inv K^(-3/2)/2
    diag.K2 <- 1/sqrt(diag.K_i)
    d.K <- -diag.K_i^(-3/2)/2 # drive of diagnal K
    K.pre <- matrix(diag.K2, p, p)
    d.K.post <- matrix(d.K, p, p, byrow = TRUE)
    F_i.tmp <- K.pre * Sigma.inv * d.K.post

    idx4P <- rep(rep(1:q_i, each = p), p)
    vecC <- matrix(P_i[, idx4P], p^2*q_i^2, p)

    F0 <- matrix(0, p, p^2)
    idx4F0 <- seq(1, p^2, p+1)
    F0[, idx4F0] <- F_i.tmp
    F1 <- matrix(F0, p, p^2, byrow = TRUE) # Just re-organize to have F' + F
    Fmat <- F0 + F1 # p-by-p^2,

    ## TODO: "Fmat" matrix is very sparse, chanage this code to make it faster. this code
    ## need to reconsider when P>5
    ## Think about the subset options,  not happy with this solution

    idx4Fmat0 <- matrix(rep(1:p^3, each = q_i), p*q_i, p^2)
    idx4Fmat <- matrix(as.vector(t(idx4Fmat0)), p*q_i^2, p^2, byrow = T)
    vecFmat <- matrix(Fmat[as.vector(idx4Fmat)], p^2*q_i^2, p)

    ## The final gradien
    out <- vecC*vecFmat

    return(out)
  }

##----------------------------------------------------------------------------------------
## Tests:: PASSED
##----------------------------------------------------------------------------------------
## diag.K_i <- 1:10
## Sigma <- matrix(rnorm(100), 10, 10)
## Sigma <- Sigma %*% t(Sigma)
## P_i <- matrix(rnorm(100), 100, 100)
## system.time(qq <- delta.vecPartiSigma4beta.inv(diag.K_i, Sigma, P_i))
