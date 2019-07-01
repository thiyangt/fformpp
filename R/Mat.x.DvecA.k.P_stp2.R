##' <description>
##'
##' <details>
##' @title
##' @param Mat "matrix"
##' @param X_i "matrix" the matrix used for obtaining P_i matrix
##' @param delta_i "matrix" the matrix for the gradient
##' @return "matrix"
##' @references
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note First version: Fri Jan 28 11:30:57 CET 2011;
##'       Current:       Fri Jan 28 11:31:06 CET 2011.
##' @export
Mat.x.DvecA.k.P_stp2 <- function(Mat, Xmats.delta.knots.i, p, q, q_i)
  {
    dim.Mat <- dim(Mat)

    ## First part
    Mat.array <- array(Mat, c(dim.Mat[1], q_i, q_i))
    Mat.list <- array2list(Mat.array, 3)

    ##  (I %x% X') D[vecX]/D knots'
    ## delta4X <- crossprod(X_i, delta.knots.i)
    ## delta4X <- Xmats.delta.knots.i

    ## dim.delta4X <- dim(delta4X)
    ## delta4X.array <- array(delta4X, c(q_i, prod(dim.delta4X)/q_i^2, q_i))
    ## FIXME: Use "Matrix" class
    delta4X.list <- Xmats.delta.knots.i
    out.tmp <- mapply("%*%", Mat.list, delta4X.list)
    out <- matrix(out.tmp, dim.Mat[1]) ## the matrix format
    return(out)
  }

##----------------------------------------------------------------------------------------
## TESTS: PASSED
##----------------------------------------------------------------------------------------
## p <- 2
## q <- 20
## q_i <- 10
## n <- 5000

## Mat <- matrix(1:(p*q*q_i^2), p*q, q_i^2)
## X_i <- matrix(rnorm(n*q_i), n, q_i)
## delta_i <- matrix(rnorm(n*q_i*q_i), n, q_i^2)

## system.time(A <- Mat.x.DvecA.k.P_stp2(Mat, X_i, delta_i, p, q, q_i))
