##' <description>
##'
##' <details>
##' @title
##' @param a
##' @param B
##' @param Sigma
##' @param P.mats
##' @param delta.xi
##' @return
##' @references
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note First version: Tue Feb 01 10:25:02 CET 2011;
##'       Current:       Tue Feb 01 10:25:08 CET 2011.
##' @export
Mat.x.DvecSigma.inv.k.XTX <- function(Mat, Sigma.inv, X.delta.knots.lst, p,
                                      q, q.i)
{
    ## dim.Mat <- dim(Mat)
    ## Calculate the first part
    Mat1 <- Mat.x.DvecA.k.XTX_stp1(Mat, Sigma.inv, p, q) # PASSED
    ## DEBUGGING
    ## Mat1.0 <- (diag1(p) %x% K(q, p)  %x% diag1(q)) %*% (matrix(Sigma.inv) %x% diag1(q^2))
    ## Mat1.dbg <- Mat%*%(Mat1.0 + K.X(q, q, Mat1.0, T))

    ## Take away first p*q1 col.
    Mat2 <- Mat1[, -(1:(q*q.i[1])), drop = FALSE]

    ## convert Mat2 to a list of q elements
    dim.Mat2 <- dim(Mat2)

    Mat2.ary <- array(Mat2, c(dim.Mat2[1], q, prod(dim.Mat2)/(dim.Mat2[1]*q)))

    Mat2.lst <- array2list(Mat2.ary, 3)

    ## merge the sublist
    X.delta.knots.lst2 <- unlist(X.delta.knots.lst, recursive = FALSE)

    out.lst <- mapply("%*%", Mat2.lst, X.delta.knots.lst2, SIMPLIFY = FALSE)

    ## combine the output,
    out.mat <- matrix(unlist(out.lst), dim.Mat2[1]) # pq-by-unknown

    ## browser()

    return(out.mat)
}

##----------------------------------------------------------------------------------------
## TESTS:
##----------------------------------------------------------------------------------------
## out <- aT.k.B.x.DvecSigma4beta.inv(a, B, Sigma.inv, P.mats, X.mats, delta.xi, P.type, p, q, q.i)
