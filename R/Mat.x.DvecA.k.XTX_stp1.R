##' Perform a dense matrix multiplying by Dev[vec[A %*% P_i]]/Dev[vecX_i'] %*% (I_q_i^2 +
##' K_qi,qi)
##'
##' Part 1.
##' @title
##' @param Mat "Matrix"
##' @param A "Matrix"
##' @return "Matrix"
##' @references Notes
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note First version: Thu Jan 27 17:41:38 CET 2011;
##'       Current:       Thu Jan 27 17:41:45 CET 2011.
##' @export
Mat.x.DvecA.k.XTX_stp1 <- function(Mat, A, p, q)
{

    dim.Mat <- dim(Mat)

    ## Mat <- diag((p*q)^2)
    ## idx4Mat <- matrix(seq(1:(p^2*q^2)), p^2*q, q, byrow = TRUE) # The auxiliary index,
                                        # see Li's notes.
    ## idx4Mat2 <- foldMat(idx4Mat, nfolds = p, byrow = TRUE) # fold the index matrix to get
                                        # the new index,  see Li's notes

    idx4Mat <- matrix(array(1:(p^2*q^2), c(q, p, p*q))
                      [, , as.vector(matrix(1:(p*q), p, q, byrow = TRUE))]
                  , , q, byrow = TRUE)

    idx4Mat2 <- foldMat(idx4Mat, nfolds = q, byrow = TRUE)

    Mat2 <- array(Mat[, as.vector(idx4Mat2)], c(dim.Mat[1], p^2, q^2)) # Reorder and
                                        # convert to an array
    A2 <- matrix(A) # Convert A to a vector

    out1.0 <- apply(X = Mat2, FUN = "%*%", MARGIN = 3, A2) # The first part
    out1 <- matrix(out1.0, dim.Mat[1]) # dimension might be dropped

    out2 <- K.X(q, q, out1, t = TRUE) ## TODO: Can K.X be a little faster?

    out <- out1 + out2
    return(out)
}
##----------------------------------------------------------------------------------------
## TESTS: PASSED
##----------------------------------------------------------------------------------------
## p <- 2
## q <- 100
## q_i <- 20
## Mat <- matrix(rnorm(p*q*p^2*q_i^2), p*q, p^2*q_i^2)
## A <- matrix(1:p^2, p)
## system.time(Mat.x.DvecA.k.P_stp1(Mat, A, p, q, q_i))
