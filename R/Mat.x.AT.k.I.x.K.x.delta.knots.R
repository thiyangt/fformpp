##' <description>
##'
##' <details>
##' @title
##' @param Mat "matrix" Sigma4beta.tilde
##' @param A
##' @param delta.knots
##' @return
##' @references
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note First version: Thu Feb 10 10:19:02 CET 2011;
##'       Current:       Thu Feb 10 10:19:08 CET 2011.
##' @export
Mat.x.AT.k.I.x.K.x.delta.knots <- function(Mat, A, delta.knots, p, q, q.i)
{

  ## browser()
  nonempty.idx.full <- matrix(1:(p*q), q, p) # The full empty index including covariates
                                        # part.

  nonempty.idx <- nonempty.idx.full[-(1:q.i[1]), , drop = FALSE]

  Mat.used <- Mat[, as.vector(t(nonempty.idx)), drop = FALSE] # use the remaining part and reorder

  Mat.ary <- array(Mat.used, c(p*q, p, sum(q.i[-1])))
  Mat.lst <- array2list(Mat.ary, 3)

  Part2.lst <- list()
  for(i in 1:length(delta.knots))
    {
      Part2.mat <- crossprod(A, delta.knots[[i]]) #  p-by-unknown   See Li' notes
      dim.Part2 <- dim(Part2.mat)
      Part2.ary <- array(Part2.mat, c(p, prod(dim.Part2)/(q.i[i+1]*p), q.i[i+1]))
      Part2.lst[[i]] <- array2list(Part2.ary, 3)
    }

  Part2.lst <- unlist(Part2.lst, recursive = FALSE)

  out.lst <- mapply("%*%", Mat.lst, Part2.lst)
  out.mat <- matrix(unlist(out.lst), p*q)

  return(out.mat)
}

##----------------------------------------------------------------------------------------
## TEST: PASSED
##----------------------------------------------------------------------------------------
## Keep this example
## n = 10
## p = 2
## q = 3

## A = matrix(1:(n*p), n, p)
## I <- diag(q)

## B <- t(A) %x% I

## C <- K.X(n, q, B, T)
## matrix(1:(p*q), q, p)
