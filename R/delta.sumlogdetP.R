##' <description>
##'
##' <details>
##' @title
##' @param P.mats
##' @param q.i
##' @param X.delta.knots.lst
##' @param priorArgs
##' @return
##' @references
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note First version: ; Current: .
##' @export
delta.sumlogdetP <- function(P.mats, q.i, Xmats.delta.knots.lst, priorArgs, n.par4knots)
  {
    knots.comp <- names(Xmats.delta.knots.lst)
    P.type <- priorArgs[["P.type"]]

    out <- NULL
    for(i in 1:length(knots.comp))
      {
        if(P.type[[i+1]] == "X'X")
          {

           vecP.mats.inv <- matrix(ginv(P.mats[[i+1]]), 1)
            q0 <- q.i[i+1]
            out1 <- vecP.mats.inv + K.X(q0, q0, vecP.mats.inv, t = TRUE)
            out1.ary <- array(out1, c(1, q0, q0))
            out1.lst <- array2list(out1.ary, 3)
            out2 <- mapply("%*%", out1.lst, Xmats.delta.knots.lst[[i]])
            out.current <- matrix(out2, 1)
          }
        else # zeros
          {
            out.current <- matrix(0, 1, n.par4knots[i])
          }
        out <- cbind(out, out.current)
      }

    return(out)
  }

##----------------------------------------------------------------------------------------
## TESTS: PASSED
##----------------------------------------------------------------------------------------
## delta.sumlogdetP(P.mats, q.i, X.delta.knots.lst, priorArgs)
