##' Convert indices from beta to b
##'
##' <details>
##' @title
##' @param p
##' @param q
##' @param q.i
##' @return
##' @references
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note First version: ; Current: .
##' @export
idx.beta2b <- function(p, q, q.i)
  {
    idx4beta <- matrix(1:(p*q), q)
    cum.q.i <- c(0, cumsum(q.i))

    out4b <- matrix(0, 1:(p*q))
    for(i in 1:length(q.i))
      {
        idx4bi <- (1 + cum.q.i[i]):cum.q.i[i+1]
        out4b[(1 + cum.q.i[i]*p):(cum.q.i[i+1]*p)] <- idx4beta[idx4bi, ]
      }
    return(out4b)
  }
