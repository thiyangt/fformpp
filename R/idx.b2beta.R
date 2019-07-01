##' Make indices from b to beta.
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
idx.b2beta <- function(p, q, q.i)
  {
    idx4b <- 1:(p*q) ## The original indices for b
    cumidx <- c(0, cumsum(q.i))

    matidx4bi <- matrix(0, q, p)
    for(i in 1:length(q.i))
      {
        idx4bi <- (1+cumidx[i]*p):(cumidx[i+1]*p)

        matidx4bi[(1+cumidx[i]):(cumidx[i+1]), ] <- matrix(idx4b[idx4bi], q.i[i])
      }
    idx4beta <- as.vector(matidx4bi)
    return(idx4beta)
  }
