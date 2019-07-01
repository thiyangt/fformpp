##' Preform X multiply Dev xi
##'
##' <details>
##' @title
##' @return
##' @references
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note First version: Thu Feb 10 11:09:46 CET 2011;
##'       Current:       Thu Feb 10 11:09:53 CET 2011.
##' @export
X.x.delta.xi <- function(X, delta.knots, q.i)
{
  knots.comp <- names(delta.knots)
  knots.comp.len <- length(knots.comp)
  q <- sum(q.i)

  out <- list()
  for(i in 1:knots.comp.len)
    {
      out.mat <- crossprod(X, delta.knots[[i]])
      out.ary <- array(out.mat, c(q, prod(dim(out.mat))/(q*q.i[i+1]) ,q.i[i+1]))
      out.lst <- array2list(out.ary, 3)
      out[[i]] <- out.lst
    }

  names(out) <- knots.comp

  return(out)
}

##----------------------------------------------------------------------------------------
## TESTS: PASSED
##----------------------------------------------------------------------------------------

## system.time(X.x.delta.xi(X.mats, delta.knots) )
