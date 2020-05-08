##' Preform X multiply Dev xi
##'
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
