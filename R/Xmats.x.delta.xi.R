##' Preform X_i multiply Dev xi
##'
##' <details>
##' @title
##' @return
##' @references
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note First version: Thu Feb 10 11:09:46 CET 2011;
##'       Current:       Thu Feb 10 11:09:53 CET 2011.
##' @export
Xmats.x.delta.xi <- function(X.mats, delta.knots, q.i, P.type)
{
  knots.comp <- names(delta.knots)
  knots.comp.len <- length(knots.comp)

  out <- list()
  for(i in 1:knots.comp.len)
    {
      ## dim.ary[[i]] <- c(q.i[i+1], prod(dim(delta.knots[[i]]))/(q.i[i+1])^2 ,q.i[i+1])
      if(P.type[i+1] == "X'X")
        {
          out.mat <- crossprod(X.mats[[i+1]], delta.knots[[i]])
          out.ary <- array(out.mat, c(q.i[i+1], prod(dim(out.mat))/(q.i[i+1])^2 ,q.i[i+1]))
          out.lst <- array2list(out.ary, 3)
          out[[i]] <- out.lst
        }
      else
        {
          out[[i]] <- NA # don't need it
        }
    }

  ## The mapply version: not so good. since only two loops and a lot of overheads
  ## out.lst <- mapply(FUN = function(x, y, z) array2list(array(crossprod(x, y), z), 3), x =
  ##                   X.mats[2:knots.comp.len], y = delta.knots, z = dim.ary, SIMPLIFY = FALSE)


  names(out) <- knots.comp

  return(out)
}

##----------------------------------------------------------------------------------------
## TESTS: PASSED
##----------------------------------------------------------------------------------------

## system.time(X.x.delta.xi(X.mats, delta.knots) )
