##' <description>
##'
##' <details>
##' @title
##' @param Mat
##' @param delta.knots
##' @param n
##' @param p
##' @param q.i
##' @return
##' @references
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note First version: ; Current: .
##' @export
Mat.delta.xi <- function(Mat, delta.knots, n, p, q.i)
{
  knots.comp.name <- names(delta.knots)
  knots.comp.len <- length(delta.knots)

  ## Remove the first n*q1 cols
  Mat2 <- Mat[, -(1:(n*q.i[1])), drop = FALSE]

  Mat2.ary <- array(Mat2, c(p^2, n, sum(q.i[-1])))
  Mat2.lst <- array2list(Mat2.ary, 3)

  ## Convert into a list
  delta.lst <- list()
  for(i in 1:knots.comp.len)
    {
      dim.delta.i <- dim(delta.knots[[i]])
      delta.i.ary <- array(delta.knots[[i]], c(n, prod(dim.delta.i)/(n*q.i[i+1]) , q.i[i+1]))
      delta.lst[[knots.comp.name[i]]] <- array2list(delta.i.ary, 3)
    }
  delta.knots.lst <- unlist(delta.lst, recursive = FALSE)

  ## the product
  out.lst <- mapply("%*%", Mat2.lst, delta.knots.lst, SIMPLIFY = FALSE)

  ## Unlist and organize the output
  out.mat <- matrix(unlist(out.lst), p^2)

##----------------------------------------------------------------------------------------
  ## DEBUGGING: PASSED
  ## xgrad <- matrix(0, n*q, 8)
  ## xgrad[(n*q.i[1]+1):(n*sum(q.i[1:2])), ] <- block.diag(list(delta.knots[[1]][, 1:4],
  ##                                                            delta.knots[[1]][, 5:8]))
  ## browser()

##----------------------------------------------------------------------------------------
  return(out.mat)
}
