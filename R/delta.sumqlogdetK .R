##' gradient for \sum q_i log|K_i|
##'
##'
##' @title
##' @param q.i
##' @param diag.K
##' @param args "list" args$subset: the subsets of the diagonal K matrix.
##' @return "matrix" 1-by-
##' @references
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note First version: Wed Jan 05 21:13:27 CET 2011;
##'       Current:       Wed Jan 05 21:13:37 CET 2011.
##' @export
delta.sumqlogdetK <- function(q.i, diag.K)
  {
    p <- length(diag.K)/length(q.i)
    q.seq <- rep(q.i, each = p)
    ## subset <- args$subset
    out0 <- q.seq/diag.K # the gradient
    out <- matrix(out0, 1) # format to a matrix
    return(out)
  }
