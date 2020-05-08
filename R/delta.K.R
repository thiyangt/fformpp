##' Give the gradient for a diagnoal matrix (K) w.r.t. its diagonal elements.
##'
##' Details
##'
##' @name delta.K
##' @title Gradient for a diagonal matrix.
##'
##' @param K "matrix".
##'         The diagonal matrix. p-by-p
##' @param power "scalar".
##'         The power for K to obtain the gradient for vec(K^(power)) w.r.t. diag(K)'.
##'
##'
##' @return "matrix".
##'         The gradient,  pp-by-p
##'
##' @references The notes.
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note First version: Sat Nov 06 16:31:01 CET 2010
##'       Current: 	 Mon Jan 10 18:02:02 CET 2011
##'
##' @export
delta.K <- function(diag.K,power)
{
  p <- length(diag.K)
  out <- matrix(0, p^2, p)
  delta.tmp <- power*diag.K^(power-1)
  idx <- seq(1, p^3, length.out = p)
  out[idx] <- delta.tmp
  return(out)
}
