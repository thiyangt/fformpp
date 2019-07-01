##' @export
sub.hessian <- function(hessian, subset)
  {
    out <- hessian[subset, , drop = FALSE][, subset, drop = FALSE]
    return(out)
  }
