#' Vech function
#'
#' Details.
#' @title calculating a vector
#' @param X data matrix
#' @param diag diagonal values
#' @return out a matrix
#' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
#' @export
vech <- function(X, diag = TRUE)
{
  out <- matrix(X[as.vector(lower.tri(X, diag = diag))])
  return(out)
}
