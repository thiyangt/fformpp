#' logit transformation
#'
#' @param x variable name
#' @param a lower boundary for transformation
#' @param b upper boundary for transformation
#' @return a vector of logit transformed feature values
#' @export
glogit <- function(x, a, b){
  epsilon <- 1e-6
  x[x < a + epsilon] <- a + epsilon
  x[x > b - epsilon] <- b - epsilon
  p <- (x - a)/(b - a)
  return(log(p/(1-p)))
}
