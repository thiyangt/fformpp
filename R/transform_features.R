#' Function to transform features
#'
#' @name funs
#' @rdname funs
#' logit transformation
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

#' @rdname funs
#' @title Function to transform features
#' @param feature feature variable you need to transform
#' @param transformation the transformation, "logit"
#' @return vector of transformed values of the feature
#' @export
transform.features <- function(feature, transformation){

  switch(transformation,
         logit={
           feature.transformed <- glogit(feature, min(feature), max(feature))
         },
         sqrt = {
           feature.transformed <- feature^(1/2)
         }
         )
  return(feature.transformed)

}
#' @example
#' a <- 1:10
#' transform.features(a, "logit")

