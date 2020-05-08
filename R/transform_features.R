#' @title Function to transform features
#' @param feature feature variable you need to transform
#' @param transformation the transformation, "logit"
#' @return vector of transformed values of the feature
#' @export
transform_features <- function(feature, transformation){

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
#' transform_features(a, "logit")

