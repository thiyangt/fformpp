#' Convert an array(or matrix) to a list
#'
#' @export
array2list <- function(X, MARGIN)
{
  dim4X <- dim(X)
  if(length(dim4X) >=  3) ## Array
  {
    dim4list <- dim4X[-MARGIN]
  }
  else if(length(dim4X) == 2)
  {
    dim4list <- matrix(c(1, dim4X[2], dim4X[1], 1), 2, 2)[, MARGIN]
  }
  else
  {
    stop("X must be dimension attributed!")
  }

  out <- lapply(apply(X, MARGIN = MARGIN, list),
                function(x, dim4list) array(unlist(x), dim4list), dim4list = dim4list)
  return(out)
}
