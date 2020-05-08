##' This function is to find the locations of subsets in the knots matrix
##'
##' @export
knots_subsetsIdx <- function(knots.subsets, k, m)
  {
    knots.matrix.idx <- matrix(1:(k*m), k, m)
    list.out <- list()
    for(i in 1:length(knots.subsets))
      {
        list.out[[i]] <- as.numeric(t(knots.matrix.idx[as.vector(knots.subsets[[i]]), ])) # no dim attribute
      }
    return(list.out)
  }
