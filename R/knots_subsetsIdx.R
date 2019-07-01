##' This function is to find the locations of subsets in the knots matrix
##'
##'
##' @name
##' @title
##' @param knots.subsets "list"
##' @param k "integer"
##'         No. of knots in the knots matrix.
##' @param m  "integer"
##'         Dimension of knots in the knots matrix.
##' @return "list"
##' @references
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note First version: Sat Oct 02 15:51:01 CEST 2010;
##'       Current:       Sat Oct 02 15:51:10 CEST 2010.
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
