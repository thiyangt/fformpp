##' Convert list of knots into a n-by-1 matrix.
##'
##' Converted byrow by default.
##' @title
##' @param knots.list List
##' @param byrow "logical"
##' @return Matrix
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note First version: Wed Feb 02 13:01:00 CET 2011;
##'       Current:       Wed Feb 02 13:01:07 CET 2011.
##' @export
knots_list2mat <- function(knots.list, byrow = TRUE)
  {
    knots.names <- names(knots.list)
    if(byrow == TRUE)
      {
        knots.list.new <- lapply(knots.list, FUN = t)
      }
    else
      {
        knots.list.new <- knots.list
      }
    out <- matrix(unlist(knots.list.new), , 1)
    names(out) <- NULL
    return(out)
  }
##----------------------------------------------------------------------------------------
## TESTS: PASSED
##----------------------------------------------------------------------------------------
## a <- list(thinplate.s = matrix(1:24, 6),
##           thinplate.a  = matrix(1:10))
## knots_list2mat(a)
