##' convert the knots matrix into the knots list
##'
##' <details>
##' @title
##' @param knots.mat
##' @return
##' @references
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note First version: ; Current: .
##' @export
knots_mat2list <- function(knots.mat, splineArgs)
  {
    out <- list()
    comp <- splineArgs$comp
    knots.remain <- knots.mat

    if("thinplate.s" %in% comp)
      {
        ks <- splineArgs$thinplate.s.dim[1]
        m <- splineArgs$thinplate.s.dim[2]

        knots4s.idx <- 1:(ks*m)
        knots4s <- knots.mat[knots4s.idx]
        knots.remain <- knots.remain[-(knots4s.idx)]
        out[["thinplate.s"]] <- matrix(knots4s, ks, m, byrow = TRUE)
      }

    if("thinplate.a" %in% comp)
      {
        ka <- sum(splineArgs$thinplate.a.locate)

        knots4a.idx <- 1:ka
        knots4a <- knots.remain[knots4a.idx]

        knots.remain <- knots.remain[-(knots4a.idx)]
        out[["thinplate.a"]] <- matrix(knots4a, , 1)
      }
    return(out)
  }
##----------------------------------------------------------------------------------------
## TESTS:
##----------------------------------------------------------------------------------------
