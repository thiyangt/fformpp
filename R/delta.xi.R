##' Give the gradient for basis function w.r.t. the knots (subsets) locations for
##' different types of splines.
##'
##' The function can be extended by adding new basis. This function gives the gradient for
##' the (vecX) w.r.t. vec(xi')'. You may need the transpose of the results if you want
##' gradient for the (vecX)' w.r.t. vec(xi'). This function is mostly used inside the
##' function of gradient, so always consider the parameters input convenience
##' when you update it.
##'
##' @name delta.xi
##' @title Gradient with respect to the knots locations.
##'
##' @param x "matrix".
##'         The data matrix,  you don't have to provide a constant column.
##' @param knots "list".
##'         knots$thinplate.s:
##'         knots$thinplate.a:
##' @param knotsArgs "list".
##'         Arguments need to pass to the function, see bellow.
##'         knotsArgs$comp: "character", is the components in the design matrix;
##'         knotsArgs$thinplate.a.locate: "character",  is the location for the additive spline
##'         if used, same as d.matrix() function.
##'         knotsArgs$thinplate.s.subset: "integer", is the subset label of the knots locations
##'         to be updated in the thinplate surface component;
##'         knotsArgs$thinplate.a.subset: "integer", is the subset label of the knots locations
##'         to be updated in the thinplate additive component;
##'
##' @return "list" the gradient parts for different knots components
##'         $thinplate.s: "matrix".
##'         The gradient w.r.t the input of args$thinplate.s.subset and other
##'         relevant values. Note this will only give the dense part of the full gradient
##'         matrix w.r.t. the subset. It can also give the full gradient matrix w.r.t. all
##'         knots location if args$KnotsSubset is set to same as the full set.
##'         $thinplate.a: "matrix". Same as above but for the additive spline part.
##'
##' @references The notes.
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note First version: Fri Sep  3 17:35:47 CEST 2010.
##'       Current: 	 Fri Jan 21 13:28:21 CET 2011.
##' @export
delta.xi <- function(x, knots, splineArgs)
{

  ## There are 4 situations need to consider:
  ## Situation 1: No knots in the model.
  ## Situation 2: surface knots in the model.
  ## Situation 3: additive knots in the model.
  ## Situation 4: both surface and additive knots in the model.

  ## The intercept is not considered here.
  ## should be considered when plug-in to the gradient part

  ## The basic idea is as follows. First reorder the data set (including the design
  ## matrix and knots matrix) so that the subset to pbe gradiented placed in front of
  ## the matrix. This can make the partition calculation much easier. Secondly, obtain
  ## the dense part of the gradient. And output it with other relevant stuff.

  ## FIXME: All knots are updated currently.

  knots.name <- names(knots) # knots structure used in the model

  ## initialize the output.
  out <- list()

  if ("thinplate.s" %in% knots.name) # thinplate surface should be included.
    {
      ## KnotsSubset <- knotsArgs$thinplate.s.subset # The index for the subset of the knots
      KnotsSubset <- 1:(splineArgs$thinplate.s.dim[1]) ## This came from my old code

      xi <- knots[["thinplate.s"]] # the knots matrix for the surface

      dim.x <- dim(x)
      n <- dim.x[1]
      m <- dim.x[2]
      k <- dim(xi)[1]
      k0 <- length(KnotsSubset) # no. of knots to be used

      ## The distances matrix from each obs. to each knot.
      ## e.g. D[i, j] is the distance for ith obs. to jth knot.
      ## browser()
      log.D <- rdist(x, xi, log = TRUE)  # n-by-k TODO: Only compute the locations for
                                        # subsets which can be faster.
      ## Take out the values for the subset
      log.D.sub <- log.D[, KnotsSubset, drop = FALSE] # n-by-k0
      xi.sub <- xi[KnotsSubset, , drop = FALSE] # k0-by-m

      ## prepare for the gradient and arrange the matrix.
      ## You can think this as an array with dim of n-by-m-by-k0
      ## But you don't have to make an array in R,  since in R even array is also a vector
      ## in the ram. This can save ram spaces and speed it up.
      log.D.sub.ary <- log.D.sub[, rep(c(1:k0), each = m), drop = FALSE] # n-by-(m*k0)
      x.ary <- matrix(x, n, m*k0) # n-by-(m*k0)
      xi.sub.ary <- matrix(t(xi.sub), n, m*k0, byrow = TRUE) # n-by-(m*k0)

      gradObsDense <- -(1+2*log.D.sub.ary)*(x.ary-xi.sub.ary)
      #browser()
      #gradObsDense[is.infinite(log.D.sub.ary)] <- 0 # take care of xlogx = NaN when x = 0.
      gradObsDense[is.na(gradObsDense)] <- 0

      ## The output
      out.mat <- gradObsDense
      ## out.ary <- array(gradObsDense, c(n, m, k0))
      ## out.lst <- array2list(out.ary, 3)

      out[["thinplate.s"]] <- out.mat
    }

  if("thinplate.a" %in% knots.name) # thinplate additive should be in the model also.
    {
      xi <- knots[["thinplate.a"]] # The knots for the surface part.
      ## KnotsSubset <- knotsArgs$thinplate.a.subset # The index for the subset of the
      ## knots
      KnotsSubset <- 1:length(xi)
      xi.location <- splineArgs$thinplate.a.locate
      ## Setup storage for the dense part of the gradient.
      dim.x <- dim(x)
      n <- dim.x[1]
      m <- dim.x[2]

      idx4x <- rep(1:m, xi.location)[KnotsSubset] # The index for x be used
                                        # w.r.t. the subset knots

      x.sub <- x[, idx4x, drop = FALSE]
      xi.sub <- matrix(xi, n, length(idx4x), byrow = TRUE)
      dist.sub <- x.sub-xi.sub

      gradObsDense <- -(1+2*log(abs(dist.sub)))*dist.sub

      ## adjust the numerical error
      gradObsDense[is.na(gradObsDense)] <- 0
      #gradObsDense[is.infinite(gradObsDense)] <- 0

      out.mat <- gradObsDense
      ## out.ary <- array(gradObsDense, c(n, 1, length(xi)))
      ## out.lst <- array2list(out.ary, 3)

      out[["thinplate.a"]] <- out.mat
    }
  return(out)
}
##----------------------------------------------------------------------------------------
## TESTS: PASSED
##----------------------------------------------------------------------------------------
## x <- matrix(1:24, 6, 4)
## knots <- list(thinplate.s = matrix(1:12, 3, 4),
##               thinplate.a = rnorm(6))

## splineArgs <- list(thinplate.s.dim = c(3, 4),
##                   thinplate.a.locate = c(2, 0, 1, 3))

## delta.xi(x, knots, splineArgs)
