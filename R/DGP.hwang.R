##' A collection of DGPs for different models from Hwang's example.
##'
##' The model can be extended.
##' @name DGP
##' @title DATA GENERATING PROCESSES
##' @param n "numeric".
##'         No. of obs. to be generated
##' @param Sigma "matrix".
##'         The variance covariance matrix.
##' @param model "character".
##'         The model name to be use.
##'         When the model is "spline",  you also need to supply paramters,  see otherArgs.
##' @param otherArgs "list".
##'         Other arguemnts need to pass to the function.
##'         otherArgs$seed: Set the seed. If this vaule is unset, use a random
##' seed.
##'         otherArgs$xi:
##'         otherArgs$withOnt:
##'         otherArgs$spline.type:
##' @param PlotData "logical"
##'         If need to plot the data out.
##' @return "list"
##' @references Hwang 2000
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note First version: Tue Mar 30 20:11:06 CEST 2010;
##'       Current:       Tue Jul 31 13:05:02 CEST 2012.
##' DEPENDS: mvtnorm
##' @export
DGP.hwang <- function(n, q, Sigma, model, otherArgs = list(seed = NA, nTesting = NA),  PlotData)
{
  ## Get user supplied seed.
  ## If the seed is set, use it when random variable is used. Otherwise, no seed set.
  seed <- otherArgs$seed
  seed.const <- 0.24

  ## The original Hwang five examples are only bivariate-covariate (q = 2).
  ## We extend it into multivariate-covariate by allowing Euclidean distance
  ## within the covariates.

  ## n uniform draws on the unit square.
  ## The testing sample is also drawn if otherArgs$nTesting is not NULL

  nTesting <- otherArgs$nTesting

  nTotalGen <- n + nTesting
  if(!is.na(seed)) set.seed(seed)
  X.origTotal <- matrix(runif(nTotalGen*q), nTotalGen, q) # [0, 1]
  X.orig1 <- X.origTotal[ ,  1:floor(q/2), drop = FALSE]
  X.orig2 <- X.origTotal[ , (floor(q/2)+1):q, drop = FALSE]

  ## Euclidean norm to reduce to the Hwang
  enorm <- function(x) sqrt(sum(x^2))
  X.rd1 <- apply(X.orig1, 1, enorm)
  X.rd2 <- apply(X.orig2, 1, enorm)
  X.rdTotal0 <- cbind(X.rd1, X.rd2)

  ## Reduce to [0, 1]
  x201 <- function(x)
    {
      xmax <- max(x)
      xmin <- min(x)
      out <- (x -xmin)/(xmax-xmin)
      return(out)
    }
  X.rdTotal <- apply(X.rdTotal0, 2, x201)


  ## Split out the testing data
  X.orig <- X.origTotal[1:n, ] ## The training set
  X.rd <- X.rdTotal[1:n, ]

  if(!is.na(nTesting))
    {
      seqTesting <- (n+1):(n+nTesting)
      xTesting.lst <- list(xTesting = X.origTotal[seqTesting, ],
                           x.rd = X.rdTotal[seqTesting, ])
    }
  else
    {
      xTesting <- list()
    }

  ## Make the grid values for printing out the surface.
  ## FIXME: maybe replace this with expand.grid().
  ## I don't like expand.grid which gives a data frame.
  n.grid <- 20
  X1.grid <- seq(from = 0, to = 1, length.out = n.grid) # n-by-1
  X2.grid <- seq(from = 0, to = 1, length.out = n.grid) # n-by-1
  X.mesh <- mesh.grid(X1.grid, X2.grid)
  ## X1.mesh <- rep(X1.grid, each = n.grid) # n*n-by-1
  ## X2.mesh <- rep(X2.grid, times = n.grid) # n*n-by-1
  ## X.mesh <- cbind(X1.mesh, X2.mesh) # n*n-by-2
  X <- rbind(X.rd, X.mesh) # (n + n*n)-by-2
  p <- dim(Sigma)[1]

  ## Generate errors
  ## If the seed is set, use it.
  if(!is.na(seed)) set.seed(seed+seed.const)
  Errors <- rmvnorm(n = n, mean = rep(0, p), sigma = Sigma)

  SurfaceMean <- surface.hwang(model = model, X = X)

  ## The final response variable
  SurfaceMean.orig <- SurfaceMean[1:n, , drop = FALSE]
  y.orig <- SurfaceMean.orig + Errors
  names(y.orig) <- NULL

  ## Computing the nonlinear factor
  NonlinFactor <- NA

  ## Signal to noise ratio
  Errors.sd <- apply(Errors, 2, sd)
  Sig2Noise <- mean(abs(SurfaceMean.orig/Errors.sd))

  ## Print the data out if p  = 1 and q = 2
  ## FIXME: use try command to avoid error when the display is not available.
  if(PlotData == TRUE && dim(y.orig)[2] == 1 && dim(X.rd)[2] == 2)
    {
      y.grid <- matrix(SurfaceMean[(n+1):(n+n.grid^2),], n.grid, n.grid)

      image(X1.grid, X2.grid, y.grid, xlab = "X1", ylab = "X2")
      filled.contour(X1.grid, X2.grid, y.grid,
                     xlab = expression(bold(x[1])),
                     ylab  = expression(bold(x[2])),
                     main = model,
                     plot.axes = {
                       contour(X1.grid, X2.grid, y.grid,
                               add = TRUE, col = "black",
                               lwd = 0.8, lty = "solid")
                     }
                     )
      ## title(toupper(model))

      ## The 3D scatter plot
      require(rgl)
      plot3d(x = X.rd[, 1], y = X.rd[, 2], z = y.orig, col = "red",
             xlab = "X1", ylab = "X2", zlab = "y")

      ## Put the surface on
      ## FIXME: merge the two figures
      rgl.open()
      ## rgl.points(x = X.rd[, 1], y = X.rd[, 2], z = y.orig, radius = 0.01)
      rgl.surface(X1.grid, X2.grid, y.grid, front = "lines", color="#CCCCFF", alpha = .8)
      ##      spheres3d(x = X.rd[, 1], y = X.rd[, 2], z = y.orig, radius = 0.01)
    }

  out <- list(Y = y.orig, x  = X.orig, SurfaceMean = SurfaceMean.orig,
              Sigma = Sigma, NonlinFactor = NonlinFactor,
              xTesting.lst = xTesting.lst)
  return(out)
}
