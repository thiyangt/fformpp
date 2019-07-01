##' Diagnosis if the spline model
##'
##' This function is only applied to the moving knots model.
##' @title
##' @param criterion
##' @param OUT.Params
##' @return
##' @references
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note First version: Tue Aug 07 16:55:03 CEST 2012;
##'       Current: Tue Aug 07 16:54:52 CEST 2012.
##' @export
FitDiagnosis <- function(x, Y, OUT.Params, Data.gen, logpost.fun.name, splineArgs,
                         splineArgs.gen,  Params_Transform, crossvalid.struc, burn.in,
                         criterion = c("LOSS", "KL", "L2"))
  {
    cat("Fit of Goodness Diagnoses >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n\n")
    out <- list()

    ## Assume no cross validation TODO: Maybe will change later
    if(length(crossvalid.struc$training)>1)
      {
        return(NULL)
      }

    n <- dim(x)[1] # no. of obs.
    p <- dim(Data.gen$B)[2]

    nIter <- dim(OUT.Params[[1]])[3] # no. of iterations
    num.burn.in <- floor(nIter*burn.in)
    nUsed <- nIter - num.burn.in

    sampleProp <- 1 # TODO: Add this as a parameter if to slow.
    sampleIdx <- sort(seq(nIter, (nIter-nUsed+1), by = -round(1/sampleProp)))

    sampleIdx.len <- length(sampleIdx)

    ## The Storage of Surface mean
    SurfaceMean <- array(NA, c(n, p, sampleIdx.len))

    X.pred <- d.matrix(x = x, knots = Data.gen[["knots"]], splineArgs.gen)
    SurfaceMean.gen0 <- X.pred %*% Data.gen[["B"]]
    SurfaceMean.gen <- SurfaceMean.gen0

    out[["details"]]$SurfaceMean.true <- SurfaceMean.gen
    out[["details"]]$OUT.Params <- OUT.Params

    ## Cal. the surface mean from predicted values
    iCross <- 1
    which.j <- 0
    for(j in sampleIdx) ## Just the likelihood function with posterior samples
      {

        Params.j <- lapply(OUT.Params, function(x) apply(x[, , j, iCross, drop =
                                                           FALSE], c(1, 2), "["))
        caller.SurfaceMean <- call(logpost.fun.name,Y = NA, x = x,
                                   Params = Params.j, callParam = list(id = "surface-mean"),
                                   priorArgs = NA, splineArgs = splineArgs,
                                   Params_Transform = Params_Transform) # X*B
        SurfaceMean.j <- eval(caller.SurfaceMean)

        which.j <- which.j + 1
        SurfaceMean[, ,which.j] <- SurfaceMean.j

        ## Simple progress bar
        ## progressbar(which.j, length(sampleIdx))
      }

    SurfaceMean.mean <- apply(SurfaceMean, c(1, 2), mean)

    out[["details"]]$SurfaceMean.pred <- SurfaceMean.mean

    ##
    if("LOSS"  %in% criterion) ## Squared error loss (Gelman 2003,  p.256)
      {
        out["LOSS"] <- mean((SurfaceMean.gen - SurfaceMean.mean)^2)
        ## TODO: Consider individual LOSS for each surface, did try,  no much difference
      }
    if("KL"  %in% criterion) ## Kullback–Leibler divergence
      {
        stop("not implemented yet!")
      }
    if("L2"  %in% criterion) ## Euclidean length (L2 distance)
      {
        stop("not implemented yet!")
      }
    return(out)
  }
