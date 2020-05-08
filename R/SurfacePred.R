##' @export
PredSurface <- function(Y, x0, OUT.FITTED, jIter, iCross, splineArgs)
{
    Params.j <- OUT.FITTED[["Params"]]
    # Params.j <- lapply(OUT.Params, function(x)
        # apply(x[, , jIter, iCross, drop = FALSE], c(1, 2), "["))

    ## diag.K <- Params.j[["shrinkages"]]
    ## Sigma <- vech2m(Params.j[["covariance"]]) #???Error
    B <- Params.j[["coefficients"]][, , jIter, iCross]
    knots.mat <- Params.j[["knots"]]
    knots.list <- knots_mat2list(knots.mat[, , jIter, iCross], splineArgs)

    X.pred <- d.matrix(x = x0, knots = knots.list, splineArgs)

    SurfaceMean <- X.pred %*% B

    return(SurfaceMean)
}

# ### Testing Rajan.Rdata
# FITTED.file <- "~/running/tsfeature_s_moving_2_plus_a_moving_2+20180322@23.57.cb16cb.Rdata"
# load(FITTED.file)
# Y.pred <- PredSurface(Y, x0 = x, OUT.FITTED, jIter = 1, iCross = 1, splineArgs)
# head(Y.pred)
# dim(Y.pred)
