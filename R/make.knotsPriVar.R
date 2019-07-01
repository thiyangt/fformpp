##' set the prior variance of the knots.
##'
##' <details>
##' @title
##' @param x
##' @param splineArgs
##' @return
##' @references
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note First version: Thu Feb 17 10:22:12 CET 2011;
##'       Current:       Thu Feb 17 10:22:19 CET 2011.
##' @export
make.knotsPriVar <- function(x, splineArgs)
{
    comp <- splineArgs$comp
    out.comp <- list()
    xTx <- crossprod(x)

    xTx.inv <- ginv(xTx)

    if("thinplate.s" %in% comp)
    {
        ks <- splineArgs$thinplate.s.dim[1]
        ## out.comp[["thinplate.s"]] <- diag(ks) %x% xTx.inv

        out.comp[["thinplate.s"]] <- Diagonal(ks) %x% xTx.inv
    }
    if("thinplate.a" %in% comp)
    {
        xTx.inv.diag <- diag(xTx.inv)
        thinplate.a.locate <- splineArgs$thinplate.a.locate
        idx4a <- rep(1:length(xTx.inv.diag), thinplate.a.locate)
        prior4a <- xTx.inv.diag[idx4a]
        out.comp[["thinplate.a"]] <- diag(prior4a, length(prior4a))
    }
    out <- block.diag(out.comp)
    return(out)
}
##----------------------------------------------------------------------------------------
## TESTS: PASSED
##----------------------------------------------------------------------------------------
