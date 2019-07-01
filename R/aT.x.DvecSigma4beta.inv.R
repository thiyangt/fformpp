##' <description>
##'
##' <details>
##' @title
##' @param aT
##' @param Sigma.inv
##' @param Xmats.delta.knots
##' @param P.type
##' @param p
##' @param q
##' @param q.i
##' @return
##' @references
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note First version: Tue Feb 01 10:25:02 CET 2011;
##'       Current:       Tue Feb 01 10:25:08 CET 2011.
##' @export
aT.x.DvecSigma4beta.inv <- function(aT, Sigma.inv, diag.K.list, Xmats.delta.knots, P.type,
                                    p, q, q.i, n.par4knots)
{
    ## Convert aT to the new order for vecSigma_b.tilde See. Lemma 1
    idx4b <- idx.beta2b(p, q, q.i)
    idx4Sigmabeta <- matrix(1:(p^2*q^2), p*q)
    idx4Sigmab <- idx4Sigmabeta[as.vector(idx4b), , drop = FALSE][, as.vector(idx4b), drop = FALSE]
    aT4b <- aT[, as.vector(idx4Sigmab), drop = FALSE]
    dim.aT4b <- dim(aT4b)

    ## basic information
    n.spline.comp <- length(Xmats.delta.knots) ## number of spline components
    name.spline.comp <- names(Xmats.delta.knots) ## the name of spline components
    q.i.cumsum<- cumsum(q.i) ## The artificial index

    ## data manipulate

    ## The output matrix
    out.lst <- list()

    ## Loop over the knots components see Li's notes
    for(i in 1:n.spline.comp)
    {
        ## Check if the prior of P_i,  and have the gradient respectively.
        if(P.type[i+1] == "X'X") ## Situation that P.i  = "X'X"
        {
            ## Calculate the "A" matrix. See Li's notes.
            idx4alphaT4b <- (p^2*q*(q.i.cumsum[i])+ 1):(p^2*q*(q.i.cumsum[i+1]))
            A.i.full <- aT4b[, idx4alphaT4b, drop = FALSE]

            idx4A.i <- matrix(FALSE, p*q, p*q.i[i+1]) # Initial index
            idxnoempty <- (p*(q.i.cumsum[i])+ 1):(p*(q.i.cumsum[i+1]))
            idx4A.i[idxnoempty, ] <- TRUE # keep these cols in A.i

            A.i <- A.i.full[, as.vector(idx4A.i), drop = FALSE] # A.i taking away
            ## Calculate the last part

            KSigmaK.inv <- ((1/sqrt(diag.K.list[[i+1]])) %d*d% Sigma.inv)

            Mat1 <- Mat.x.DvecA.k.P_stp1(A.i, KSigmaK.inv, p, q, q.i[i+1])
            Mat2 <- Mat.x.DvecA.k.P_stp2(Mat1, Xmats.delta.knots[[i]], p, q, q.i[i+1])

            out.lst[[name.spline.comp[i]]] <- Mat2
        }
        else ## return zero matrix
        {
            out.lst[[name.spline.comp[i]]] <- matrix(0,  dim.aT4b[1], n.par4knots[i])
        }
    }
    ## combine the output,
    out.mat <- matrix(unlist(out.lst), dim.aT4b[1]) # dim.aT4b[1]-by-unknown

    return(out.mat)
}

##----------------------------------------------------------------------------------------
## TESTS:
##----------------------------------------------------------------------------------------
## out <- aT4b.k.B.x.DvecSigma4beta.inv(a, B, Sigma.inv, P.mats, X.mats, delta.xi, P.type, p, q, q.i)
