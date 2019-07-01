##' Compute a matrix to pre-multiply the gradient for vec beta's inverse variance
##'
##' <details>
##' @title
##' @param Mat
##' @param delta4K.list
##' @param args "list" args$subset: subset of the shrinkage
##' @return
##' @references
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note First version: Wed Jan 19 16:37:25 CET 2011;
##'       Current:       Wed Jan 19 16:37:34 CET 2011.
##' @export
Mat.x.delta.vecSigma4beta.inv <- function(Mat, delta4K.list, p, q, q.i)
  {

    ##  Make indices checking
    idx4b <- idx.beta2b(p, q, q.i)
    idx4Sigmabeta <- matrix(1:(p^2*q^2), p*q)
    idx4Sigmab <- idx4Sigmabeta[as.vector(idx4b), , drop = FALSE][, as.vector(idx4b), drop = FALSE]

    ## Reorder the Mat matrix
    Mat <- Mat[, as.vector(idx4Sigmab), drop = FALSE]

    ## p <- dim(delta4K.list[[1]])[2]
    nb0 <- length(delta4K.list) #. number of blocks used
    q0 <- sqrt(sapply(delta4K.list, nrow)/p^2) # no. row/col of the P matrices
    ## q <- sum(q0)
    pq0 <- p*q0 # no. of row/col for each block of the Sigma4beta matrix
    ## gap0 <- sum(pq0) - pq0 # The empty space between each dense blocks
    head0 <- c(0, cumsum(pq0)[1:(nb0-1)]) # how many empty space shoud be used ahead.

    ## orgnize the output matrix
    out <- matrix(0, nrow(Mat), p*nb0)

    idx1.end <- 0
    idx2.end <- 0
    for(i in 1:nb0) ## TODO: Use mapply? maybe a little faster. This part is very tricky.
      ## Don't touch
      ## it unless you know what you are doing
      {
        idx0 <- (mseq(init.seq = 1:(p*q0[i]), init.length = p*q,
                      length.out = p*q0[i]) + head0[i])
        ## partition the pre-multiply matrix

        idx1 <- (idx1.end+1):(idx1.end + p^2*q*q0[i])
        idx1.end <- idx1[length(idx1)]

        MatParti <- Mat[, idx1[idx0], drop = FALSE]

        ## output matrix
        idx2 <- (idx2.end+1):(idx2.end + p)
        idx2.end <- idx2[length(idx2)]

        ## out[, idx2] <- MatParti %*% diag.K.parti
        out[, idx2] <- MatParti %*% delta4K.list[[i]]
      }
    return(out)
  }

##----------------------------------------------------------------------------------------
## Tests:: PASSED
##----------------------------------------------------------------------------------------
## p <- 2
## q <- 10 + 100 + 100
## ## Sigma <- matrix(1, 1, 1)
## Sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
## diag.K.list <- list(c(10, 20), c(30, 40), c(50, 60))
## ## diag.K.list <- list(c(10), c(30), c(50))

## P.mats <- list(matrix(rnorm(100), 10), matrix(rnorm(10000), 100), matrix(rnorm(10000), 100))
## delta4K.list <- mapply(FUN = delta.vecPartiSigma4beta.inv, diag.K_i = diag.K.list, Sigma =
##                        list(Sigma), P_i = P.mats, SIMPLIFY = FALSE) # get the gradient
##                                         # for each part of the model w.r.t. K
## Mat <- matrix(rnorm(p*q*p^2*q^2), p*q, p^2*q^2)

## ## Rprof()
## system.time(out <- delta.vecSigma4beta.inv(Mat, delta4K.list))
## ## Rprof(NULL)
## ## summaryRprof()
