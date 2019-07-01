##' Calculate the variance for the prior of coefficients (beta) in the moving knots model
##'
##' <details>
##' @title
##' @param diag.K "matrix" diagonal part of the matrix K
##' @param Sigma "matrix" variance for the error term
##' @param P.mats
##' @param inverse "logical" if TRUE,  will return the inverse variance matrix.
##' @param P "list" containing the P matrices.
##' @return "matrix"
##' @references
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note First version: Tue Feb 08 09:26:58 CET 2011;
##'       Current:       Thu Feb 24 12:50:04 CET 2011.
##' @export
Sigma4betaFun <- function(diag.K, Sigma, P.mats, inverse = FALSE)
{
  ## Step1: calculate Sigma4b,
  ## Step2: Sigma4b -> Sigma4beta

  p <- dim(Sigma)[1]
  q.n <- length(P.mats) # no. ob blocks used
  q.i <- sapply(P.mats, function(x) dim(x)[1]) # the size for each P matrics
  q <- sum(q.i)

  q1 <- c(1, p*q.i) #for index convenient
  p1 <- p*q.n
  k.init.idx <- seq(1, p1, by = p)


  out4b <- Matrix(0, p*q, p*q) # the size of the output matrix Sigma4b

  if(inverse  == FALSE) # variance
  {

      P.mats.inv <- lapply(P.mats, ginv)
      for(i in 1:q.n) # TODO: Avoid the loop
        {
          idx0 <- k.init.idx[i]:(k.init.idx[i]+p-1) # index for K matrix
          idx <- sum(q1[1:i]):sum(p*q.i[1:i]) # index for the outcome matrix
          diag.K.2 <- sqrt(diag.K[idx0])
          ## K.2 <- diag(diag.K.2, ncol = length(diag.K.2))
          Sigma.P <- (diag.K.2 %d*d% Sigma) %x% P.mats.inv[[i]]
          out4b[idx, idx] <- Sigma.P
        }
    }
  else if(inverse == TRUE) # inverse variance
    {

      for(i in 1:q.n) # TODO: Avoid the loop
        {
          idx0 <- k.init.idx[i]:(k.init.idx[i]+p-1) # index for K matrix
          idx <- sum(q1[1:i]):sum(p*q.i[1:i]) # index for the outcome matrix
          diag.K.inv2 <- 1/sqrt(diag.K[idx0])
          ## browser()
          ## K.inv2 <- diag(diag.K.inv2, ncol = length(diag.K.inv2))

          Sigma.inv <- ginv(Sigma)
          Sigma.P.inv <- (diag.K.inv2 %d*d% Sigma.inv) %x% P.mats[[i]]

          out4b[idx, idx] <- Sigma.P.inv
        }
    }


  ## Sigma4b -> Sigma4beta
  idxb2beta <- as.vector(idx.b2beta(p, q, q.i))

  out4beta <- out4b[idxb2beta, , drop = FALSE][, idxb2beta, drop = FALSE]

  return(out4beta)
}

##----------------------------------------------------------------------------------------
## Tests::PASSED
##----------------------------------------------------------------------------------------
## diag.K <- c(100, 200, 300, 400, 500, 600)
## Sigma <- matrix(c(1, 0.5, 0.5, 1), 2)
## P.mats <- list(crossprod(matrix(rnorm(20), 10)), crossprod(matrix(rnorm(30), 10)),
##                crossprod(matrix(rnorm(40), 10)))

## a = Sigma.beta(diag.K, Sigma, P.mats, inverse = FALSE)
## b = Sigma.beta(diag.K, Sigma, P.mats, inverse = TRUE)
