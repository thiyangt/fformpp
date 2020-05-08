##' Generate the P matrices for the moving knots model
##'
##' @export
P.matrix <- function(X, q.i, priorArgs)
{
    P.type <- priorArgs$P.type
    q1 <- c(1, q.i) #for index convenient when type is "X'X"

    out.P <- list()
    out.X <- list()

    for(i in 1:length(q.i))
    {
        ## Diagonal matrix
        if(P.type[i]  == "identity")
        {
            ## out.P[[i]] <- diag1(q.i[i])

            out.P[[i]] <- Diagonal(q.i[i])

            out.X[[i]] <- NA # don't need it
        }
        else if((P.type[i]  == "X'X"))
        {
            idx <- sum(q1[1:i]):sum(q.i[1:i])
            X4P <- X[, idx, drop = FALSE] # X matrix for P
            out.P[[i]] <- crossprod(X4P)
            out.X[[i]] <- X4P
        }
        else
        {
            stop("Wrong type for the P matrices!")
        }
    }

    out <- list(P = out.P, X = out.X)
    return(out)
}

##----------------------------------------------------------------------------------------
## TESTS: PASSED
##----------------------------------------------------------------------------------------

## a <- matrix(rnorm(24), 6, 4)
## q <- c(2, 1, 1)
## P.type <- c("X'X", "X'X", "X'X")

## P.matrix(a, q, list(P.type = P.type))
