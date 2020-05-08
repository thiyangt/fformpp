##' A simple sequence generater
##'
##' @export
mseq <- function(init.seq, init.length, length.out)
{

  nrow0 <- length(init.seq)
  if(init.length<nrow0)
    {
      stop("Not correct sequence setup!")
    }
  seq0 <- matrix(init.seq, nrow0, length.out)
  seq1.tmp <- seq(from = 0, by = init.length, length.out = length.out)
  seq1 <- matrix(seq1.tmp, nrow0, length.out, byrow = TRUE)

  out <- as.vector(seq0+seq1)

  return(out)
}
