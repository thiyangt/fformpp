##' D[vecf(g(x)]/D[vecg(x)'] %*% D[vec g(x)]/D[vec x']
##'
##' Just the chain rule
##' @title
##' @param grad "matrix" gradient
##' @param Param
##' @param link
##' @param par "matrix" should always be a matrix
##' @return
##' @references
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note First version: Sun Feb 20 15:03:34 CET 2011;
##'       Current:       Sun Feb 20 15:03:43 CET 2011.
##' @export
grad.x.deriv_link <- function(grad, Param, link)
  {
    if(tolower(link) == "identity")
      {
        out <- grad # unchanged
      }
    else if(tolower(link) == "log")
      {
        dim.grad <- dim(grad)

        chain2 <- matrix(Param,dim.grad[1],dim.grad[2], byrow = TRUE)
        out <- grad*chain2
      }
    return(out)
  }
