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
