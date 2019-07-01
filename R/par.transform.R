##' Description
##'
##' Details.
##' @name
##' @title
##' @param par
##' @param method
##' @return
##' @references
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note First version: ; Current: .
##' @export
par.transform <- function(par, method)
  {
    if(tolower(method) == "identity")
      {
        out <- par
      }
    else if(tolower(method) == "log")
      {
        out <- exp(par)
      }
    return(out)
  }
