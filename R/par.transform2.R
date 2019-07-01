##' Transform the parameters from the original scale to the new scale (the final scale) wrt the linkages
##'
##' Reverse step of par.transfrom()
##' @name
##' @title
##' @param par
##' @param method
##' @return
##' @references
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note First version: ; Current: .
##' @export
par.transform2 <- function(par, method)
  {
    if(tolower(method) == "identity")
      {
        out <- par
      }
    else if(tolower(method) == "log")
      {
        out <- log(par)
      }
    return(out)
  }
