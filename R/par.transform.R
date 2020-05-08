##' Description
##'
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
