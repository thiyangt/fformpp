##' Transform the parameters from the original scale to the new scale (the final scale) wrt the linkages
##'
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
