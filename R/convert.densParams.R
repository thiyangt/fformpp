##' Convert the parameters corresponding to different transformations
##'
##' Details.
##' @name
##' @title
##' @param mean
##' @param var
##' @return
##' @references
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note First version: Wed Feb 02 16:09:43 CET 2011;
##'       Current: Wed Feb 02 16:09:48 CET 2011.
##' @export
convert.densParams <- function(mean, var, linkage)
  {
    if(linkage == "log")
      {

        s2 <- log(var/mean^2+1)
        m <- log(mean) - s2/2
        out <- c(m, s2)
      }
    else if (linkage  == "identity")
      {
        out <- c(mean, var)
      }
    return(out)
  }
