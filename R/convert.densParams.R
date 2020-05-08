##' Convert the parameters corresponding to different transformations
##'
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
