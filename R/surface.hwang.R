##' @export
surface.hwang <- function(model, X, ...)
  {
    if(tolower(model) == "simple") ## only univariate response case
      {
        SurfaceMean <- 10.391*((X[,1, drop = FALSE]-0.4)*(X[,2, drop = FALSE]-0.6)+0.36)
      }
    else if(tolower(model)  == "radial") ## only univariate response case
      {
        r2 <- (X[,1, drop = FALSE]-0.5)^2+(X[,2, drop = FALSE]-0.5)^2
        SurfaceMean <- 24.234*(r2*(0.75-r2))
      }
    else if(tolower(model) == "harmonic") ## only univariate response case
      {
        XTilde <- X-0.5
        SurfaceMean <- 42.659*(0.1+XTilde[,1, drop = FALSE]
                               *(0.05+XTilde[,1, drop = FALSE]^4
                                 -10*XTilde[,1, drop = FALSE]^2*XTilde[,2, drop = FALSE]^2
                                 +5*XTilde[,2, drop = FALSE]^4))
      }
    else if(tolower(model) == "additive") ## only univariate response case
      {
        SurfaceMean <- 1.3356*(1.5*(1-X[,1, drop = FALSE])
                               + exp(2*X[,1, drop = FALSE]-1)
                               *sin(3*pi*(X[,1, drop = FALSE]-0.6)^2)
                               +exp(3*(X[,2, drop = FALSE]-0.5))
                               *sin(4*pi*(X[,2, drop = FALSE]-0.9)^2))
      }
    else if(tolower(model) == "interaction") ## only univariate response case
      {
        SurfaceMean <- 1.9*(1.35 + exp(X[,1, drop = FALSE])
                            *sin(13*(X[,1, drop = FALSE]-0.6)^2)
                            *exp(-X[,2, drop = FALSE])*sin(7*X[,2, drop = FALSE]))
      }

    return(SurfaceMean)

  }
