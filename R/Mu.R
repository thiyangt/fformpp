## This will calculate the mean with respect to links
## E(Y) = Mu, Mu = g^-1(XB)
## See the GLM textbook for details

Mu <- function(X,B,link)
  {
    if(link == "identity")
      {
        Mu.out <- X%*%B
      }
    if(link == "inverse")
      {
        Mu.out <- 1/(X%*%B)
      }
    if(link == "log")
      {
        
      }
    if(link == "logit")
      {

      }
    return(Mu.out)
  }
