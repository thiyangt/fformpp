##' Densities for priors.
##'
##' A collection of density for common priors.
##' 
##' @name density_prior
##' @title Densities for priors
##' 
##' @param vecB "Matrix".
##'         q-by-1, the prameters to be added with a prior.
##' @param prior_type "Charactor".
##'         Types of prior
##' @param log "Logical".
##'         IF TRUE(default), the logarithm density is returned. 
##' @param ...
##'         Addtional arguments w.r.t. different priors.
##' 
##' @return  "Matrix" p-by-1. 
##' 
##' @references g-Prior: Zellner(1986)
##' @author Feng Li, Dept. of Statistics, Stockholm University, Sweden.
##' 
##' @note First version: Wed Mar 31 10:50:41 CEST 2010;
##'       Current:       Wed Sep 15 10:36:55 CEST 2010.
##' 
##' DEPENDS: mvtnorm
##' TODO: This function is used for the conditional model,  not the conjugate case.
density_prior <- function(vecB, prior_type, log = TRUE, ...)
{
  require("mvtnorm") ## Multivariate normal and t desities
  ## Do not use "library" here.     
  if (prior_type == "g-prior") ## g-prior
    { 
      
    }
  
  
}
