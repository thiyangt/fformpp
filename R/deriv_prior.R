##' A collection of gradient for common priors.
##'
##' The parameters after "..." should be matched exactly.
##'
##' @name deriv_prior
##' @title Gradient for priors
##'
##' @param B "matrix".
##'         The paramter that need to be added with a prior. The gradient and hessian are
##'         calculated conditional on B. B should be always an one-column matrix,
##' @param priorArgs "list".
##'         priorArgs$prior_type: when prior_type is set to "mvnorm", you have to provide
##'         priorArgs$mean: "matrix", the mean of parameter, mu0 should be always an
##'         one-column matrix;
##'         priorArgs$covariance: "matrix", the covariance matrix. A g-prior can be
##'         constructed by setting it to X'X, where X is the covariates matrix.;
##'         priorArgs$shrinkage: "numeric", the shrinkage for the covariance.
##'
##' @return "list". The gradient and hessian matrix, see below.
##' \item   {gradObsPri}
##'         {"matrix". One-colunm.}
##'
##' \item   {hessObsPri}
##'         {"matrix". A squre matrix. Dimension same as prior_type$Sigma0.}
##'
##' @references
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note First version: Tue Mar 30 09:33:23 CEST 2010;
##'       Current:       Wed Sep 15 14:39:01 CEST 2010.
##' TODO:
##' @export
deriv_prior <- function(B, priorArgs, hessMethod)
{
  if (tolower(priorArgs$prior_type) == "mvnorm") # vecB ~ N(mean, shrinkage*covariance)
    {
      mean <- priorArgs$mean
      covariance <- priorArgs$covariance
      shrinkage <- priorArgs$shrinkage

      gradient.out <- (- 1/shrinkage * ginv(covariance) %*% (B-mean))  # TODO:
      ## if(is(gradient.out, "try-error")) browser()
      if(tolower(hessMethod) == "exact")
      {
          hessian.out <- - 1/shrinkage * ginv(covariance)
      }
      else
        {
            hessian.out = NA
        }

    }#

  out <- list(gradObsPri = gradient.out, hessObsPri = hessian.out)
  return(out)
}
