##' A collection of log densities for common priors.
##'
##' The parameters after "..." should be matched exactly.
##'
##' @name log_prior
##' @title logarithm density for priors
##'
##' @param B "matrix".
##'         The paramter that need to be added with a prior. The prior density is
##'         calculated conditional on B. B should be always an one-column matrix,
##' @param priorArgs "list".
##'         when
##'         priorArgs$prior_type: is set to "mvnorm", you have to provide
##'         priorArgs$mean: "matrix", the mean of parameter, mu0 should be always an
##'         one-column matrix;
##'         priorArgs$covariance: "matrix", the covariance matrix. A g-prior can be
##'         constructed by setting it to X'X, where X is the covariates matrix.;
##'         priorArgs$shrinkage: "numeric", the shrinkage for the covariance.
##'
##' @return "scalor". The log density of priors.
##'
##' @references
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note First version: Tue Mar 30 09:33:23 CEST 2010;
##'       Current:       Wed Sep 15 14:39:01 CEST 2010.
##'
##' DEPENDS: mvtnorm::dmvnorm
##' TODO:
##' @export
log_prior <- function(B, priorArgs)
{
  if (tolower(priorArgs$prior_type) == "mvnorm") # vecB ~ N(mu0, c0*Sigma0)
    {

      mean <- priorArgs$mean
      covariance <- priorArgs$covariance
      shrinkage <- priorArgs$shrinkage

      log.density <- dmvnorm(t(B), mean, covariance*shrinkage, log = TRUE)

    }# if (tolower(prior_type) == "mvnorm"

  return(log.density)
}
