% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make.knots.R
\name{make.knots}
\alias{make.knots}
\title{Locate the knots in m-space.}
\usage{
make.knots(x, method, splineArgs)
}
\arguments{
\item{x}{"matrix".
n-by-m. Covariates matrix *without* intercept.}

\item{method}{"character".
Method to be use in the knots locating method. Currently value are "k-means",
"mahalanobis-eball" which came from Villani et al (2009) and "es" for equal
spaced sample quantile with single covariate and "random".}

\item{n.knots}{"integer".
Number of knots used.}

\item{args}{"list".
Other arguments need to pass to the function w.r.t different "method". When
method is "mahalanobis-eball", you need to provide:
args$RadiusShrink: "numeric", the radus shrinkage for the ball.}
}
\value{
"list".
        A list with knots locations for given numbers of knots.
}
\description{
Setup knots for spline models.
}
\details{
This function can be used in the initial values for splines
}
\note{
First version: Wed Mar 10 14:03:31  CET 2010;
      Current:       Thu Sep 16 13:56:37 CEST 2010.
}
\references{
Appendix C. in Villani et al (2009)
}
\author{
Feng Li, Department of Statistics, Stockholm University, Sweden.
}
