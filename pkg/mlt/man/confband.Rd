\name{confband}
\alias{confband}
\alias{confband.mlt}
\title{
  Confidence bands
}
\description{
  Confidence bands for transformation, distribution, survivor or
  cumulative hazard functions
}
\usage{
confband(object, newdata, level = 0.95, ...)
\method{confband}{mlt}(object, newdata, level = 0.95, 
       type = c("trafo", "distribution", "survivor", "cumhazard"), n = 20, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class \code{\link{mlt}}}
  \item{newdata}{a data frame}
  \item{level}{the confidence level}
  \item{type}{the function to compute the confidence band for}
  \item{n}{number of grid points the function is evaluated at}
  \item{\dots}{additional arguments to \code{\link[multcomp]{confint.glht}}}
}
\details{
  The function is evaluated at \code{n} grid points and simultaneous
  confidence intervals are then interpolated in order to construct the band.
}
\value{
  For each row in \code{newdata} the function and corresponding confidence
  band evaluated at the \code{n} grid points is returned.
}
