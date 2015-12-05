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
       type = c("trafo", "distribution", "survivor", "cumhazard"), 
       n = 20, ncheat = n, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class \code{\link{mlt}}}
  \item{newdata}{a data frame}
  \item{level}{the confidence level}
  \item{type}{the function to compute the confidence band for}
  \item{n}{number of grid points the function is evaluated at}
  \item{ncheat}{number of grid points the function is evaluated at when
               using the quantile obtained for \code{n} number of grid points}
  \item{\dots}{additional arguments to \code{\link[multcomp]{confint.glht}}}
}
\details{
  The function is evaluated at \code{n} grid points and simultaneous
  confidence intervals are then interpolated in order to construct the band.

  A smoother band can be obtained by setting \code{ncheat} to something larger 
  than \code{n}: The quantile is obtained for \code{n} grid points but
  the number of evaluated grid points \code{ncheat} can be much larger at no 
  additional cost. Technically, the nominal level is not maintained in
  this case but the deviation will be small for reasonably large \code{n}.
}
\value{
  For each row in \code{newdata} the function and corresponding confidence
  band evaluated at the \code{n} (or \code{ncheat}) grid points is returned.
}
