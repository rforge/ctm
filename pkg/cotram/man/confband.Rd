\name{confband}
\alias{confband}
\alias{confband.cotram}
\title{
  Confidence Bands
}
\description{
  Confidence bands for transformation, distribution, survivor or
  cumulative hazard functions
}
\usage{
confband(object, newdata, level = 0.95, ...)
\method{confband}{cotram}(object, newdata,  level = 0.95, 
       type = c("trafo", "distribution", "survivor", "cumhazard"), 
       smooth = FALSE, q = NULL, K = 20, cheat = K, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class \code{\link{cotram}}.}
  \item{newdata}{a data frame of observations.}
  \item{level}{the confidence level.}
  \item{type}{the function to compute the confidence band for.}
  \item{smooth}{logical; if \code{TRUE} a smoothed function of
              \code{type} is returned.}
  \item{q}{quantiles at which to evaluate the model.}
  \item{K}{number of grid points the function is evaluated at
          (in the absence of \code{q} and \code{smooth = TRUE}).}
  \item{cheat}{number of grid points the function is evaluated at when
               using the quantile obtained for \code{K} grid points
               (in the absence of \code{q} and \code{smooth = TRUE}).}
  \item{\dots}{additional arguments to \code{\link[multcomp]{confint.glht}}.}
}
\details{
  The function is evaluated at the count response or at \code{K} grid points
  and simultaneous confidence intervals are then interpolated in order to
  construct the band.
}
\value{
  For each row in \code{newdata} the function and corresponding confidence
  band evaluated at the count response (or \code{K} or \code{cheat} grid points)
  is returned.
}
\examples{

  data("birds", package = "TH.data")
  
  ### fit count transformation model with cloglog link
  m_birds <- cotram(SG5 ~ AOT + AFS + GST + DBH + DWC + LOG, data = birds,
                    method = "cloglog")
  
  ### compute asymptotic confidence bands for the distribution function
  ### for the first oberservation
  confband(m_birds, newdata = birds[1, ], type = "distribution")

}
