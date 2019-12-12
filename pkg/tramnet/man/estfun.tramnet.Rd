
\name{estfun.tramnet}
\alias{estfun.tramnet}
\title{estfun method for class \code{"tramnet"}}
\usage{
\method{estfun}{tramnet}(object, parm = coef(object, with_baseline =
  TRUE, tol = 0), w = NULL, newdata, ...)
}
\arguments{
\item{object}{object of class \code{"tramnet"}}

\item{parm}{parameters for evaluating the score}

\item{w}{weights}

\item{newdata}{data on which to compute the score contributions}

\item{...}{additional arguments to \code{estfun}}
}
\description{
estfun method for class \code{"tramnet"}
}
\author{
Lucas Kook
}
