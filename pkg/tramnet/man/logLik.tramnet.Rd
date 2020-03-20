
\name{logLik.tramnet}
\alias{logLik.tramnet}
\title{logLik method for class \code{"tramnet"}}
\usage{
\method{logLik}{tramnet}(object, parm = coef(object, tol = 0,
  with_baseline = TRUE), w = NULL, newdata, ...)
}
\arguments{
\item{object}{object of class \code{"tramnet"}}

\item{parm}{parameters to evaluate the log likelihood at}

\item{w}{weights}

\item{newdata}{data to evaluate the log likelihood at}

\item{...}{Additional arguments to logLik}
}
\description{
logLik method for class \code{"tramnet"}
}
\author{
Lucas Kook, Torsten Hothorn
}
\value{
returns potentially weighted (\code{w}) log-likelihood based on 
\code{object} evaluated at parameters \code{parm} and data \code{newdata}
}