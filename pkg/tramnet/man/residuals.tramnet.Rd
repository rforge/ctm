
\name{residuals.tramnet}
\alias{residuals.tramnet}
\title{residuals method for class \code{"tramnet"}}
\usage{
\method{residuals}{tramnet}(object, parm = coef(object, tol = 0,
  with_baseline = TRUE), w = NULL, newdata, ...)
}
\arguments{
\item{object}{object of class \code{"tramnet"}}

\item{parm}{parameters to evaluate score at}

\item{w}{weights}

\item{newdata}{data to evaluate score at}

\item{...}{additional arguments to \code{residuals}}
}
\description{
residuals method for class \code{"tramnet"}
}
\author{
Lucas Kook
}
\value{
Returns a numeric vector of residuals for each row in \code{newdata}
}