
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
\code{estfun} method for class \code{"tramnet"} which computes
the score contributions w.r.t. each model parameter.
}
\author{
Lucas Kook
}
\value{Matrix of score contributions w.r.t. model parameters
evaluated at \code{parm}
}
