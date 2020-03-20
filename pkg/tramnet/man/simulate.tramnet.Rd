
\name{simulate.tramnet}
\alias{simulate.tramnet}
\title{simulate method for class \code{"tramnet"}}
\usage{
\method{simulate}{tramnet}(object, nsim = 1, seed = NULL,
  newdata = .get_tramnet_data(object), bysim = TRUE, ...)
}
\arguments{
\item{object}{object of class \code{"tramnet"}}

\item{nsim}{number of simulation}

\item{seed}{random number generator seed}

\item{newdata}{new data to simulate from}

\item{bysim}{see \code{\link[mlt]{simulate.ctm}}}

\item{...}{Additional arguments to \code{simulate.ctm}}
}
\description{
simulate method for class \code{"tramnet"}
}
\author{
Lucas Kook
}
\value{
returns a \code{list} of \code{data.frames} containing parametric bootstrap
samples based on the data supplied in \code{newdata}
}