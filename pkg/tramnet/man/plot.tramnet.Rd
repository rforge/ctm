
\name{plot.tramnet}
\alias{plot.tramnet}
\title{plot method for class \code{"tramnet"}}
\usage{
\method{plot}{tramnet}(x, newdata, type = c("distribution", "survivor",
  "density", "logdensity", "hazard", "loghazard", "cumhazard", "quantile",
  "trafo"), q = NULL, prob = 1:(K - 1)/K, K = 50, col = rgb(0.1,
  0.1, 0.1, 0.1), lty = 1, add = FALSE, ...)
}
\arguments{
\item{x}{object of class \code{"tramnet"}}

\item{newdata}{data used to predict and plot}

\item{type}{type of plot produced}

\item{q}{vector of quantiles}

\item{prob}{vector of probabilities}

\item{K}{number of data points to plot}

\item{col}{see \code{\link{plot}}}

\item{lty}{see \code{\link{plot}}}

\item{add}{see \code{\link{plot}}}

\item{...}{additional options to \code{plot}}
}
\description{
plot method for class \code{"tramnet"}
}
\author{
Lucas Kook
}
\value{None}