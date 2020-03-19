
\name{coef.tramnet}
\alias{coef.tramnet}
\title{coef method for class \code{"tramnet"}}
\usage{
\method{coef}{tramnet}(object, with_baseline = FALSE, tol = 1e-06, ...)
}
\arguments{
\item{object}{object of class \code{"tramnet"}}

\item{with_baseline}{If \code{TRUE}, also prints coefficients
for the baseline transformation}

\item{tol}{tolerance when an estimate should be considered 0
and not returned (default: \code{1e-6})}

\item{...}{Additional arguments to coef}
}
\description{
coef method for class \code{"tramnet"}
}
\author{
Torsten Hothorn, Lucas Kook
}
\value{Numeric vector containing the model shift 
parameter estimates
}
