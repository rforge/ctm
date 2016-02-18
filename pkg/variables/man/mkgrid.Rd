\name{mkgrid}
\alias{mkgrid}
\alias{mkgrid.continuous_var}
\title{
  Generate Grids of Observations
}
\description{
  Make a grid of values
}
\usage{
mkgrid(object, ...)
\method{mkgrid}{continuous_var}(object, n = 2, ...)
}
\arguments{
  \item{object}{an object of class \code{var} or \code{vars}}
  \item{n}{number of grid points for a continous variable}
  \item{\dots}{additional arguments}
}
\details{
  The function returns a names list of values for each variable.
}
