\name{mmlt}
\alias{mmlt}
\title{
  Multivariate Conditional Transformation Models
}
\description{
  A proof-of-concept implementation of multivariate conditional
transformation models
}
\usage{
mmlt(..., formula, data, theta = NULL, control.outer = list(trace = FALSE), scale = FALSE)
}
\arguments{
  \item{\dots}{marginal transformation models, one for each response
}
  \item{formula}{a model formula describing a model for the dependency
structure via the lambda parameters
}
  \item{data}{a data.frame
}
  \item{theta}{an optional vector of starting values
}
  \item{control.outer}{a list controlling \code{\link[alabama]{auglag}}
}
  \item{scale}{logical; parameters are not scaled prior to optimisation by default
}
}
\details{
  The function implements multivariate conditional transformation models
  as described by Klein et al (2019). The response is assumed absolutely
  continuous at the moment, discrete versions will be added later.
}
\value{
An object of class \code{mmlt} with \code{coef} and \code{predict}
methods.
}
\references{

  Nadja Klein, Torsten Hothorn, Thomas Kneib (2019),
  Multivariate Conditional Transformation Models.
  <arxiv:1906.03151>
}
\keyword{models}
