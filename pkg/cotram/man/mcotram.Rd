\name{mcotram}
\alias{mcotram}
\title{
  Multivariate Count Conditional Transformation Models
}
\description{
  A proof-of-concept implementation of multivariate conditional
transformation models for count data
}
\usage{
mcotram(..., formula = ~ 1, data, theta = NULL, diag = FALSE,
        control.outer = list(trace = FALSE), scale = FALSE)
}
\arguments{
  \item{\dots}{marginal count transformation models, one for each response
}
  \item{formula}{a model formula describing a model for the dependency
structure via the lambda parameters. The default is set to \code{~ 1} for constant lambdas.
}
  \item{data}{a data.frame
}
  \item{theta}{an optional vector of starting values
}
  \item{diag}{logical; by default, diagonal elements are constrained to 1, corresponding to \code{diag = FALSE}.
}
  \item{control.outer}{a list controlling \code{\link[alabama]{auglag}}
}
  %\item{scale}{logical; parameters are not scaled prior to optimisation by default
%}
}
\details{
  The function implements multivariate count conditional transformation models.
  The response is assumed to be a vector of counts.
}
\value{
An object of class \code{mcotram} and \code{mmlt} with \code{coef} and 
\code{predict} methods.
}
\references{

%  Nadja Klein, Torsten Hothorn, Thomas Kneib (2019),
%  Multivariate Conditional Transformation Models.
%  <arxiv:1906.03151>

}
% \examples{
%
% }
\keyword{models}
