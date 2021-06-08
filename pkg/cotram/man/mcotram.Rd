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
mcotram(..., formula = ~ 1, data, theta = NULL,
        control.outer = list(trace = FALSE), tol = sqrt(.Machine$double.eps), 
        dofit = TRUE)
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
  \item{control.outer}{a list controlling \code{\link[alabama]{auglag}}
}
  \item{tol}{tolerance}
  \item{dofit}{logical; parameters are fitted by default, otherwise a list
               with log-likelihood and score function is returned
}
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

  Nadja Klein, Torsten Hothorn, Luisa Barbanti, Thomas Kneib (2020),
  Multivariate Conditional Transformation Models. \emph{Scandinavian Journal
  of Statistics}, \doi{10.1111/sjos.12501}.

}
% \examples{
%
% }
\keyword{models}
