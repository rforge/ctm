\name{tbm}
\alias{tbm}
\title{
  Transformation Boosting Machines
}
\description{
  Likelihood-based Boosting for Transformation Models
}
\usage{
tbm(model, formula, data = list(), weights = NULL, 
    gradient = c("ctm", "shift"), baselearner = c("bbs", "bols"), ...)
}
\arguments{
  \item{model}{
  An object of class \code{\link[mlt]{mlt}}.
}
  \item{formula}{
  A model formula.
}
  \item{data}{
  An optional data frame.
}
  \item{weights}{
  An optional vector of weights.
}
  \item{gradient}{
  A specification of the gradients to be fitted; either the
  gradient for the baseline transformation parameters (\code{ctm})
  or the gradient for a scalar shift parameter.
}
  \item{baselearner}{
  A specification of the baselearners to be used.
}
  \item{\dots}{
  Additional arguments to \code{\link[mboost]{mboost}}.
}
}
\details{
  Tbd.
}
\value{
  An object of class \code{tbm} with corresponding \code{predict},
  \code{logLik} and \code{coef} methods.
}
\keyword{regression}
