\name{ctmboost}
\alias{ctmboost}
\title{
  Likelihood-based Boosting for Conditional Transformation Models
}
\description{
  Employs maximisation of the likelihood for estimation of conditional
  transformation models
}
\usage{
ctmboost(model, formula, data = list(), weights = NULL, method = quote(mboost::mboost), ...)
}
\arguments{
  \item{model}{an object of class \code{mlt} as returned by
\code{\link{mlt}[mlt]}.
}
  \item{formula}{a model formula describing how the parameters of
\code{model} depend on explanatory variables, see
\code{\link[mboost]{mboost}}.
}
  \item{data}{an optional data frame of observations.
}
  \item{weights}{an optional vector of weights.
}
  \item{method}{a call to \code{\link[mboost]{mboost}},
\code{\link[mboost]{gamboost}}, or
\code{\link[mboost]{blackboost}}.
}
  \item{\dots}{additional arguments to \code{method}.
}
}
\details{
  The parameters of \code{model} depend on explanatory variables in a
  possibly structured way (see Hothorn, 2018). The number of boosting
  iterations is a hyperparameter which needs careful tuning.
}
\value{An object of class \code{ctmboost} with \code{predict} and
\code{logLik} methods.}
\references{

  Torsten Hothorn (2018). Transformation Boosting Machines.

}
\examples{

  if (require("TH.data")) {
      data("bodyfat", package = "TH.data")

  }

}
\keyword{models}
\keyword{nolinear}
