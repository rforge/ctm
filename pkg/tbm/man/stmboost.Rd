\name{stmboost}
\alias{stmboost}
\title{
  Likelihood Boosting for Shift Transformation Models  
}
\description{
  Employs maximisation of the likelihood for estimation of shift
  transformation models
}
\usage{
stmboost(model, formula, data = list(), weights = NULL, 
         method = quote(mboost::mboost), mltargs = list(), ...)
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
  \item{mltargs}{a list with arguments to be passed to
\code{\link[mlt]{mlt}}.}
  \item{\dots}{additional arguments to \code{method}.
}
}
\details{
  The parameters of \code{model} depend on explanatory variables in a
  possibly structured additive way (see Hothorn, 2020). The number of boosting
  iterations is a hyperparameter which needs careful tuning.
}
\value{An object of class \code{stmboost} with \code{predict} and
\code{logLik} methods.}
\references{

  Torsten Hothorn (2020). Transformation Boosting Machines.
  \emph{Statistics and Computing}, \bold{30}, 141--152.

}
\examples{

  if (require("TH.data") && require("tram")) {
      data("bodyfat", package = "TH.data")

      ### estimate unconditional model
      m_mlt <- BoxCox(DEXfat ~ 1, data = bodyfat, prob = c(.1, .99))
      ### get corresponding in-sample log-likelihood
      logLik(m_mlt)

      ### estimate conditional transformation model
      bm <- stmboost(m_mlt, formula = DEXfat ~ ., data = bodyfat,
                     method = quote(mboost::mboost))
      ### in-sample log-likelihood (NEEDS TUNING OF mstop!)
      logLik(bm)

      ### evaluate conditional densities for two observations
      predict(bm, newdata = bodyfat[1:2,], type = "density")
  }

}
\keyword{models}
\keyword{nonlinear}
