\name{Aareg}
\alias{Aareg}
\title{
  Aalen Additive Hazards Model
}
\description{
  Aalen model with fully parameterised hazard function
}
\usage{
Aareg(formula, data, subset, weights, offset, cluster, na.action = na.omit, ...)
}
\arguments{
  \item{formula}{an object of class \code{"formula"}: a symbolic description 
          of the model structure to be
          fitted.  The details of model specification are given under
          \code{\link{tram}} and in the package vignette.
  }
  \item{data}{an optional data frame, list or environment (or object
          coercible by \code{as.data.frame} to a data frame) containing the
          variables in the model.  If not found in \code{data}, the
          variables are taken from \code{environment(formula)}.
  }
  \item{subset}{an optional vector specifying a subset of observations to be
          used in the fitting process.
  }
  \item{weights}{an optional vector of weights to be used in the fitting
          process.  Should be \code{NULL} or a numeric vector. If present,
          the weighted log-likelihood is maximised.
  }
  \item{offset}{this can be used to specify an _a priori_ known component to
          be included in the linear predictor during fitting.  This
          should be \code{NULL} or a numeric vector of length equal to the
          number of cases.
  }
  \item{cluster}{optional factor with a cluster ID employed for computing
          clustered covariances.
  }
  \item{na.action}{a function which indicates what should happen when the data
          contain \code{NA}s.  The default is set to \code{na.omit}.
  }
  \item{\dots}{additional arguments to \code{\link{tram}}.}
}
\details{

This function allows simultaneous
estimation of the cumulative hazard parameterised by a Bernstein polynomial. 
The model is typically fitted with time-varying coefficients, all types of random
censoring and trunction are allowed.

The responses is bounded (\code{bounds = c(0, Inf)}) when specified as a
\code{Surv} object. Otherwise, \code{bounds} can be specified via
\code{\dots}.

}
\value{
An object of class \code{Aareg}, with corresponding \code{coef},
\code{vcov}, \code{logLik}, \code{estfun}, \code{summary}, 
\code{print}, \code{plot} and \code{predict} methods.
}
\references{

  Torsten Hothorn, Lisa Moest, Peter Buehlmann (2018), Most Likely
  Transformations, \emph{Scandinavian Journal of Statistics}, \bold{45}(1),
  110--134, \doi{10.1111/sjos.12291}.

}
\examples{

  data("GBSG2", package = "TH.data")
  library("survival")
  GBSG2$time <- as.numeric(GBSG2$time)
  GBSG2$y <- with(GBSG2, Surv(time, cens))

  ### Cox proportional hazards model
  m1 <- Coxph(y ~ horTh, data = GBSG2, support = c(1, 1500))
  logLik(m1)

  ### Aalen additive hazards model with time-varying effects
  m2 <- Aareg(y | horTh ~ 1, data = GBSG2, support = c(1, 1500))
  logLik(m2)

  ### compare the hazard functions
  nd <- data.frame(horTh = unique(GBSG2$horTh))
  col <- 1:2
  lty <- 1:2
  plot(as.mlt(m1), newdata = nd, type = "hazard", 
       col = col, lty = lty[1], xlab = "time")
  plot(as.mlt(m2), newdata = nd, type = "hazard", 
       col = col, lty = 2, add = TRUE)
  legend("topright", col = rep(col, each = 2), 
         lty = rep(1:2), bty = "n",
         legend = paste(rep(paste("horTh:", 
                                  levels(nd$horTh)), each = 2),
                        rep(c("Cox", "Aalen"), 2)))

}
\keyword{models}
\keyword{regression}
\keyword{survival}
\keyword{smooth}
