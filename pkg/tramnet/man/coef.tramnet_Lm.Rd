
\name{coef.tramnet_Lm}
\alias{coef.tramnet_Lm}
\title{coef method for class \code{"tramnet_Lm"}}
\usage{
\method{coef}{tramnet_Lm}(object, with_baseline = FALSE, tol = 1e-06,
  as.lm = FALSE, ...)
}
\arguments{
\item{object}{object of class \code{"tramnet_Lm"}}

\item{with_baseline}{If \code{TRUE}, also prints coefficients
for the baseline transformation}

\item{tol}{tolerance when an estimate should be considered 0
and not returned (default: \code{1e-6})}

\item{as.lm}{If TRUE parameters are rescaled to the usual parametrization of lm}

\item{...}{Additional arguments to coef}
}
\description{
coef method for class \code{"tramnet_Lm"}
}
\examples{
\donttest{
data(cars)
m0 <- Lm(dist ~ 1, data = cars)
x <- as.matrix(cars[, "speed", drop = FALSE])
mt <- tramnet(m0, x = x, alpha = 0, lambda = 0, check_dcp = FALSE)
coef(mt)
coef(mt, with_baseline = TRUE)
coef(mt, as.lm = TRUE)
coef(lm(dist ~ speed, data = cars))
}
}
\author{
Torsten Hothorn, Lucas Kook
}
\value{Numeric vector containing the linear model shift parameter 
estimates
}
