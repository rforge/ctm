
\name{tramnet}
\alias{tramnet}
\title{Regularised Transformation Models}
\usage{
tramnet(model, x, lambda, alpha, constraints = NULL, ...)
}
\arguments{
\item{model}{an object of class \code{"tram"} as returned by any of the
modelling functions from package tram.}

\item{x}{a numeric matrix, where each row corresponds to the same row in the
\code{data} argument used to fit \code{model}}

\item{lambda}{a positive penalty parameter for the whole penalty function}

\item{alpha}{a mixing parameter (between zero and one) defining the fraction
between absolute and quadratic penalty terms}

\item{constraints}{an optional list containing a matrix of linear inequality 
    contraints on the regression coefficients and a vector specifying the rhs
    of the inequality}

\item{\dots}{additional parameters to \code{\link[CVXR]{solve}}}
}
\value{
An object of class \code{"tramnet"} with \code{coef}, \code{logLik}, 
\code{summary}, \code{simulate}, \code{residuals} and \code{plot} methods
}
\description{
Partially penalized and constrained transformation models,
including Cox models and continuous outcome logistic regression.
The methodology is described in the \code{tramnet} vignette
accompanying this package.
}
\examples{
\donttest{
library("penalized")
library("survival")
## --- Comparison with penalized
data("nki70", package = "penalized")
nki70$resp <- with(nki70, Surv(time, event))
x <- scale(model.matrix( ~ 0 + DIAPH3 + NUSAP1 + TSPYL5 + C20orf46,
                        data = nki70))
fit <- penalized(response = resp, penalized = x, lambda1 = 1, lambda2 = 0,
                 standardize = FALSE, data = nki70)
y <- Coxph(resp ~ 1, data = nki70, order = 10, log_first = TRUE)
fit2 <- tramnet(y, x, lambda = 1, alpha = 1) ## L1 only
coef(fit)
coef(fit2)
}
}
\author{
Lucas Kook, Balint Tamasi, Sandra Sigfried
}
\keyword{model}
