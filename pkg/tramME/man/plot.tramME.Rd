% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.R
\name{plot.tramME}
\alias{plot.tramME}
\title{Plotting method for tramME objects}
\usage{
\method{plot}{tramME}(x, newdata = NULL, ranef = NULL, ...)
}
\arguments{
\item{x}{A tramME object}

\item{newdata}{an optional data frame of observations}

\item{ranef}{Vector of random effects}

\item{...}{Additional arguments, passed to \code{\link[mlt]{plot.mlt}}.}
}
\value{
A numeric matrix of the predicted values invisibly
}
\description{
Plot the conditional distribution evaluated at a grid of possible response
values and a set of covariate and random effects values on a specified scale.
}
\examples{
data("sleepstudy", package = "lme4")
fit <- BoxCoxME(Reaction ~ Days + (Days | Subject), data = sleepstudy)
plot(fit, K = 100, type = "density")
}