
\name{prof_lambda}
\alias{prof_lambda}
\title{Profiling tuning parameters}
\usage{
prof_lambda(model, min_lambda = 0, max_lambda = 15, nprof = 5,
  as.lm = FALSE)
}
\arguments{
\item{model}{model of class \code{"tramnet"}}

\item{min_lambda}{minimal value of lambda (default = \code{0})}

\item{max_lambda}{maximal value of lambda (default = \code{15})}

\item{nprof}{number of profiling steps (default = \code{5})}

\item{as.lm}{return scaled coefficients for class \code{"tramnet_Lm"}}
}
\description{
Computes the regularization path of all coefficients for a single tuning
parameter, lambda, over a sequence of values.
}
\author{
Lucas Kook
}
\value{Object of class \code{"prof_lambda"} which contains the regularization
path of all coefficients and the log-likelihood over the penalty parameter
lambda}
\examples{
\donttest{
library("tramnet")
library("survival")

data("nki70", package = "penalized")
nki70$resp <- with(nki70, Surv(time, event))
x <- scale(model.matrix( ~ 0 + DIAPH3 + NUSAP1 + TSPYL5 + C20orf46, data = nki70))
y <- Coxph(resp ~ 1, data = nki70, order = 10, log_first = TRUE)
fit <- tramnet(y, x, lambda = 0, alpha = 1)
pfl <- prof_lambda(fit)
plot_path(pfl)
}
}