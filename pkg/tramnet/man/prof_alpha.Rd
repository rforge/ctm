
\name{prof_alpha}
\alias{prof_alpha}
\title{Profiling tuning parameters}
\usage{
prof_alpha(model, min_alpha = 0, max_alpha = 1, nprof = 5,
  as.lm = FALSE)
}
\arguments{
\item{model}{model of class tramnet}

\item{min_alpha}{minimal value of alpha (default = \code{0})}

\item{max_alpha}{maximal value of alpha (default = \code{15})}

\item{nprof}{number of profiling steps (default = \code{5})}

\item{as.lm}{return scaled coefficients for class \code{"tramnet_Lm"}}
}
\description{
Computes the regularization path of all coefficients for a single tuning, alpha,
parameter over a sequence of values.
}
\author{
Lucas Kook
}
\value{Object of class \code{"prof_alpha"} which contains the regularization
path of all coefficients and the log-likelihood over the mixing parameter
alpha}
\examples{
\donttest{
library("tramnet")
library("survival")

data("nki70", package = "penalized")
nki70$resp <- with(nki70, Surv(time, event))
x <- scale(model.matrix( ~ 0 + DIAPH3 + NUSAP1 + TSPYL5 + C20orf46, data = nki70))
y <- Coxph(resp ~ 1, data = nki70, order = 10, log_first = TRUE)
fit <- tramnet(y, x, lambda = 1, alpha = 1)
pfa <- prof_alpha(fit)
plot_path(pfa)
}
}
