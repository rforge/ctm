
\name{plot_path}
\alias{plot_path}
\title{Plot regularization paths for \code{"prof_*"} classes}
\usage{
plot_path(object, plot_logLik = FALSE, ...)
}
\arguments{
\item{object}{object of class \code{"prof_alpha"} or \code{"prof_lambda"}}

\item{plot_logLik}{Whether logLik trajectory should be plotted (default: \code{FALSE})}

\item{...}{additional arguments to \code{\link{plot}}}
}
\description{
Plot regularization paths and optionally log-likelihood trajectories of objects 
of class \code{"prof_alpha"} and \code{"prof_lambda"}. Coefficient names are
automatically added to the plot.
}
\examples{
\donttest{
library("tramnet")
library("survival")

data("nki70", package = "penalized")
nki70$resp <- with(nki70, Surv(time, event))
x <- scale(model.matrix( ~ 0 + DIAPH3 + NUSAP1 + TSPYL5 + C20orf46, data = nki70))
y <- Coxph(resp ~ 1, data = nki70, order = 10, log_first = TRUE)
fit1 <- tramnet(y, x, lambda = 0, alpha = 1)
pfl <- prof_lambda(fit1)
plot_path(pfl)
fit2 <- tramnet(y, x, lambda = 1, alpha = 1)
pfa <- prof_alpha(fit2)
plot_path(pfa)
}
}
\author{
Lucas Kook
}
\value{None}
