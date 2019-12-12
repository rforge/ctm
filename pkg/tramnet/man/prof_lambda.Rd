
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
Profiles a single tuning parameter for a \code{"tramnet"} model over
    a sequence of values and returns the trajectory of the coefficients
    and the in-sample logLik.
}
\author{
Lucas Kook
}
