
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
Profiles a single tuning parameter for a tramnet model over
    a sequence of values and returns the trajectory of the coefficients
    and the in-sample logLik.
}
\author{
Lucas Kook
}
