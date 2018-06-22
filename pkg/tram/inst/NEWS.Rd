
\name{NEWS}
\title{NEWS file for the tram package}

\section{Changes in Version 0.2-2 (2018-06-22)}{
  \itemize{
    \item{Use default \code{bounds = c(0, Inf)} for \code{Surv} objects (for example
          in \code{Coxph}.}
    \item{Update to \pkg{mlt} 1.0-2.}
  }
}

\section{Changes in Version 0.2-1 (2018-04-20)}{
  \itemize{
    \item{Likelihood ratio tests with fixed parameters gave an error.}
    \item{Remove strata from linear predictor (experimental).}
    \item{Try to implement a better guess for the \code{add} argument
          later defining the axes corresponding to the response
          in plots.}
    \item{All models accept \code{add} and \code{bounds} arguments.}
  }
}

\section{Changes in Version 0.2-0 (2018-02-25)}{
  \itemize{
    \item{\pkg{tram} first published on CRAN.}
  }
}
