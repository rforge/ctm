\name{ctm-package}
\alias{ctm-package}
\docType{package}
\title{
  ctm: Conditional Transformation Models
}
\description{
  An experimental implementation of conditional transformation
  models for the semiparametric estimation of conditional distribution
  functions. Contains example analyses and a simulation study.
}
\details{
\tabular{ll}{
Package: \tab ctm\cr
Type: \tab Package\cr
Version: \tab 0.0-3\cr
Date: \tab 2012-09-04\cr
License: \tab GPL-2\cr
}

  This package is highly experimental! DO NOT USE this package
  for anything beyond demonstration and theoretical experiments.

  The package contains a very rough implementation of conditional 
  transformation models as introduced by Hothorn, Kneib, and 
  Buehlmann (2012). This class of semiparametric regression models
  fits conditional transformation functions that can be used
  to evaluate the conditional distribution of the response
  given explanatory variables.

  The applications presented in Hothorn, Kneib, and 
  Buehlmann (2013) can be reproduced by running `make' 
  in the \code{inst/applications} directory.

  The empirical evaluation reported on in the paper
  can be reproduced by running R on the file
  \code{sim.R} in the directory \code{inst/empeval}.

}
\author{
  Torsten Hothorn \email{Torsten.Hothorn@R-project.org}
}

\references{
  Torsten Hothorn, Thomas Kneib, and Peter Buehlmann (2013),
  Conditional Transformation Models. \emph{Journal of the Royal Statistical 
  Society: Series B (Statistical Methodology)}, accepted 2012-11-27. Preprint
  available from \url{http://arxiv.org/abs/1201.5786v1}.
}

\keyword{ package }
\keyword{ smooth }
\keyword{ nonparametric }
\keyword{ models }
