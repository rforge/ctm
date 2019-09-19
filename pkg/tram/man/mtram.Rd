\name{mtram}
\alias{mtram}
\title{
  Transformation Models for Clustered Data
}
\description{
  Marginally interpretable transformation models for clustered data.
  Highly experimental, use at your own risk.
}
\usage{
mtram(object, formula, data, standardise = FALSE, 
      grd = SparseGrid::createSparseGrid(type = "KPU", 
                dimension = length(rt$cnms[[1]]), k = 10), 
      Hessian = FALSE, ...)
}
\arguments{
  \item{object}{
  A \code{tram} object.
}
  \item{formula}{
  A formula specifying the random effects.
}
  \item{data}{
  A data frame. 
}
  \item{standardise}{
  Two types of models can be estimated: M1 (with \code{standardise = FALSE})
  corresponds to a marginal distribution function without direct
  interpretation of the fixed effects, M2 (with \code{standardise = TRUE})
  allows a marginal interpretation of scaled fixed effects as 
  log-odds or log-hazard ratios (depending on \code{object}). See Hothorn
  (2019).
}
  \item{grd}{
  A sparse grid used for numerical integration to get the likelihood.
}
  \item{Hessian}{
  A logical, if \code{TRUE}, the hessian is computed and returned.
}
  \item{\dots}{
  Additional argument.
}
}
\details{
  A Gaussian copula with a correlation structure obtained from a random
  intercept or random intercept / random slope model (that is, clustered or
  longitudinal data can by modelled only) is used to capture the
  correlations whereas the marginal distributions are described by a 
  transformation model. The methodology is described in Hothorn (2019) 
  and examples are given in the \code{mtram} package vignette.

  This is a proof-of-concept implementation and still highly experimental.
  Only \code{coef()} and \code{logLik()} methods are available at the
moment.
}
\value{
  An object of class \code{tram} with \code{coef()} and \code{logLik()}
  methods.
}
\references{

  Torsten Hothorn (2019). Marginally Interpretable Parametric Linear
  Transformation Models for Clustered Observations. Technical Report.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
