\name{mcotram}
\alias{mcotram}
\title{
  Multivariate Count Conditional Transformation Models
}
\description{
  A proof-of-concept implementation of multivariate conditional
transformation models for count data
}
\usage{
mcotram(..., formula = ~ 1, data, theta = NULL,
        control.outer = list(trace = FALSE), tol = sqrt(.Machine$double.eps), 
        dofit = TRUE)
}
\arguments{
  \item{\dots}{marginal count transformation models, one for each response
}
  \item{formula}{a model formula describing a model for the dependency
structure via the lambda parameters. The default is set to \code{~ 1} for constant lambdas.
}
  \item{data}{a data.frame
}
  \item{theta}{an optional vector of starting values
}
  \item{control.outer}{a list controlling \code{\link[alabama]{auglag}}
}
  \item{tol}{tolerance}
  \item{dofit}{logical; parameters are fitted by default, otherwise a list
               with log-likelihood and score function is returned
}
}
\details{
  The function implements multivariate count conditional transformation models.
  The response is assumed to be a vector of counts.
}
\value{
An object of class \code{mcotram} and \code{mmlt} with \code{coef} and 
\code{predict} methods.
}
\references{

  Nadja Klein, Torsten Hothorn, Luisa Barbanti, Thomas Kneib (2020),
  Multivariate Conditional Transformation Models. \emph{Scandinavian Journal
  of Statistics}, \doi{10.1111/sjos.12501}.

}
\examples{
  data("birds", package = "TH.data")
  ## fit conditional marginal count transformation models
  m_SG4 <- cotram(SG4 ~ AOT + AFS + GST + DBH + DWC + LOG, 
                  data = birds, method = "probit", 
                  add = c(-0.5, 10))
  m_SG5 <- cotram(SG5 ~ AOT + AFS + GST + DBH + DWC + LOG, 
                  data = birds, method = "probit", 
                  add = c(-0.5, 10))

  ## fit multivariate count conditional transformation model
  ## discrete approximation
  m_SG4_SG5 <- mcotram(m_SG4, m_SG5, formula = ~ 1, data = birds)
  m_SG5_SG4 <- mcotram(m_SG5, m_SG4, formula = ~ 1, data = birds)

  ## lambda defining the Cholesky of the precision matrix,
  ## with standard error
  coef(m_SG4_SG5, newdata = birds[1,], type = "Lambda")
  V <- vcov(m_SG4_SG5)["SG5.SG4.(Intercept)", "SG5.SG4.(Intercept)"]
  (se <- sqrt(V))

  ## linear correlation, ie Pearson correlation of SG4 and SG5 after
  ## transformation to bivariate normality
  (r <- coef(m_SG4_SG5, type = "Corr")[1])

  ## Spearman's rho (rank correlation), can be computed easily 
  ## for Gaussian copula as
  (rs <- 6 * asin(r / 2) / pi)
}
\keyword{models}
