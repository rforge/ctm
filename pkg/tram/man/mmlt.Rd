\name{mmlt}
\alias{mmlt}
\title{
  Multivariate Conditional Transformation Models
}
\description{
  A proof-of-concept implementation of multivariate conditional
transformation models
}
\usage{
mmlt(..., formula = ~ 1, data, theta = NULL, 
     control.outer = list(trace = FALSE), scale = FALSE)
}
\arguments{
  \item{\dots}{marginal transformation models, one for each response
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
  \item{scale}{logical; parameters are not scaled prior to optimisation by default
}
}
\details{
  The function implements multivariate conditional transformation models
  as described by Klein et al (2019). The response is assumed absolutely
  continuous at the moment, discrete versions will be added later.

  Below is a simple example for an unconditional bivariate distribution.
  See \code{demo("undernutrition", package = "tram")} for a conditional
  three-variate example.
}
\value{
An object of class \code{mmlt} with \code{coef} and \code{predict}
methods.
}
\references{

  Nadja Klein, Torsten Hothorn, Thomas Kneib (2019),
  Multivariate Conditional Transformation Models.
  <arxiv:1906.03151>

}
\examples{

  data("cars")

  ### fit unconditional bivariate distribution of speed and distance to stop
  ## fit unconditional marginal transformation models
  m_speed <- BoxCox(speed ~ 1, data = cars, support = c(4, 25))
  m_dist <- BoxCox(dist ~ 1, data = cars, support = c(0, 120))

  ## fit multivariate unconditional transformation model
  m_speed_dist <- mmlt(m_speed, m_dist, formula = ~ 1, data = cars)

  ## "latent" correlation, ie correlation of speed and dist after
  ## transformation to bivariate normality
  coef(m_speed_dist, newdata = cars[1,], type = "Corr")

  ## evaluate joint and marginal densities (needs to be more user-friendly)
  nd <- expand.grid(c(nd_s <- mkgrid(m_speed, 100), nd_d <- mkgrid(m_dist, 100)))
  nd$hs <- predict(m_speed_dist, newdata = nd, marginal = 1L)
  nd$hps <- predict(m_speed_dist, newdata = nd, marginal = 1L, 
                    deriv = c("speed" = 1))
  nd$hd <- predict(m_speed_dist, newdata = nd, marginal = 2L)
  nd$hpd <- predict(m_speed_dist, newdata = nd, marginal = 2L, 
                    deriv = c("dist" = 1))

  ## joint density
  nd$d <- with(nd, 
               dnorm(hs) * 
               dnorm(coef(m_speed_dist)["dist.sped.(Intercept)"] * hs + hd) * 
               hps * hpd)

  ## plot joint density
  library("lattice")
  p1 <- contourplot(d ~ speed + dist, data = nd, 
      panel = function(...) {
          panel.contourplot(...)
          panel.xyplot(cars$speed, cars$dist, pch = 19)
  }, main = "Joint density")

  ## plot marginal densities
  nd_s <- as.data.frame(nd_s)
  nd_s$d <- predict(m_speed_dist, newdata = nd_s, type = "density")
  p2 <- xyplot(d ~ speed, data = nd_s, type = "l", 
               ylab = "Marginal density")

  nd_d <- as.data.frame(nd_d)
  nd_d$d <- predict(m_speed_dist, newdata = nd_d, marginal = 2L, 
                    type = "density")
  p3 <- xyplot(dist ~ d, data = nd_d, type = "l", xlab = "Marginal density")

  ### do the plot (not quite perfect yet)
  library("gridExtra")
  laymat <- matrix(c(2, 1, 1, 2, 1, 1, 4, 3, 3), nrow = 3)
  plot(arrangeGrob(p1, p2, p3, layout_matrix = laymat))

  ### NOTE: marginal densities are NOT normal, nor is the joint
  ### distribution. The non-normal shape comes from the data-driven 
  ### transformation of both variables to joint normality in this model.

}
\keyword{models}
