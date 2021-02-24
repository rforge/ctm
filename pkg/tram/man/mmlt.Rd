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
mmlt(..., formula = ~ 1, data, theta = NULL, diag = FALSE,
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
  \item{diag}{logical; by default, diagonal elements are constrained to 1, corresponding to \code{diag = FALSE}.
}
}
\details{
  The function implements multivariate conditional transformation models
  as described by Klein et al (2020). The response is assumed absolutely
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

  Nadja Klein, Torsten Hothorn, Luisa Barbanti, Thomas Kneib (2020),
  Multivariate Conditional Transformation Models.
  <arxiv:1906.03151>

}
\examples{

  data("cars")

  ### fit unconditional bivariate distribution of speed and distance to stop
  ## fit unconditional marginal transformation models
  m_speed <- BoxCox(speed ~ 1, data = cars, support = ss <- c(4, 25), 
                    add = c(-5, 5))
  m_dist <- BoxCox(dist ~ 1, data = cars, support = sd <- c(0, 120), 
                   add = c(-5, 5))

  ## fit multivariate unconditional transformation model
  m_speed_dist <- mmlt(m_speed, m_dist, formula = ~ 1, data = cars)

  ## lambda defining the Cholesky of the precision matrix,
  ## with standard error
  coef(m_speed_dist, newdata = cars[1,], type = "Lambda")
  sqrt(vcov(m_speed_dist)["dist.sped.(Intercept)", 
                          "dist.sped.(Intercept)"])

  ## linear correlation, ie Pearson correlation of speed and dist after
  ## transformation to bivariate normality
  (r <- coef(m_speed_dist, newdata = cars[1,], type = "Corr"))
  
  ## Spearman's rho (rank correlation), can be computed easily 
  ## for Gaussian copula as
  (rs <- 6 * asin(r / 2) / pi)

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

  ## compute marginal densities
  nd_s <- as.data.frame(nd_s)
  nd_s$d <- predict(m_speed_dist, newdata = nd_s, type = "density")
  nd_d <- as.data.frame(nd_d)
  nd_d$d <- predict(m_speed_dist, newdata = nd_d, marginal = 2L, 
                    type = "density")

  ## plot bivariate and marginal distribution
  col1 <- rgb(.1, .1, .1, .9)
  col2 <- rgb(.1, .1, .1, .5)
  w <- c(.8, .2)
  layout(matrix(c(2, 1, 4, 3), nrow = 2), width = w, height = rev(w))
  par(mai = c(1, 1, 0, 0) * par("mai"))
  sp <- unique(nd$speed)
  di <- unique(nd$dist)
  d <- matrix(nd$d, nrow = length(sp))
  contour(sp, di, d, xlab = "Speed (in mph)", ylab = "Distance (in ft)", xlim = ss, ylim = sd,
          col = col1)
  points(cars$speed, cars$dist, pch = 19, col = col2)
  mai <- par("mai")
  par(mai = c(0, 1, 0, 1) * mai)
  plot(d ~ speed, data = nd_s, xlim = ss, type = "n", axes = FALSE, 
       xlab = "", ylab = "")
  polygon(nd_s$speed, nd_s$d, col = col2, border = FALSE)
  par(mai = c(1, 0, 1, 0) * mai)
  plot(dist ~ d, data = nd_d, ylim = sd, type = "n", axes = FALSE, 
       xlab = "", ylab = "")
  polygon(nd_d$d, nd_d$dist, col = col2, border = FALSE)

  ### NOTE: marginal densities are NOT normal, nor is the joint
  ### distribution. The non-normal shape comes from the data-driven 
  ### transformation of both variables to joint normality in this model.

}
\keyword{models}
