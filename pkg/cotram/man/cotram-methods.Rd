\name{cotram-methods}
\alias{predict.cotram}
\alias{plot.cotram}
\alias{as.mlt.cotram}
\alias{logLik.cotram}
\title{
  Methods for Count Transformation Models
}
\description{
  Methods for objects inheriting from class cotram
}
\usage{
\method{predict}{cotram}(object, newdata = model.frame(object),
        type = c("lp", "trafo", "distribution", "survivor", "density", 
                 "logdensity", "hazard", "loghazard", "cumhazard", 
                 "logcumhazard", "odds", "logodds", "quantile"),
        smooth = FALSE, q = NULL, K = 20, prob = 1:(10-1)/10, ...)
\method{plot}{cotram}(x, newdata, type = c("distribution", "survivor","density",
                          "logdensity", "cumhazard", "quantile", "trafo"),
        confidence = c("none", "band"), level = 0.95, 
        smooth = FALSE, q = NULL, K = 20, cheat = K, prob = 1:(10-1)/10,
        col = "black", fill = "lightgrey",
        lty = 1, lwd = 1, add = FALSE, ...) 
\method{as.mlt}{cotram}(object)
\method{logLik}{cotram}(object, parm = coef(as.mlt(object), fixed = FALSE), newdata, ...)
}
\arguments{
  \item{object, x}{a fitted linear count transformation model inheriting
                from class \code{cotram}.}
  \item{newdata}{an optional data frame of observations.}
  \item{parm}{model parameters.}
  \item{type}{type of prediction, current options include
              linear predictors (\code{"lp"}, of \code{x} variables in the
              formula \code{y ~ x}), transformation functions
              (\code{"trafo"}) or distribution functions on the
              scale of the cdf (\code{"distribution"}), survivor function,
              density function, log-density function, cumulative
              hazard function or quantile function.}
  \item{confidence}{whether to plot a confidence band (see \code{\link[cotram]{confband}}).}
  \item{level}{the confidence level.}
  \item{smooth}{logical; if \code{TRUE} a smoothed function of \code{type} is returned.}
  \item{q}{quantiles at which to evaluate the model.}
  \item{prob}{probabilities for the evaluation of the quantile function \code{(type = "quantile")}.}
  \item{K}{number of grid points the function is evaluated at
          (for \code{smooth = TRUE} and in absence of \code{q}).}
  \item{cheat}{number of grid points the function is evaluated at when
               using the quantile obtained for \code{K} grid points
               (in the absence of \code{q} and \code{smooth = TRUE}).}
  \item{col}{color for the lines to plot.}
  \item{fill}{color for the confidence band.}
  \item{lty}{line type for the lines to plot.}
  \item{lwd}{line width.}
  \item{add}{logical; indicating if a new plot shall be generated (the default).}
  \item{\dots}{additional arguments to the underlying methods for \code{predict}.}
}
\details{
	\code{predict} and \code{plot} can be used to inspect the model on
  different scales.
}
\seealso{
  \code{\link[cotram]{predict.cotram}}, \code{\link[cotram]{confband.cotram}},
  \code{\link[tram]{tram-methods}}, \code{\link[mlt]{mlt-methods}}, \code{\link[mlt]{plot.ctm}}
}
\examples{
	
  data("birds", package = "TH.data")
  
  ### fit count transformation model with cloglog link
  m_birds <- cotram(SG5 ~ AOT + AFS + GST + DBH + DWC + LOG, data = birds,
                    method = "cloglog")
  logLik(m_birds)

  ### classical likelihood inference
  summary(m_birds)

  ### coefficients of the linear predictor (discrete hazard ratios)
  exp(-coef(m_birds))

  ### compute predicted median along with 10% and 90% quantile for the first
  ### three observations
  nd <- birds[1:3,]
  predict(m_birds, newdata = nd, type = "quantile", prob = c(.1, .5, .9), 
          smooth = TRUE)

  ### plot the predicted distribution for these observations
  plot(m_birds, newdata = nd, type = "distribution",
       col = c("skyblue", "grey", "seagreen"))
    
}
