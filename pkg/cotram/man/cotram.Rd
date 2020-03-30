\name{cotram}
\alias{cotram}
\title{
	Count Transformation Models
}
\description{
	Likelihood-based count transformation models for fully parameterised discrete conditional
	distribution functions. The link function governing the interpretation of the predictor 
	can be chosen and results in discrete hazard ratios, odds ratios, reverse time
  hazard ratios or conditional expectation of transformed counts.
}
\usage{
cotram(formula, data, method = c("logit", "cloglog", "loglog", "probit"),
       log_first = TRUE, plus_one = log_first, prob = 0.9,	
       subset, weights, offset, cluster, na.action = na.omit, ...)
}
\arguments{
  \item{formula}{an object of class \code{"formula"}: a symbolic description 
          of the model structure to be
          fitted.  The details of model specification are given under
          \code{\link{tram}} and in the package vignette.
  }
  \item{data}{an optional data frame, list or environment (or object
          coercible by \code{as.data.frame} to a data frame) containing the
          variables in the model.  If not found in \code{data}, the
          variables are taken from \code{environment(formula)}.
  }
  \item{method}{character specifying the choice of the link function,
  			mapping the transformation function into probabilities. Available choices
  			include the logit, complementary log-log, log-log or probit link. 
  			The different link functions govern the interpretation of the linear
  			predictor. Details of the interpretation can be found in the 
  			package vignette.
  }
  \item{prob}{probability giving the quantile of the response defining the upper limit of the 
  	support of a smooth Bernstein polynomial.
  }
  \item{log_first}{logical; if \code{TRUE}, a Bernstein polynomial is
  			defined on the log-scale.
  }
  \item{plus_one}{logical; if \code{TRUE}, a Bernstein polynomial of (y + 1)
  			is defined.
  }
  \item{subset}{an optional vector specifying a subset of observations to be
          used in the fitting process.
  }
  \item{weights}{an optional vector of weights to be used in the fitting
          process.  Should be \code{NULL} or a numeric vector. If present,
          the weighted log-likelihood is maximised.
  }
  \item{offset}{this can be used to specify an _a priori_ known component to
          be included in the linear predictor during fitting.  This
          should be \code{NULL} or a numeric vector of length equal to the
          number of cases.
  }
  \item{cluster}{optional factor with a cluster ID employed for computing
          clustered covariances.
  }
  \item{na.action}{a function which indicates what should happen when the data
          contain \code{NA}s.  The default is set to \code{na.omit}.
  }
  \item{\dots}{additional arguments to \code{\link{tram}}.}
}
\details{
	Likelihood-based estimation of a fully parameterised conditional discrete
	distribution function for count data, while ensuring interpretability of
	the linear predictors. The models are defined with a negative shift term 
	relating positive predictors to larger values of the conditional mean. 
	For the model with logistic or cloglog link \code{exp(-coef())}
	is the multiplicative change of discrete odds-ratios or hazard ratios. For 
	the model with loglog link \code{exp(coef())} is the multiplicative change of
	the reverse time hazard ratios. Applying a transformation model with probit link
  \code{coef()} gives the conditional expectation of the transformed link,
  with transformation function estimated from data.
}
\value{
An object of class \code{cotram} and \code{tram}, with corresponding \code{coef},
\code{vcov}, \code{logLik}, \code{summary}, 
\code{print}, \code{plot} and \code{predict} methods.
}
\references{

  Sandra Siegfried, Torsten Hothorn (2020), Count Transformation Models,
  \emph{Methods in Ecology and Evolution}, \doi{10.1111/2041-210X.13383}.
  
  Torsten Hothorn, Lisa Moest, Peter Buehlmann (2018), Most Likely
  Transformations, \emph{Scandinavian Journal of Statistics}, \bold{45}(1),
  110--134, \doi{10.1111/sjos.12291}.
  
  Torsten Hothorn (2020), Most Likely Transformations: The mlt Package, 
  \emph{Journal of Statistical Software}, \bold{92}(1), 1--68,
  \doi{10.18637/jss.v092.i01}.

}
\examples{

  data("birds", package = "TH.data")
  cotram(SG5 ~ AOT + AFS + GST + DBH + DWC + LOG, data = birds)

}
\keyword{models}
\keyword{regression}
