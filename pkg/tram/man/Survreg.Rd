\name{Survreg}
\alias{Survreg}
\title{
Parametric Survival Models
}
\description{
Weibull, log-normal, log-logistic and other parametric models for survival
analysis
}
\usage{
Survreg(formula, data, subset, weights, offset, cluster, na.action = na.omit, 
        dist = c("weibull", "logistic", "gaussian", "exponential", "rayleigh", 
                 "loggaussian", "lognormal", "loglogistic"), scale = 0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{an object of class \code{"formula"}: a symbolic description 
          of the model structure to be
          fitted.  The details of model specification are given under
          ‘Details’.
  }
  \item{data}{an optional data frame, list or environment (or object
          coercible by \code{as.data.frame} to a data frame) containing the
          variables in the model.  If not found in \code{data}, the
          variables are taken from \code{environment(formula)}.
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
          contain \code{NA}s.  The default is set by the \code{na.action} setting
          of \code{options}, and is \code{na.fail} if that is unset.
  }
  \item{dist}{character defining the conditional distribution of the
    (not necessarily positive) response, current choices include 
    Weibull, logistic, normal, exponential, Rayleigh, log-normal (same as
    log-gaussian), or log-logistic.
  }
  \item{scale}{a fixed value for the scale parameter(s).
  }
  \item{\dots}{additional arguments to \code{\link[tram]{tram}}.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (formula, data, subset, weights, offset, cluster, na.action = na.omit, 
    dist = c("weibull", "logistic", "gaussian", "exponential", 
        "rayleigh", "loggaussian", "lognormal", "loglogistic"), 
    scale = 0, ...) 
{
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action", "weights", 
        "offset", "cluster"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- quote(tram_data)
    td <- eval(mf, parent.frame())
    stopifnot(inherits(td$response, "Surv") || inherits(td$response, 
        "response") || is.numeric(td$response))
    dist <- match.arg(dist)
    distribution <- switch(dist, logistic = "Logistic", loglogistic = "Logistic", 
        gaussian = "Normal", loggaussian = "Normal", lognormal = "Normal", 
        weibull = "MinExtrVal", exponential = "MinExtrVal", rayleigh = "MinExtrVal")
    transformation <- switch(dist, logistic = "linear", loglogistic = "logarithmic", 
        gaussian = "linear", loggaussian = "logarithmic", lognormal = "logarithmic", 
        weibull = "logarithmic", exponential = "logarithmic", 
        rayleigh = "logarithmic")
    ret <- tram(td, transformation = transformation, distribution = distribution, 
        negative = TRUE, model_only = TRUE, ...)
    rname <- names(td$mf)[1]
    cfnm <- names(coef(ret))
    scalecf <- grep(rname, cfnm, fixed = TRUE)
    if (dist == "exponential") 
        scale <- 1
    if (dist == "rayleigh") 
        scale <- 0.5
    if (scale > 0) {
        fixed <- rep(1/scale, length(scalecf))
        names(fixed) <- cfnm[scalecf]
        ret <- tram(td, transformation = transformation, distribution = distribution, 
            negative = TRUE, fixed = fixed, ...)
    }
    else {
        ret <- tram(td, transformation = transformation, distribution = distribution, 
            negative = TRUE, ...)
    }
    if (!inherits(ret, "mlt")) 
        return(ret)
    if (length(scalecf) > 0) 
        ret$invscale <- 1/coef(as.mlt(ret))[scalecf]
    ret$call <- match.call(expand.dots = TRUE)
    .simpleCap <- function(x) {
        s <- strsplit(x, " ")[[1]]
        paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "", 
            collapse = " ")
    }
    ret$tram <- paste(ifelse(is.null(td$terms$s), "", "(Stratified)"), 
        .simpleCap(dist), "Linear Regression Model")
    class(ret) <- c("Survreg", class(ret))
    ret
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
