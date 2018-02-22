\name{BoxCox}
\alias{BoxCox}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  (Almost) Box-Cox Models
}
\description{
  Non-normal linear regression inspired by Box-Cox models
}
\usage{
BoxCox(formula, data, subset, weights, offset, cluster, na.action = na.omit, ...)
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
          contain \code{NA}s.  The default is set to \code{na.omit}.
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
    ...) 
{
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action", "weights", 
        "offset", "cluster"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- quote(tram_data)
    td <- eval(mf, parent.frame())
    stopifnot(inherits(td$response, "Surv") || inherits(td$response, 
        "response") || is.numeric(td$response))
    ret <- tram(td, transformation = "smooth", distribution = "Normal", 
        negative = TRUE, ...)
    if (!inherits(ret, "mlt")) 
        return(ret)
    ret$call <- match.call(expand.dots = TRUE)
    ret$tram <- paste(ifelse(is.null(td$terms$s), "", "(Stratified)"), 
        "Transformed Normal Linear Regression Model")
    class(ret) <- c("BoxCox", class(ret))
    ret
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
