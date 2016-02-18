\name{as.basis}
\alias{as.basis}
\alias{as.basis.formula}
\alias{as.basis.factor_var}
\alias{as.basis.ordered_var}
\alias{as.basis.factor}
\alias{as.basis.ordered}
\title{
    Convert Formula or Factor to Basis Function
}
\description{
    Convert a formula or factor to basis functions
}
\usage{
as.basis(object, ...)
\method{as.basis}{formula}(object, data = NULL, remove_intercept = FALSE, 
         ui = NULL, ci = NULL, negative = FALSE, scale = FALSE, ...)
\method{as.basis}{factor_var}(object, ...)
\method{as.basis}{ordered_var}(object, ...)
\method{as.basis}{factor}(object, ...)
\method{as.basis}{ordered}(object, ...)
}
\arguments{
  \item{object}{a formula or an object of class \code{factor}, \code{factor_var},
                \code{ordered} or \code{ordered_var}}
  \item{data}{either a \code{\link[variables]{vars}} object or a \code{data.frame}}
  \item{remove_intercept}{a logical indicating if any intercept term shall be removed}
  \item{ui}{a matrix defining constraints}
  \item{ci}{a vector defining constraints}
  \item{negative}{a logical indicating negative basis functions}
  \item{scale}{a logical indicating a scaling of each column of the
               model matrix to the unit interval (based on observations in \code{data})}
  \item{\dots}{additional arguments to \code{\link{model.matrix}}, for
               example contrasts}
}
\details{
  \code{as.basis} returns a function for the evaluation of
  the basis functions with corresponding \code{model.matrix} and \code{predict}
  methods.

  Unordered factors (classes \code{factor} and \code{factor_var}) use a dummy coding and
  ordered factor (classes \code{ordered} or \code{ordered_var}) lead to a treatment contrast
  to the last level and removal of the intercept term with monotonicity constraint.
  Additional arguments (\code{\dots}) are ignored for ordered factors.

  Linear constraints on parameters \code{parm} are defined by \code{ui \%*\% parm >= ci}.
}
\examples{

  ## define variables and basis functions
  v <- c(numeric_var("x"), factor_var("y", levels = LETTERS[1:3]))
  fb <- as.basis(~ x + y, data = v, remove_intercept = TRUE, negative = TRUE,
                 contrasts.arg = list(y = "contr.sum"))

  ## evaluate basis functions
  model.matrix(fb, data = as.data.frame(v, n = 10))
  ## basically the same as (but wo intercept and times -1)
  model.matrix(~ x + y, data = as.data.frame(v, n = 10))

  ### factor
  xf <- gl(3, 1)
  model.matrix(as.basis(xf), data = data.frame(xf = xf))

  ### ordered
  xf <- gl(3, 1, ordered = TRUE)
  model.matrix(as.basis(xf), data = data.frame(xf = unique(xf)))

}
