\name{as.basis}
\alias{as.basis}
\alias{as.basis.formula}
\title{
    Formula to basis function
}
\description{
    Convert a formula to basis functions
}
\usage{
as.basis(object, ...)
\method{as.basis}{formula}(object, data = NULL, remove_intercept = FALSE, 
                ui = NULL, ci = NULL, negative = FALSE, ...)
}
\arguments{
  \item{object}{a formula}
  \item{data}{either a \code{\link[variables]{vars}} object or a \code{data.frame}}
  \item{remove_intercept}{a logical indicating if any intercept term shall be removed}
  \item{ui}{a matrix defining constraints}
  \item{ci}{a vector defining constraints}
  \item{negative}{a logical indicating negative basis functions}
  \item{\dots}{additional arguments to \code{\link{model.matrix}}, for
               example contrasts}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\examples{

  ## define variables and basis functions
  v <- c(numeric_var("x"), factor_var("y", levels = LETTERS[1:3]))
  fb <- as.basis(~ x + y, data = v, remove_intercept = TRUE, negative = TRUE,
                 contrasts.arg = list(y = "contr.sum"))

  ## evaluate basis functions
  model.matrix(fb, data = as.data.frame(v, n = 10)))
  ## basically the same as (but wo intercept and times -1)
  model.matrix(~ x + y, data = as.data.frame(v, n = 10))

}
