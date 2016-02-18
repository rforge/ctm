\name{predict.basis}
\alias{predict.basis}
\alias{predict.cbind_bases}
\alias{predict.box_bases}
\title{
  Evaluate Basis Functions
}
\description{
  Evaluate basis functions and compute the function defined by the corresponding basis
}
\usage{
\method{predict}{basis}(object, newdata, coef, dim = !is.data.frame(newdata), ...)
\method{predict}{cbind_bases}(object, newdata, coef, dim = !is.data.frame(newdata), 
         terms = names(object), ...)
\method{predict}{box_bases}(object, newdata, coef, dim = !is.data.frame(newdata), ...)
}
\arguments{
  \item{object}{a \code{basis} or \code{bases} object}
  \item{newdata}{a \code{list} or \code{data.frame}}
  \item{coef}{a vector of coefficients}
  \item{dim}{either a logical indicating that the dimensions shall be
             obtained from the \code{bases} object or an integer vector
             with the corresponding dimensions (the latter option being very
             experimental} 
  \item{terms}{a character vector defining the elements of a \code{cbind_bases}
               object to be evaluated}
  \item{\dots}{additional arguments}
}
\details{
  \code{predict} evaluates the basis functions and multiplies them with \code{coef}.
  There is no need to expand multiple variables as \code{predict} uses array models
  (Currie et al, 2006) to compute the corresponding predictions efficiently.
}
\references{

  Ian D. Currie, Maria Durban, Paul H. C. Eilers, P. H. C. (2006), Generalized Linear 
  Array Models with Applications to Multidimensional Smoothing, 
  \emph{Journal of the Royal Statistical Society, Series B: Methodology}, \bold{68}(2), 
  259--280.

}
\examples{

  ### set-up a Bernstein polynomial
  xv <- numeric_var("x", support = c(1, pi))
  bb <- Bernstein_basis(xv, order = 3, ui = "increasing")
  ## and treatment contrasts for a factor at three levels
  fb <- as.basis(~ g, data = factor_var("g", levels = LETTERS[1:3]))

  ### join them: we get one intercept and two deviation _functions_
  bfb <- b(bern = bb, f = fb)

  ### generate data + coefficients
  x <- mkgrid(bfb, n = 10)
  cf <- c(1, 2, 2.5, 2.6)
  cf <- c(cf, cf + 1, cf + 2)

  ### evaluate predictions for all combinations in x (a list!)
  predict(bfb, newdata = x, coef = cf)

  ## same but slower
  matrix(predict(bfb, newdata = expand.grid(x), coef = cf), ncol = 3)

}
