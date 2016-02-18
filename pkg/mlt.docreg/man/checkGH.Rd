\name{checkGH}
\alias{checkGH}
\title{
  Check Gradient and Hessian
}
\description{
  Check analytic gradient and hessian against their numerical approximation
}
\usage{
checkGH(object)
}
\arguments{
  \item{object}{an object of class \code{mlt} as returned by \code{\link[mlt]{mlt}}}
}
\details{
  Computes the gradient and hessian analytically and numerically and compares the
  two. Warnings are suppressed for the numerical approximations as they are 
  not aware of model constraints which might get violated.
}
